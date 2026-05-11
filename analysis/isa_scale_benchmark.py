"""Scale benchmark: nwt_substrate.isa across backends.

Compares chemistry aromaticity-RE compute across four paths:
  1. legacy networkx path (smiles_resonance_energy per molecule)
  2. ISA numpy backend (smiles_resonance_energy_batch backend=numpy)
  3. ISA torch_cpu backend
  4. ISA torch_cuda backend  (if CUDA available)

For each: time the END-TO-END pipeline AND the kernel-only step
(once parse + ring extraction is paid).

Reports speedup curves. Verifies all four paths produce identical
RE values (correctness).
"""
import time
import numpy as np

import nwt_substrate.chemistry as chem
import nwt_substrate.isa as isa


# Reference set with weights to scale up
BASE_SMILES = [
    "c1ccccc1",                                  # benzene
    "c1ccc2ccccc2c1",                            # naphthalene
    "c1ccc2cc3ccccc3cc2c1",                      # anthracene
    "c1cc2ccc3cccc4ccc(c1)c2c34",                # pyrene
    "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",   # coronene
    f"c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61"
        "-c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",  # dicoronylene
    "c1ccncc1",                                  # pyridine
    "c1cc[nH]c1",                                # pyrrole
    "c1ccoc1",                                   # furan
    "Cc1ccccc1",                                 # toluene
]


def make_batch(n: int) -> list[str]:
    """Cycle the base SMILES set up to size n."""
    out = []
    while len(out) < n:
        out.extend(BASE_SMILES)
    return out[:n]


def bench_legacy(smiles_list: list[str]) -> tuple[float, np.ndarray]:
    """Per-molecule legacy networkx path."""
    re = np.zeros(len(smiles_list))
    t0 = time.perf_counter()
    for i, smi in enumerate(smiles_list):
        re[i] = chem.smiles_resonance_energy(smi)
    t = time.perf_counter() - t0
    return t, re


def bench_isa(smiles_list: list[str], backend: str) -> tuple[float, float, np.ndarray]:
    """ISA path. Returns (total_time, kernel_only_time, re_array).

    kernel_only_time is the time for the substrate-native step
    (so(7) trace contraction + RE assembly), excluding SMILES parse.
    """
    # Step 1: parse + extract (per-molecule, can't be GPU-batched)
    t0 = time.perf_counter()
    sys_adj = []
    sys_pi = []
    sys_nu = []
    owner = []
    for k, smi in enumerate(smiles_list):
        result = chem.smiles_to_aromaticity(smi)
        if result.classification not in ("aromatic", "mobius_aromatic"):
            continue
        for s in result.aromatic_ring_systems:
            adj, n_used = chem.ring_system_to_so7_adjacency(s)
            sys_adj.append(adj)
            sys_pi.append(s.n_pi_pairs)
            sys_nu.append(n_used)
            owner.append(k)
    t_parse = time.perf_counter() - t0

    adj_batch = np.stack(sys_adj) if sys_adj else np.zeros((0, 7, 7))
    pi_batch = np.asarray(sys_pi)
    nu_batch = np.asarray(sys_nu)
    owner_arr = np.asarray(owner)

    # Step 2: substrate-native kernel
    # For torch_cuda: include a CUDA synchronize for honest timing
    if backend == "torch_cuda":
        import torch
        torch.cuda.synchronize()
    t0 = time.perf_counter()
    re_per_sys = isa.aromaticity_score(
        adj_batch, pi_batch,
        n_vertices_used=nu_batch,
        backend=backend,
    )
    if backend == "torch_cuda":
        import torch
        torch.cuda.synchronize()
    t_kernel = time.perf_counter() - t0

    re_per_mol = np.zeros(len(smiles_list))
    np.add.at(re_per_mol, owner_arr, re_per_sys)

    return t_parse + t_kernel, t_kernel, re_per_mol


def warmup_cuda():
    """Run a small batch on CUDA to warm up the device."""
    if "torch_cuda" in isa.available_backends():
        A = np.zeros((10, 7, 7))
        isa.graphs_to_invariants(A, backend="torch_cuda")
        import torch
        torch.cuda.synchronize()


def report_row(label: str, t_total: float, t_kernel: float | None,
                 n: int, ref_re: np.ndarray, this_re: np.ndarray,
                 t_legacy: float):
    correct = "✓" if np.allclose(ref_re, this_re) else "✗"
    spk_total = t_legacy / t_total
    if t_kernel is not None:
        kernel_str = f"{t_kernel*1000:>10.2f}"
        spk_kernel = t_legacy / t_kernel
        kernel_speedup = f"{spk_kernel:>9.1f}×"
    else:
        kernel_str = " " * 10 + "—"
        kernel_speedup = "        —"
    print(f"  {label:<18} {t_total*1000:>10.2f}  {kernel_str}  "
          f"{spk_total:>9.1f}×  {kernel_speedup}   {correct}")


def main():
    print("=" * 80)
    print("nwt_substrate.isa scale benchmark")
    print("=" * 80)
    print(f"Backends available: {isa.available_backends()}")

    warmup_cuda()

    for n in [100, 1000, 10_000]:
        print(f"\n--- N = {n:>6} molecules ---")
        smiles_list = make_batch(n)

        # Legacy
        t_legacy, re_legacy = bench_legacy(smiles_list)

        print(f"  {'path':<18} {'total (ms)':>10}  {'kernel (ms)':>10}  "
              f"{'total spk':>10}  {'kernel spk':>10}    OK")
        print("  " + "-" * 75)

        # Header row using legacy as 1×
        print(f"  {'legacy networkx':<18} {t_legacy*1000:>10.2f}  "
              f"{' '*10}—  {1.0:>9.1f}×  {' '*8}—   ✓")

        for backend in ["numpy", "torch_cpu", "torch_cuda"]:
            if backend not in isa.available_backends():
                continue
            t_total, t_kernel, re_isa = bench_isa(smiles_list, backend)
            label = f"isa.{backend}"
            report_row(label, t_total, t_kernel, n, re_legacy, re_isa, t_legacy)

    # Pure kernel benchmark (skip parse) at very large scale
    print(f"\n--- Pure substrate-kernel benchmark (no parse) ---")
    print("  Batched einsum on (N, 7, 7) — substrate ISA in isolation")
    print()
    for n in [1_000, 10_000, 100_000, 1_000_000]:
        # Build random {0,1} adjacency (symmetric) batch
        rng = np.random.default_rng(42)
        upper = rng.integers(0, 2, size=(n, 7, 7))
        A = np.triu(upper, k=1)
        A = A + A.transpose(0, 2, 1)

        timings = {}
        for backend in ["numpy", "torch_cpu", "torch_cuda"]:
            if backend not in isa.available_backends():
                continue
            if backend == "torch_cuda":
                import torch
                torch.cuda.synchronize()
            t0 = time.perf_counter()
            isa.graphs_to_invariants(A, backend=backend, return_numpy=False)
            if backend == "torch_cuda":
                import torch
                torch.cuda.synchronize()
            timings[backend] = time.perf_counter() - t0

        print(f"  N = {n:>8}", end="  ")
        for backend, t in timings.items():
            per_mol_ns = t / n * 1e9
            print(f"{backend}: {t*1000:>7.2f} ms ({per_mol_ns:>6.1f} ns/mol)", end="  ")
        print()


if __name__ == "__main__":
    main()
