"""Tests for the nwt.heron quantum-circuit shim."""

import pytest

import nwt_substrate.heron as heron


# Skip everything if qiskit not available
pytestmark = pytest.mark.skipif(not heron.HAS_QISKIT,
                                  reason="qiskit not installed")


def test_k7_graph_state_structure():
    """K_7 prep: 7 H + 21 CZ = 28 gates on 7 qubits."""
    qc = heron.k7_graph_state()
    assert qc.num_qubits == 7
    info = heron.circuit_summary(qc)
    assert info["gate_counts"]["h"] == 7
    assert info["gate_counts"]["cz"] == 21
    assert info["n_gates"] == 28


def test_stabilizer_measurement_uses_8_qubits():
    """Stabilizer measurement needs an ancilla -> 8 qubits."""
    qc = heron.stabilizer_measurement(vertex=3)
    assert qc.num_qubits == 8
    assert qc.num_clbits == 1


def test_stabilizer_measurement_has_correct_structure():
    """For S_v we need: H on ancilla, CX(anc, v), 6 CZs (other vertices), H, measure."""
    qc = heron.stabilizer_measurement(vertex=2)
    info = heron.circuit_summary(qc)
    # 2 H on ancilla, 1 CX, 6 CZ, 1 measurement
    assert info["gate_counts"]["h"] == 2
    assert info["gate_counts"]["cx"] == 1
    assert info["gate_counts"]["cz"] == 6
    assert info["gate_counts"]["measure"] == 1


def test_full_self_test_circuit():
    """Full self-test: prep + 7 stabilizer measurements -> 7 classical bits."""
    qc = heron.full_k7_state_prep_with_measurement()
    assert qc.num_qubits == 8
    assert qc.num_clbits == 7


def test_x_basis_tomography():
    """X-basis tomography: prep, rotate, measure all 7 qubits."""
    qc = heron.entanglement_tomography_x_basis()
    assert qc.num_qubits == 7
    assert qc.num_clbits == 7
    info = heron.circuit_summary(qc)
    # 14 H (7 prep + 7 rotate to X) + 21 CZ + 7 measurement
    assert info["gate_counts"]["h"] == 14
    assert info["gate_counts"]["cz"] == 21
    assert info["gate_counts"]["measure"] == 7


# ----- Experiment registry -----

def test_experiment_registry_has_11_entries():
    assert len(heron.EXPERIMENTS) == 11


def test_summary_counts_run_and_proposed():
    s = heron.summary()
    assert s["total"] == 11
    assert s["by_status"]["run"] >= 1
    assert s["by_status"]["proposed"] >= 1


def test_lookup_experiment_11_sidereal():
    e11 = heron.experiment(11)
    assert e11.number == 11
    assert e11.status == "proposed"
    assert "sidereal" in e11.name.lower() or "sidereal" in e11.description.lower()


def test_lookup_experiment_by_number():
    e9 = heron.experiment(9)
    assert e9.number == 9
    assert "neutron" in e9.name.lower()


def test_filter_by_status():
    run = heron.list_experiments("run")
    proposed = heron.list_experiments("proposed")
    assert all(e.status == "run" for e in run)
    assert all(e.status == "proposed" for e in proposed)


def test_unknown_experiment_raises():
    with pytest.raises(KeyError):
        heron.experiment(999)


# ----- Script export -----

def test_export_single_experiment(tmp_path):
    """export_experiment_script writes a self-contained, compilable .py file."""
    out = tmp_path / "exp01.py"
    heron.export_experiment_script(1, out)
    assert out.exists()
    assert out.stat().st_size > 1000  # has substantive content
    code = out.read_text()
    # Required boilerplate
    assert "QISKIT_IBM_TOKEN" in code
    assert "qiskit_ibm_runtime" in code
    assert "ibm_marrakesh" in code
    # The K_7 prep should be inlined
    assert "qc.h(" in code
    assert "qc.cz(" in code
    # Compiles to valid Python
    import py_compile
    py_compile.compile(str(out), doraise=True)


def test_export_all_experiments(tmp_path):
    """export_all_experiments writes one file per registered experiment."""
    paths = heron.export_all_experiments(tmp_path)
    assert len(paths) == len(heron.EXPERIMENTS)
    import py_compile
    for p in paths:
        py_compile.compile(str(p), doraise=True)


def test_exported_script_contains_experiment_metadata(tmp_path):
    """Generated script's docstring carries experiment number, name, paper."""
    out = tmp_path / "exp09.py"
    heron.export_experiment_script(9, out)
    code = out.read_text()
    assert "Experiment 9" in code
    assert "neutron decay" in code.lower()
    assert "Paper 20" in code or "paper 20" in code.lower()
    assert "proposed" in code.lower()
