"""Make `python -m nwt_substrate.benchmarks` runnable.

Usage:
    python -m nwt_substrate.benchmarks                  # full output
    python -m nwt_substrate.benchmarks --summary        # just the summary line
    python -m nwt_substrate.benchmarks --json           # machine-readable JSON
"""

from __future__ import annotations
import argparse
import json
import sys

from . import run_all


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run nwt_substrate.benchmarks comparing substrate algebra "
                    "to traditional methods.")
    parser.add_argument("--summary", action="store_true",
                        help="Print only the summary (no per-benchmark table)")
    parser.add_argument("--json", action="store_true",
                        help="Output machine-readable JSON")
    parser.add_argument("--max-time-ms", type=float, default=None,
                        help="Regression check: fail if total time exceeds this (ms)")
    args = parser.parse_args(argv)

    results = run_all(verbose=not (args.summary or args.json))
    total_us = sum(r.substrate_time_us for r in results)

    if args.json:
        payload = {
            "n_benchmarks": len(results),
            "total_substrate_time_us": total_us,
            "total_substrate_time_ms": total_us / 1000,
            "results": [
                {
                    "name": r.name,
                    "substrate_time_us": r.substrate_time_us,
                    "substrate_value": r.substrate_value,
                    "substrate_accuracy": r.substrate_accuracy,
                    "traditional_method": r.traditional_method,
                    "traditional_cost": r.traditional_cost,
                    "speedup_factor": r.speedup_factor_str,
                    "notes": r.notes,
                }
                for r in results
            ],
        }
        print(json.dumps(payload, indent=2))
    elif args.summary:
        print(f"{len(results)} benchmarks in {total_us:.0f} μs "
              f"({total_us / 1000:.2f} ms)")

    if args.max_time_ms is not None and total_us / 1000 > args.max_time_ms:
        print(f"REGRESSION: total time {total_us / 1000:.2f} ms > "
              f"{args.max_time_ms:.2f} ms", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
