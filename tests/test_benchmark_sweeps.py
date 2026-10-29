"""The benchmark emits reusable measurements for a small finite sweep."""

import json
import subprocess
import sys
from pathlib import Path


def test_sweep_benchmark_reports_measurements():
    script = Path(__file__).resolve().parents[1] / "tools/benchmark_sweeps.py"

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            "--shape",
            "3",
            "3",
            "3",
            "--wavelengths",
            "2",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    report = json.loads(completed.stdout)

    assert report["shape"] == [3, 3, 3]

    assert report["geometry_preparations"] == 1

    assert len(report["solve_seconds"]) == 2

    assert all(value > 0 for value in report["solve_seconds"])

    assert report["peak_traced_bytes"] > 0
