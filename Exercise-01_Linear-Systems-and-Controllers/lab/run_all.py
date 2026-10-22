"""Run every Exercise 1 lab experiment in order and regenerate all figures and GIFs.

Usage:  python run_all.py
"""
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
EXPERIMENTS = [
    "exp1_eigenvalues_and_modes.py",
    "exp2_drone_open_loop.py",
    "exp3_controllability.py",
    "exp4_observability.py",
    "exp5_state_feedback.py",
    "exp6_tracking_error.py",
    "exp7_prefilter_vs_integrator.py",
]
for script in EXPERIMENTS:
    runpy.run_path(str(HERE / script), run_name="__main__")
print("\nAll experiments finished. Figures are in", HERE / "figures")
