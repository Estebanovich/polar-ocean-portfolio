"""Machine-specific locations, overridable with environment variables.

EXP1_RUNS     run-output root of Experiment 1 (default: the external run disk)
EXP1_BUILDS   MITgcm build root of Experiment 1
MITGCM_ROOT   MITgcm source checkout (read-only)
"""
from __future__ import annotations

import os
from pathlib import Path

EXP1_RUNS = Path(os.environ.get("EXP1_RUNS", "/Volumes/LaCie/portfolio_runs/exp1"))
EXP1_BUILDS = Path(os.environ.get("EXP1_BUILDS", "/Volumes/ESTEBAN_MIT/portfolio_builds/exp1"))
MITGCM_ROOT = Path(os.environ.get("MITGCM_ROOT", "/Volumes/ESTEBAN_MIT/MITgcm"))
