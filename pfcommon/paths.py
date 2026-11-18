"""Machine-specific locations, overridable with environment variables.

EXP1_RUNS     run-output root of Experiment 1 (default: the external run disk)
EXP1_BUILDS   MITgcm build root of Experiment 1
EXP2_RUNS     run-output root of Experiment 2 (ISOMIP+ cavity)
EXP2_BUILDS   MITgcm build root of Experiment 2
EXP4_RUNS     run-output root of Experiment 4 (Arctic BayIW, Van Mijenfjorden)
EXP4_BUILDS   MITgcm build root of Experiment 4
EXP5_RUNS     run-output root of Experiment 5 (adjoint, Antarctic slope channel)
EXP5_BUILDS   MITgcm build root of Experiment 5
EXP6_RUNS     run-output root of Experiment 6 (Southern Ocean carbon channel)
EXP6_BUILDS   MITgcm build root of Experiment 6
TAPENADE_HOME Tapenade install used for the MITgcm adjoint (genmake2 -tap)
MITGCM_ROOT   MITgcm source checkout (read-only)
"""
from __future__ import annotations

import os
from pathlib import Path

EXP1_RUNS = Path(os.environ.get("EXP1_RUNS", "/Volumes/LaCie/portfolio_runs/exp1"))
EXP1_BUILDS = Path(os.environ.get("EXP1_BUILDS", "/Volumes/ESTEBAN_MIT/portfolio_builds/exp1"))
EXP2_RUNS = Path(os.environ.get("EXP2_RUNS", "/Volumes/LaCie/portfolio_runs/exp2"))
EXP2_BUILDS = Path(os.environ.get("EXP2_BUILDS", "/Volumes/ESTEBAN_MIT/portfolio_builds/exp2"))
EXP4_RUNS = Path(os.environ.get("EXP4_RUNS", "/Volumes/LaCie/portfolio_runs/exp4"))
EXP4_BUILDS = Path(os.environ.get("EXP4_BUILDS", "/Volumes/ESTEBAN_MIT/portfolio_builds/exp4"))
EXP5_RUNS = Path(os.environ.get("EXP5_RUNS", "/Volumes/LaCie/portfolio_runs/exp5"))
EXP5_BUILDS = Path(os.environ.get("EXP5_BUILDS", "/Volumes/ESTEBAN_MIT/portfolio_builds/exp5"))
EXP6_RUNS = Path(os.environ.get("EXP6_RUNS", "/Volumes/LaCie/portfolio_runs/exp6"))
EXP6_BUILDS = Path(os.environ.get("EXP6_BUILDS", "/Volumes/ESTEBAN_MIT/portfolio_builds/exp6"))
TAPENADE_HOME = Path(os.environ.get(
    "TAPENADE_HOME", "/Volumes/ESTEBAN_MIT/tools/tapenade-e89660a5921a1ccf68b4afc331b3b0910878a179"))
MITGCM_ROOT = Path(os.environ.get("MITGCM_ROOT", "/Volumes/ESTEBAN_MIT/MITgcm"))
