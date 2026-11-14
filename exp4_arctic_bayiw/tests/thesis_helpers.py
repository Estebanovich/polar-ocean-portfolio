"""Shared test helper: the frozen thesis wind_patch.py (thesis_copy/BayIW_Rectan/input)."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

WIND_PATCH = Path(__file__).resolve().parents[1] / "thesis_copy" / "BayIW_Rectan" / "input" / "wind_patch.py"


def thesis_wind_patch_module():
    """Import the verbatim thesis wind_patch.py (non-interactive matplotlib backend)."""
    os.environ.setdefault("MPLBACKEND", "Agg")
    spec = importlib.util.spec_from_file_location("thesis_wind_patch", WIND_PATCH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
