"""Melt-emulator method interface and registry.

``REGISTRY`` maps method ids ("M0".."M7") to zero-argument factories returning an
unfitted ``Method``. Importing this package registers every method whose module
is importable without torch: M1-M3 (physics) here, and later tasks add M0, M4-M7
by importing their modules at the bottom of this file. Torch-based methods (M6,
M7) must never be imported eagerly: their factories import their module only when
called, so the pyocean suite (no torch) can import this package.
"""
from __future__ import annotations

from typing import Callable, Protocol, runtime_checkable

import numpy as np

from exp3_melt_emulator.data import Run


@runtime_checkable
class Method(Protocol):
    """Melt emulator. ``predict`` returns a (ny, nx) map, NaN outside run.shelf_mask."""
    name: str

    def fit(self, runs: list[Run]) -> None: ...

    def predict(self, run: Run) -> np.ndarray: ...


@runtime_checkable
class ProbMethod(Method, Protocol):
    """Method with a predictive distribution: (mean, sigma) maps, each (ny, nx)."""

    def predict_dist(self, run: Run) -> tuple[np.ndarray, np.ndarray]: ...


REGISTRY: dict[str, Callable[[], Method]] = {}

from exp3_melt_emulator.methods import physics as _physics  # noqa: E402

REGISTRY["M1"] = _physics.LinearLocal
REGISTRY["M2"] = _physics.QuadLocal
REGISTRY["M3"] = _physics.QuadNonLocal


def _lazy(cls_name: str) -> Callable[[], Method]:
    """Factory importing methods.stats (needs scikit-learn) only when called."""
    def factory() -> Method:
        from exp3_melt_emulator.methods import stats
        return getattr(stats, cls_name)()
    factory.__name__ = f"make_{cls_name}"
    return factory


REGISTRY["M0"] = _lazy("MeanMap")
REGISTRY["M4"] = _lazy("CellRegression")
REGISTRY["M5"] = _lazy("GPonEOF")
