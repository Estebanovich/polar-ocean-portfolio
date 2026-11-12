"""Tests for the cache keys of exp4_arctic_bayiw.analysis.thesis_contrast (no model output)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from exp4_arctic_bayiw.analysis import thesis_contrast as tc


def _src(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "src.py"
    p.write_text(text)
    return p


def test_cache_key_depends_on_run_and_source_content(tmp_path):
    src = _src(tmp_path, "a = 1\n")
    k1 = tc.cache_file(tmp_path, "library_kf", "R1", [src])
    assert k1 == tc.cache_file(tmp_path, "library_kf", "R1", [src])
    assert k1 != tc.cache_file(tmp_path, "library_kf", "R1k", [src])     # another run
    src.write_text("a = 2\n")
    assert k1 != tc.cache_file(tmp_path, "library_kf", "R1", [src])      # changed source


def test_cached_misses_on_new_key_and_honours_recompute(tmp_path):
    calls = []

    def compute():
        calls.append(1)
        return {"x": np.arange(3.0)}

    src = _src(tmp_path, "a = 1\n")
    f = tc.cache_file(tmp_path, "library_kf", "R1", [src])
    tc.cached(f, compute)
    tc.cached(f, compute)
    assert len(calls) == 1                                   # hit
    tc.cached(f, compute, recompute=True)
    assert len(calls) == 2                                   # forced
    tc.cached(tc.cache_file(tmp_path, "library_kf", "R1k", [src]), compute)
    assert len(calls) == 3                                   # other run: miss
    src.write_text("a = 3\n")
    r = tc.cached(tc.cache_file(tmp_path, "library_kf", "R1", [src]), compute)
    assert len(calls) == 4 and np.array_equal(r["x"], np.arange(3.0))   # changed source: miss


def test_every_cache_is_keyed():
    keys = tc.CACHE_SOURCES
    assert set(keys) == {"thesis_nash", "thesis_kf", "thesis_diss", "library_kf"}
    assert all(len(v) > 0 for v in keys.values())
