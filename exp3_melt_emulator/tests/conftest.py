"""Shared fixtures. No global importorskip: modules needing torch/sklearn skip themselves."""
import pytest

from exp3_melt_emulator.tests.synth import synthetic_runs


@pytest.fixture
def synth_runs():
    return synthetic_runs()
