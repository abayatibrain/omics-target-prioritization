"""Shared fixtures for the test suite."""

from __future__ import annotations

import pytest

from omics_target_prioritization.simulate import SimConfig, SimulatedLocus, simulate_locus


@pytest.fixture
def sim_config() -> SimConfig:
    """A small, fast, deterministic simulation configuration."""
    return SimConfig(seed=7, n_snps=120, n_neighbour_genes=3)


@pytest.fixture
def sim_locus(sim_config: SimConfig) -> SimulatedLocus:
    """A simulated locus with a known causal gene."""
    return simulate_locus(sim_config)
