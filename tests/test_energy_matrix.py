"""Tests for larest.energy_matrix energy descriptors."""

from __future__ import annotations

import numpy as np
import pytest

from larest.energy_matrix import (
    INF_LABEL,
    R2_COLUMN,
    compute_energy_matrices,
    compute_energy_matrix,
    extrapolate_to_infinite_length,
    format_energy_matrix,
    low_r2_warning,
)

# energies in J/mol
_MONOMER = -100_000.0
_INITIATOR = -40_000.0


class TestExtrapolate:
    def test_exact_linear_in_inverse_length(self):
        lengths = [2, 3, 4]
        values = [-10.0 + 6.0 / n for n in lengths]
        intercept, r2 = extrapolate_to_infinite_length(lengths, values)
        assert intercept == pytest.approx(-10.0)
        assert r2 == pytest.approx(1.0)

    def test_ignores_nan(self):
        intercept, _ = extrapolate_to_infinite_length(
            [2, 3, 4],
            [-7.0, np.nan, -8.5],
        )
        assert intercept == pytest.approx(-10.0)

    def test_single_point_not_extrapolated(self):
        value, r2 = extrapolate_to_infinite_length([2], [-5.0])
        assert np.isnan(value)
        assert r2 is None

    def test_no_points(self):
        value, r2 = extrapolate_to_infinite_length([], [])
        assert np.isnan(value)
        assert r2 is None


class TestComputeEnergyMatrixROR:
    def setup_method(self):
        self.polymers = {1: -150_000.0, 2: -255_000.0, 3: -358_000.0}
        self.matrix = compute_energy_matrix(
            _MONOMER,
            self.polymers,
            _INITIATOR,
            "ROR",
        )

    def test_diagonal_is_reaction_energy(self):
        for n, p in self.polymers.items():
            expected = (p - n * _MONOMER - _INITIATOR) / 1000
            assert self.matrix.loc[n, n] == pytest.approx(expected)

    def test_column_zero_is_normalised(self):
        for n in self.polymers:
            assert self.matrix.loc[n, 0] == pytest.approx(self.matrix.loc[n, n] / n)

    def test_incremental_energy(self):
        # IE^3_1 = (P3 - P1) / 2 - M
        expected = ((-358_000.0 + 150_000.0) / 2 - _MONOMER) / 1000
        assert self.matrix.loc[3, 1] == pytest.approx(expected)
        # IE^n_(n-1) is the energy of a single insertion
        assert self.matrix.loc[3, 2] == pytest.approx(
            (-358_000.0 + 255_000.0 - _MONOMER) / 1000,
        )

    def test_upper_triangle_is_empty(self):
        assert np.isnan(self.matrix.loc[1, 2])
        assert np.isnan(self.matrix.loc[2, 3])

    def test_extrapolated_row(self):
        assert not np.isnan(self.matrix.loc[INF_LABEL, 0])
        assert 0.0 <= self.matrix.loc[INF_LABEL, R2_COLUMN] <= 1.0


class TestComputeEnergyMatrixRER:
    def test_loop_of_length_one_is_monomer(self):
        polymers = {2: -205_000.0, 3: -310_000.0}
        matrix = compute_energy_matrix(_MONOMER, polymers, None, "RER")
        # IE^2_1 uses the monomer as the cyclic chain of degree 1, so equals E^2_2
        assert matrix.loc[2, 1] == pytest.approx(matrix.loc[2, 2])
        assert matrix.loc[2, 2] == pytest.approx((-205_000.0 - 2 * _MONOMER) / 1000)

    def test_initiator_ignored(self):
        polymers = {2: -205_000.0}
        with_init = compute_energy_matrix(_MONOMER, polymers, _INITIATOR, "RER")
        without = compute_energy_matrix(_MONOMER, polymers, None, "RER")
        assert with_init.loc[2, 2] == pytest.approx(without.loc[2, 2])

    def test_missing_intermediate_length(self):
        matrix = compute_energy_matrix(
            _MONOMER,
            {2: -205_000.0, 4: -412_000.0},
            None,
            "RER",
        )
        assert np.isnan(matrix.loc[4, 3])
        assert not np.isnan(matrix.loc[4, 2])

    def test_missing_values_give_nan(self):
        matrix = compute_energy_matrix(_MONOMER, {2: None}, None, "RER")
        assert np.isnan(matrix.loc[2, 0])


class TestHelpers:
    def test_compute_energy_matrices(self):
        matrices = compute_energy_matrices(
            monomer={"E": _MONOMER, "H": _MONOMER, "G": None},
            polymers={2: {"E": -205_000.0, "H": -205_000.0, "G": None}},
            initiator=None,
            reaction_type="RER",
            params=["E", "H", "G"],
        )
        assert set(matrices) == {"E", "H", "G"}
        assert np.isnan(matrices["G"].loc[2, 0])

    def test_low_r2_warning(self):
        noisy = compute_energy_matrix(
            _MONOMER,
            {2: -205_000.0, 3: -310_000.0, 4: -412_000.0},
            None,
            "RER",
        )
        assert noisy.loc[INF_LABEL, R2_COLUMN] < 0.8
        assert low_r2_warning(noisy)
        exact = compute_energy_matrix(
            _MONOMER,
            {n: n * _MONOMER - 2_000.0 * n + 3_000.0 for n in (2, 3, 4)},
            None,
            "RER",
        )
        assert not low_r2_warning(exact)

    def test_two_points_have_no_r2(self):
        _, r2 = extrapolate_to_infinite_length([2, 3], [-5.0, -6.0])
        assert r2 is None

    def test_format(self):
        matrix = compute_energy_matrix(
            _MONOMER,
            {2: -205_000.0, 3: -310_000.0, 4: -412_000.0},
            None,
            "RER",
        )
        text = format_energy_matrix(matrix)
        assert "m=0" in text
        assert "∞" in text
        assert "R²" in text
