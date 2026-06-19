"""Tests for filtering module: unfolding and matched filtering."""

import numpy as np
import pytest
from cosmic_string_signals.filtering import (
    unfold,
    matched_filter_1d,
    matched_filter_2d,
)


class TestUnfold:
    """Verify array unfolding in both orientations."""

    def test_horizontal_unfold_length(self):
        grid = np.ones((5, 10))
        result = unfold(grid, "horizontal")
        assert len(result) == 50

    def test_vertical_unfold_length(self):
        grid = np.ones((5, 10))
        result = unfold(grid, "vertical")
        assert len(result) == 50

    def test_horizontal_is_row_major(self):
        grid = np.array([[1, 2], [3, 4]])
        result = unfold(grid, "horizontal")
        assert result == [1, 2, 3, 4]

    def test_vertical_is_column_major(self):
        grid = np.array([[1, 2], [3, 4]])
        result = unfold(grid, "vertical")
        assert result == [1, 3, 2, 4]

    def test_does_not_mutate_input(self):
        grid = np.array([[1.0, 2.0], [3.0, 4.0]])
        original = grid.copy()
        unfold(grid, "horizontal")
        np.testing.assert_array_equal(grid, original)


class TestMatchedFilter1D:
    """Verify 1D cross-correlation properties."""

    def test_autocorrelation_peak_at_center(self):
        template = np.array([0, 0, 1, 0, 0], dtype=float)
        result = matched_filter_1d(template, template)
        peak_idx = np.argmax(result)
        assert peak_idx == len(result) // 2

    def test_output_length(self):
        a = np.ones(10)
        b = np.ones(15)
        result = matched_filter_1d(a, b)
        assert len(result) == len(a) + len(b) - 1

    def test_correlation_with_noise_lower_than_self(self):
        np.random.seed(42)
        template = np.sin(np.linspace(0, 2 * np.pi, 100))
        noise = np.random.randn(100) * 0.1
        self_corr = matched_filter_1d(template, template)
        noise_corr = matched_filter_1d(template, noise)
        assert np.max(np.abs(self_corr)) > np.max(np.abs(noise_corr))


class TestMatchedFilter2D:
    """Verify 2D convolution properties."""

    def test_output_shape(self):
        a = np.ones((5, 5))
        b = np.ones((3, 3))
        result = matched_filter_2d(a, b)
        assert result.shape == (7, 7)

    def test_self_convolution_peak(self):
        template = np.zeros((5, 5))
        template[2, 2] = 1.0
        result = matched_filter_2d(template, template)
        peak = np.unravel_index(np.argmax(result), result.shape)
        assert peak == (4, 4)
