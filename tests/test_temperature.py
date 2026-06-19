"""Tests for temperature module: brightness temperature models."""

import numpy as np
import pytest
from cosmic_string_signals.temperature import (
    photon_cmb_temperature,
    kinetic_temperature,
    gas_temperature,
    brightness_temperature,
)


class TestPhotonCMBTemperature:
    """Verify CMB photon temperature scaling."""

    def test_present_day_temperature(self):
        # At z=0, T_gamma should be approximately 2.725 K
        assert abs(photon_cmb_temperature(0) - 2.72548) < 1e-5

    def test_scales_with_redshift(self):
        z = 10.0
        expected = 2.72548 * (1 + z)
        assert abs(photon_cmb_temperature(z) - expected) < 1e-10

    def test_monotonically_increasing(self):
        temps = [photon_cmb_temperature(z) for z in range(100)]
        assert all(temps[i] < temps[i + 1] for i in range(len(temps) - 1))


class TestKineticTemperature:
    """Verify kinetic temperature model."""

    def test_positive_at_recombination(self):
        t = kinetic_temperature(1100, 3e-7, 0.7)
        assert t > 0

    def test_decreases_with_redshift(self):
        # T_kin ~ 1/(z+1), so should decrease with increasing z
        t1 = kinetic_temperature(10, 3e-7, 0.7)
        t2 = kinetic_temperature(100, 3e-7, 0.7)
        assert t1 > t2


class TestGasTemperature:
    """Verify gas temperature scaling."""

    def test_positive(self):
        assert gas_temperature(10) > 0

    def test_quadratic_scaling(self):
        z = 5.0
        expected = 0.02 * (z + 1) ** 2
        assert abs(gas_temperature(z) - expected) < 1e-14


class TestBrightnessTemperature:
    """Verify brightness temperature computation."""

    def test_returns_finite(self):
        result = brightness_temperature(30, 3e-7, 0.7)
        assert np.isfinite(result)

    def test_returns_float(self):
        result = brightness_temperature(50, 3e-7, 0.5)
        assert isinstance(result, (float, np.floating))
