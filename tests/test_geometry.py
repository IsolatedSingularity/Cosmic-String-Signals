"""Tests for geometry module: rotations, wake construction, hull checks."""

import numpy as np
import pytest
from cosmic_string_signals.geometry import (
    rotation_x,
    rotation_y,
    rotation_z,
    apply_rotation,
    build_wake_wedge,
    point_in_hull,
)


class TestRotationMatrices:
    """Verify SO(3) rotation matrix properties."""

    def test_rotation_x_identity(self):
        R = rotation_x(0.0)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-15)

    def test_rotation_y_identity(self):
        R = rotation_y(0.0)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-15)

    def test_rotation_z_identity(self):
        R = rotation_z(0.0)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-15)

    def test_rotation_x_orthogonal(self):
        R = rotation_x(np.pi / 4)
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-15)

    def test_rotation_y_orthogonal(self):
        R = rotation_y(np.pi / 3)
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-15)

    def test_rotation_z_orthogonal(self):
        R = rotation_z(np.pi / 6)
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-15)

    def test_rotation_determinant_is_one(self):
        for angle in [0.1, 0.5, 1.0, np.pi]:
            assert abs(np.linalg.det(rotation_x(angle)) - 1.0) < 1e-14
            assert abs(np.linalg.det(rotation_y(angle)) - 1.0) < 1e-14
            assert abs(np.linalg.det(rotation_z(angle)) - 1.0) < 1e-14

    def test_rotation_2pi_returns_identity(self):
        for rot_fn in [rotation_x, rotation_y, rotation_z]:
            R = rot_fn(2 * np.pi)
            np.testing.assert_allclose(R, np.eye(3), atol=1e-14)


class TestApplyRotation:
    """Verify composite rotation application."""

    def test_zero_rotation_preserves_points(self):
        pts = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        result = apply_rotation(pts, (0.0, 0.0, 0.0))
        np.testing.assert_allclose(result, pts, atol=1e-14)

    def test_rotation_preserves_centroid(self):
        pts = np.random.randn(6, 3)
        result = apply_rotation(pts, (0.5, 1.0, 1.5))
        np.testing.assert_allclose(
            np.mean(result, axis=0), np.mean(pts, axis=0), atol=1e-14
        )

    def test_rotation_preserves_distances(self):
        pts = np.random.randn(6, 3)
        result = apply_rotation(pts, (0.3, 0.7, 1.1))
        center = np.mean(pts, axis=0)
        original_dists = np.linalg.norm(pts - center, axis=1)
        rotated_dists = np.linalg.norm(result - center, axis=1)
        np.testing.assert_allclose(original_dists, rotated_dists, atol=1e-13)


class TestBuildWakeWedge:
    """Verify wake wedge geometry construction."""

    def test_wedge_has_six_vertices(self):
        tip = np.array([0.0, 0.0, 0.0])
        wedge = build_wake_wedge(tip, np.pi / 3, 1.0, 2.0)
        assert wedge.shape == (6, 3)

    def test_wedge_tip_is_first_vertex(self):
        tip = np.array([1.0, 2.0, 3.0])
        wedge = build_wake_wedge(tip, np.pi / 4, 0.5, 1.0)
        np.testing.assert_array_equal(wedge[0], tip)

    def test_wedge_projected_points_offset_by_length(self):
        tip = np.array([0.0, 0.0, 0.0])
        length = 5.0
        wedge = build_wake_wedge(tip, np.pi / 3, 1.0, length)
        # Projected points (indices 3, 4, 5) should have z = base_z + length
        for i in range(3):
            assert abs(wedge[i + 3][2] - (wedge[i][2] + length)) < 1e-14


class TestPointInHull:
    """Verify convex hull membership checks."""

    def test_interior_point_detected(self):
        cube = np.array([
            [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1],
        ], dtype=float)
        assert point_in_hull(cube, np.array([[0.5, 0.5, 0.5]]))

    def test_exterior_point_rejected(self):
        cube = np.array([
            [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1],
        ], dtype=float)
        assert not point_in_hull(cube, np.array([[5.0, 5.0, 5.0]]))
