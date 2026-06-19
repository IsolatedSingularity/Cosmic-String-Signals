"""Wake geometry construction and SO(3) coordinate transforms."""

import numpy as np
from scipy.spatial import ConvexHull


def rotation_x(theta: float) -> np.ndarray:
    """Rotation matrix about the x-axis by angle theta."""
    return np.array([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta), np.cos(theta)],
    ])


def rotation_y(theta: float) -> np.ndarray:
    """Rotation matrix about the y-axis by angle theta."""
    return np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)],
    ])


def rotation_z(theta: float) -> np.ndarray:
    """Rotation matrix about the z-axis by angle theta."""
    return np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ])


def apply_rotation(wedge: np.ndarray, angles: tuple[float, float, float]) -> np.ndarray:
    """Apply SO(3) rotation (Euler angles) to a set of 3D points.

    Parameters
    ----------
    wedge : np.ndarray
        Array of shape (N, 3) representing points in R^3.
    angles : tuple of float
        Euler angles (alpha, beta, gamma) for rotations about x, y, z axes.

    Returns
    -------
    np.ndarray
        Rotated points with the same shape as input.
    """
    center = np.mean(wedge, axis=0)
    shifted = wedge - center
    rotated = shifted @ rotation_x(angles[0]).T
    rotated = rotated @ rotation_y(angles[1]).T
    rotated = rotated @ rotation_z(angles[2]).T
    return rotated + center


def build_wake_wedge(
    tip: np.ndarray,
    deficit_angle: float,
    wake_depth: float,
    wake_length: float,
) -> np.ndarray:
    """Construct a 6-vertex wake wedge geometry in physical space.

    Parameters
    ----------
    tip : np.ndarray
        3D coordinates of the wake tip point.
    deficit_angle : float
        Opening angle of the wake wedge.
    wake_depth : float
        Radial depth of the wake in Mpc.
    wake_length : float
        Length of the wake along the string direction in Mpc.

    Returns
    -------
    np.ndarray
        Array of shape (6, 3) defining the wake wedge vertices.
    """
    end_points = np.array([
        [tip[0] + wake_depth * np.cos(deficit_angle / 2),
         tip[1] + wake_depth * np.sin(deficit_angle / 2),
         tip[2]],
        [tip[0] + wake_depth * np.cos(deficit_angle / 2),
         tip[1] - wake_depth * np.sin(deficit_angle / 2),
         tip[2]],
    ])
    projected = np.array([
        tip + [0, 0, wake_length],
        end_points[0] + [0, 0, wake_length],
        end_points[1] + [0, 0, wake_length],
    ])
    return np.array([tip, end_points[0], end_points[1],
                     projected[0], projected[1], projected[2]])


def point_in_hull(hull_points: np.ndarray, test_points: np.ndarray) -> bool:
    """Check whether test points lie inside the convex hull of hull_points.

    Uses the vertex-stability method: if adding the test points does not
    change the hull vertices, all test points are interior.

    Parameters
    ----------
    hull_points : np.ndarray
        Array of shape (M, 3) defining the convex hull.
    test_points : np.ndarray
        Array of shape (K, 3) of points to test.

    Returns
    -------
    bool
        True if all test points are inside the hull.
    """
    hull = ConvexHull(hull_points)
    combined = np.append(hull_points, test_points, axis=0)
    new_hull = ConvexHull(combined)
    return list(hull.vertices) == list(new_hull.vertices)
