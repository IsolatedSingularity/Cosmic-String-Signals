"""Matched filtering and signal extraction utilities."""

import numpy as np
from scipy import signal


def unfold(grid: np.ndarray, orientation: str = "horizontal") -> list[float]:
    """Unfold a 2D array into a 1D sequence for correlation analysis.

    Parameters
    ----------
    grid : np.ndarray
        2D array to unfold.
    orientation : str
        Either 'horizontal' (row-major) or 'vertical' (column-major).

    Returns
    -------
    list of float
        Flattened values in the specified traversal order.
    """
    copy = np.copy(grid)
    if orientation == "vertical":
        copy = copy.T
    return copy.ravel().tolist()


def matched_filter_1d(template: np.ndarray, data: np.ndarray) -> np.ndarray:
    """Compute 1D cross-correlation between a template and data signal.

    Parameters
    ----------
    template : np.ndarray
        1D reference signal (the expected wake profile).
    data : np.ndarray
        1D observed signal (may contain noise).

    Returns
    -------
    np.ndarray
        Full cross-correlation output.
    """
    return np.correlate(template, data, mode="full")


def matched_filter_2d(template: np.ndarray, data: np.ndarray) -> np.ndarray:
    """Compute 2D convolution between a template and data map.

    Parameters
    ----------
    template : np.ndarray
        2D reference map.
    data : np.ndarray
        2D observed map.

    Returns
    -------
    np.ndarray
        Full 2D convolution output.
    """
    return signal.convolve2d(template, data, mode="full")
