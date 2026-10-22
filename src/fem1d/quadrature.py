"""Gauss-Legendre numerical quadrature on reference and physical 1D intervals."""

from __future__ import annotations
import numpy as np


class GaussLegendre1D:
    """Provides Gauss-Legendre quadrature points and weights.

    Reference domain is [-1, 1].
    """

    @staticmethod
    def get_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
        """Returns (points, weights) on [-1, 1] for a given number of quadrature points.

        Args:
            order: Number of Gauss quadrature points (1 to 5).

        Returns:
            points: Array of Gauss points in [-1, 1].
            weights: Array of corresponding weights summing to 2.
        """
        if order == 1:
            points = np.array([0.0])
            weights = np.array([2.0])
        elif order == 2:
            pt = 1.0 / np.sqrt(3.0)
            points = np.array([-pt, pt])
            weights = np.array([1.0, 1.0])
        elif order == 3:
            pt = np.sqrt(3.0 / 5.0)
            points = np.array([-pt, 0.0, pt])
            weights = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
        elif order == 4:
            pt1 = np.sqrt((3.0 - 2.0 * np.sqrt(6.0 / 5.0)) / 7.0)
            pt2 = np.sqrt((3.0 + 2.0 * np.sqrt(6.0 / 5.0)) / 7.0)
            w1 = (18.0 + np.sqrt(30.0)) / 36.0
            w2 = (18.0 - np.sqrt(30.0)) / 36.0
            points = np.array([-pt2, -pt1, pt1, pt2])
            weights = np.array([w2, w1, w1, w2])
        elif order == 5:
            pt1 = (1.0 / 3.0) * np.sqrt(5.0 - 2.0 * np.sqrt(10.0 / 7.0))
            pt2 = (1.0 / 3.0) * np.sqrt(5.0 + 2.0 * np.sqrt(10.0 / 7.0))
            w1 = (322.0 + 13.0 * np.sqrt(70.0)) / 900.0
            w2 = (322.0 - 13.0 * np.sqrt(70.0)) / 900.0
            points = np.array([-pt2, -pt1, 0.0, pt1, pt2])
            weights = np.array([w2, w1, 128.0 / 225.0, w1, w2])
        else:
            # Fallback to numpy roots of Legendre polynomials for arbitrary orders
            points, weights = np.polynomial.legendre.leggauss(order)

        return points, weights

    @classmethod
    def transform_to_interval(
        cls, x_left: float, x_right: float, order: int = 3
    ) -> tuple[np.ndarray, np.ndarray]:
        """Maps reference quadrature points and weights from [-1, 1] to [x_left, x_right].

        Args:
            x_left: Left boundary of physical element.
            x_right: Right boundary of physical element.
            order: Quadrature order.

        Returns:
            x_physical: Coordinates in [x_left, x_right].
            w_physical: Scaled weights such that sum(w_physical) = x_right - x_left.
        """
        ref_points, ref_weights = cls.get_rule(order)
        h = x_right - x_left
        x_physical = 0.5 * (x_left + x_right) + 0.5 * h * ref_points
        w_physical = 0.5 * h * ref_weights
        return x_physical, w_physical
