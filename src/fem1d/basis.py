"""1D Lagrange P1 linear basis functions and derivatives."""

from __future__ import annotations
import numpy as np


class LinearBasis1D:
    """Standard continuous piecewise linear (P1 / hat) shape functions on an element [x_0, x_1]."""

    @staticmethod
    def eval_reference(xi: np.ndarray) -> np.ndarray:
        """Evaluates shape functions on the reference element xi in [-1, 1].

        Returns:
            phi: Array of shape (len(xi), 2), where
                 phi[:, 0] = (1 - xi) / 2
                 phi[:, 1] = (1 + xi) / 2
        """
        xi = np.asarray(xi, dtype=float)
        phi_0 = 0.5 * (1.0 - xi)
        phi_1 = 0.5 * (1.0 + xi)
        return np.column_stack((phi_0, phi_1))

    @staticmethod
    def eval_reference_deriv(xi: np.ndarray) -> np.ndarray:
        """Evaluates d(phi)/d(xi) on reference element xi in [-1, 1].

        Returns:
            dphi_dxi: Array of shape (len(xi), 2), where
                      dphi_dxi[:, 0] = -0.5
                      dphi_dxi[:, 1] = +0.5
        """
        xi = np.asarray(xi, dtype=float)
        dphi_0 = np.full_like(xi, -0.5)
        dphi_1 = np.full_like(xi, 0.5)
        return np.column_stack((dphi_0, dphi_1))

    @staticmethod
    def eval_physical(x: np.ndarray, x_left: float, x_right: float) -> np.ndarray:
        """Evaluates physical shape functions on [x_left, x_right].

        phi_0(x) = (x_right - x) / h
        phi_1(x) = (x - x_left) / h
        """
        x = np.asarray(x, dtype=float)
        h = x_right - x_left
        phi_0 = (x_right - x) / h
        phi_1 = (x - x_left) / h
        return np.column_stack((phi_0, phi_1))

    @staticmethod
    def eval_physical_deriv(x: np.ndarray, x_left: float, x_right: float) -> np.ndarray:
        """Evaluates d(phi)/dx on [x_left, x_right].

        d(phi_0)/dx = -1 / h
        d(phi_1)/dx = +1 / h
        """
        x = np.asarray(x, dtype=float)
        h = x_right - x_left
        dphi_0 = np.full_like(x, -1.0 / h)
        dphi_1 = np.full_like(x, 1.0 / h)
        return np.column_stack((dphi_0, dphi_1))
