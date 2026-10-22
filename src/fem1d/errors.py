"""Error norm calculations for 1D Finite Element solutions."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable
import numpy as np

from .mesh import Mesh1D
from .quadrature import GaussLegendre1D
from .basis import LinearBasis1D


@dataclass
class FEMErrors:
    """Stores error metrics for a FEM solution."""
    l2_error: float
    h1_semi_error: float
    h1_error: float
    linf_nodal_error: float
    h_max: float


def compute_errors(
    mesh: Mesh1D,
    u_fem: np.ndarray,
    u_exact: Callable[[np.ndarray], np.ndarray],
    u_prime_exact: Callable[[np.ndarray], np.ndarray],
    quad_order: int = 5,
) -> FEMErrors:
    """Computes L2, H1-semi, full H1, and maximum nodal errors between the

    FEM solution u_fem and the continuous analytical solution u_exact.

    Integration over each element e = [x_left, x_right] uses Gauss quadrature:
        ||u - u_h||_{L^2}^2 = sum_e int_e |u(x) - u_h(x)|^2 dx
        |u - u_h|_{H^1}^2   = sum_e int_e |u'(x) - u_h'(x)|^2 dx

    Args:
        mesh: 1D finite element mesh.
        u_fem: Array of nodal values of the FEM solution, shape (num_nodes,).
        u_exact: Analytical solution function u(x).
        u_prime_exact: Analytical derivative function u'(x).
        quad_order: Number of Gauss quadrature points per element.

    Returns:
        FEMErrors object with computed norms.
    """
    total_l2_sq = 0.0
    total_h1_semi_sq = 0.0

    for e in range(mesh.num_elements):
        node_left, node_right = mesh.elements[e]
        x_left, x_right = mesh.get_element_coords(e)
        u_elem = np.array([u_fem[node_left], u_fem[node_right]])

        # Quadrature points and weights on this element
        x_q, w_q = GaussLegendre1D.transform_to_interval(x_left, x_right, order=quad_order)

        # Basis functions and derivatives at quadrature points
        phi_q = LinearBasis1D.eval_physical(x_q, x_left, x_right)         # (N_q, 2)
        dphi_q = LinearBasis1D.eval_physical_deriv(x_q, x_left, x_right)  # (N_q, 2)

        # Evaluate FEM approximation and derivative at quad points
        u_h_q = phi_q @ u_elem
        du_h_q = dphi_q @ u_elem

        # Exact values at quad points
        u_ex_q = u_exact(x_q)
        du_ex_q = u_prime_exact(x_q)

        # Integrate squared errors
        total_l2_sq += np.sum(w_q * (u_ex_q - u_h_q) ** 2)
        total_h1_semi_sq += np.sum(w_q * (du_ex_q - du_h_q) ** 2)

    l2_error = float(np.sqrt(total_l2_sq))
    h1_semi_error = float(np.sqrt(total_h1_semi_sq))
    h1_error = float(np.sqrt(total_l2_sq + total_h1_semi_sq))

    # Nodal maximum error
    nodal_exact = u_exact(mesh.nodes)
    linf_nodal_error = float(np.max(np.abs(nodal_exact - u_fem)))

    return FEMErrors(
        l2_error=l2_error,
        h1_semi_error=h1_semi_error,
        h1_error=h1_error,
        linf_nodal_error=linf_nodal_error,
        h_max=mesh.h_max,
    )
