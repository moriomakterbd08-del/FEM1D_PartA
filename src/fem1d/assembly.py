"""Element-by-element assembly of stiffness matrix and load vector in 1D."""

from __future__ import annotations
from typing import Callable
import numpy as np
import scipy.sparse as sp

from .mesh import Mesh1D
from .quadrature import GaussLegendre1D
from .basis import LinearBasis1D


def compute_element_stiffness(h_e: float) -> np.ndarray:
    """Computes the 2x2 element stiffness matrix for -u'' on element of length h_e:

    K^e = (1 / h_e) * [[ 1, -1],
                       [-1,  1]]
    """
    return (1.0 / h_e) * np.array([[1.0, -1.0], [-1.0, 1.0]])


def compute_element_load(
    f: Callable[[np.ndarray], np.ndarray],
    x_left: float,
    x_right: float,
    quad_order: int = 4,
) -> np.ndarray:
    """Computes the 2x1 element load vector:

    F_a^e = int_{x_left}^{x_right} f(x) * phi_a^e(x) dx
    using Gauss-Legendre quadrature.
    """
    x_quad, w_quad = GaussLegendre1D.transform_to_interval(x_left, x_right, order=quad_order)
    f_vals = f(x_quad)
    phi_vals = LinearBasis1D.eval_physical(x_quad, x_left, x_right)  # (N_quad, 2)

    # Vectorized integration: sum_q w_q * f(x_q) * phi_{a}(x_q)
    F_elem = np.sum(w_quad[:, None] * f_vals[:, None] * phi_vals, axis=0)
    return F_elem


def assemble_poisson_system(
    mesh: Mesh1D,
    f: Callable[[np.ndarray], np.ndarray],
    quad_order: int = 4,
) -> tuple[sp.csr_matrix, np.ndarray]:
    """Assembles the global stiffness matrix K and load vector F element by element.

    Args:
        mesh: 1D finite element mesh.
        f: Source function f(x).
        quad_order: Order of Gauss quadrature for load vector integration.

    Returns:
        K: Global stiffness matrix as a scipy.sparse.csr_matrix of shape (num_nodes, num_nodes).
        F: Global load vector of shape (num_nodes,).
    """
    n_nodes = mesh.num_nodes
    n_elems = mesh.num_elements

    # Use lil_matrix for efficient incremental assembly
    K = sp.lil_matrix((n_nodes, n_nodes), dtype=float)
    F = np.zeros(n_nodes, dtype=float)

    for e in range(n_elems):
        elem_nodes = mesh.elements[e]
        x_left, x_right = mesh.get_element_coords(e)
        h_e = x_right - x_left

        # Element matrices
        K_e = compute_element_stiffness(h_e)
        F_e = compute_element_load(f, x_left, x_right, quad_order=quad_order)

        # Assemble into global system
        for i_local in range(2):
            i_global = elem_nodes[i_local]
            F[i_global] += F_e[i_local]
            for j_local in range(2):
                j_global = elem_nodes[j_local]
                K[i_global, j_global] += K_e[i_local, j_local]

    return K.tocsr(), F
