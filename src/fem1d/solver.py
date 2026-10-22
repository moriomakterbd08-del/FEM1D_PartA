"""Linear solver for 1D Poisson FEM problems."""

from __future__ import annotations
from typing import Callable
import numpy as np
import scipy.sparse.linalg as spla

from .mesh import Mesh1D
from .assembly import assemble_poisson_system
from .boundary import apply_homogeneous_dirichlet, reconstruct_full_solution


def solve_poisson_1d(
    mesh: Mesh1D,
    f: Callable[[np.ndarray], np.ndarray],
    quad_order: int = 4,
) -> np.ndarray:
    """Solves -u'' = f on [0, 1] with u(0) = u(1) = 0 using continuous P1 FEM.

    Args:
        mesh: 1D mesh partition of [0, 1].
        f: Source function f(x).
        quad_order: Gauss quadrature order for load vector integration.

    Returns:
        u_full: Solution vector of nodal values, shape (mesh.num_nodes,).
    """
    K_global, F_global = assemble_poisson_system(mesh, f, quad_order=quad_order)
    K_reduced, F_reduced, _ = apply_homogeneous_dirichlet(K_global, F_global)

    # Solve the reduced symmetric positive definite tridiagonal linear system
    u_reduced = spla.spsolve(K_reduced, F_reduced)

    # Reconstruct full nodal vector with boundary zeros
    u_full = reconstruct_full_solution(u_reduced, mesh.num_nodes, u_left=0.0, u_right=0.0)
    return u_full
