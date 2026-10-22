"""Application of boundary conditions to 1D FEM linear systems."""

from __future__ import annotations
import numpy as np
import scipy.sparse as sp


def apply_homogeneous_dirichlet(
    K: sp.csr_matrix, F: np.ndarray
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray]:
    """Applies homogeneous Dirichlet boundary conditions u(0) = 0, u(1) = 0

    by eliminating the boundary degrees of freedom (first and last nodes).

    Args:
        K: Full assembled stiffness matrix of shape (N+1, N+1).
        F: Full assembled load vector of shape (N+1,).

    Returns:
        K_reduced: Reduced stiffness matrix for interior nodes of shape (N-1, N-1).
        F_reduced: Reduced load vector for interior nodes of shape (N-1,).
        interior_indices: Array of global node indices [1, 2, ..., N-1].
    """
    n_nodes = F.shape[0]
    interior_indices = np.arange(1, n_nodes - 1)

    # Slice sparse matrix
    K_reduced = K[interior_indices, :][:, interior_indices].tocsr()
    F_reduced = F[interior_indices]

    return K_reduced, F_reduced, interior_indices


def reconstruct_full_solution(
    u_reduced: np.ndarray,
    n_nodes: int,
    u_left: float = 0.0,
    u_right: float = 0.0,
) -> np.ndarray:
    """Reconstructs the full solution vector U of size n_nodes including boundary values."""
    u_full = np.zeros(n_nodes, dtype=float)
    u_full[0] = u_left
    u_full[-1] = u_right
    u_full[1:-1] = u_reduced
    return u_full
