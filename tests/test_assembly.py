"""Unit tests for element and global FEM matrix assembly."""

import pytest
import numpy as np
import scipy.sparse as sp
from fem1d.mesh import Mesh1D
from fem1d.assembly import compute_element_stiffness, assemble_poisson_system
from fem1d.boundary import apply_homogeneous_dirichlet


def test_element_stiffness():
    h = 0.2
    K_e = compute_element_stiffness(h)
    expected = 5.0 * np.array([[1.0, -1.0], [-1.0, 1.0]])
    assert np.allclose(K_e, expected)
    # Row sums of element stiffness must be zero (conservation / constant nullspace)
    assert np.allclose(np.sum(K_e, axis=1), 0.0)


def test_global_stiffness_uniform():
    n_elems = 4
    mesh = Mesh1D.uniform(0.0, 1.0, n_elems)
    h = 0.25

    # Constant f = 0 to inspect K only
    f_zero = lambda x: np.zeros_like(x)
    K, _ = assemble_poisson_system(mesh, f_zero)
    K_dense = K.toarray()

    # Diagonal entries: interior should be 2/h, boundaries 1/h
    assert K_dense[0, 0] == pytest.approx(1.0 / h)
    assert K_dense[-1, -1] == pytest.approx(1.0 / h)
    for i in range(1, n_elems):
        assert K_dense[i, i] == pytest.approx(2.0 / h)
        assert K_dense[i, i - 1] == pytest.approx(-1.0 / h)
        assert K_dense[i, i + 1] == pytest.approx(-1.0 / h)

    # Symmetry
    assert np.allclose(K_dense, K_dense.T)


def test_reduced_system_positive_definiteness():
    mesh = Mesh1D.uniform(0.0, 1.0, 8)
    f = lambda x: np.ones_like(x)
    K, F = assemble_poisson_system(mesh, f)
    K_red, F_red, idxs = apply_homogeneous_dirichlet(K, F)

    assert K_red.shape == (7, 7)
    assert len(F_red) == 7
    eigvals = np.linalg.eigvalsh(K_red.toarray())
    # All eigenvalues of reduced stiffness matrix must be strictly positive (coercivity)
    assert np.all(eigvals > 0.0)
