"""Regression tests for theoretical convergence rates and nodal exactness."""

import pytest
import numpy as np
from fem1d.mesh import Mesh1D
from fem1d.solver import solve_poisson_1d
from fem1d.errors import compute_errors
from fem1d.manufactured import get_sin_benchmark, get_exponential_benchmark, get_quadratic_benchmark


@pytest.mark.parametrize("benchmark_fn", [get_sin_benchmark, get_exponential_benchmark])
def test_convergence_rates(benchmark_fn):
    bench = benchmark_fn()
    n_elements_list = [8, 16, 32, 64, 128]

    h_values = []
    l2_errors = []
    h1_errors = []

    for n_elems in n_elements_list:
        mesh = Mesh1D.uniform(0.0, 1.0, n_elems)
        u_fem = solve_poisson_1d(mesh, bench.f, quad_order=5)
        err = compute_errors(mesh, u_fem, bench.u, bench.u_prime, quad_order=5)

        h_values.append(err.h_max)
        l2_errors.append(err.l2_error)
        h1_errors.append(err.h1_semi_error)

    # Compute observed convergence rates: rate = log(e_k / e_{k-1}) / log(h_k / h_{k-1})
    log_h = np.log(h_values)
    log_l2 = np.log(l2_errors)
    log_h1 = np.log(h1_errors)

    # Linear regression slopes
    poly_l2 = np.polyfit(log_h, log_l2, 1)
    poly_h1 = np.polyfit(log_h, log_h1, 1)
    rate_l2 = poly_l2[0]
    rate_h1 = poly_h1[0]

    # Céa's lemma predicts rate 1 in H1, Aubin-Nitsche predicts rate 2 in L2
    assert rate_l2 >= 1.95, f"Expected L2 rate ~2, got {rate_l2:.3f}"
    assert rate_h1 >= 0.98, f"Expected H1 rate ~1, got {rate_h1:.3f}"


def test_nodal_exactness_quadratic():
    """For u(x) = x(1-x), f(x) = 2 is constant.

    In 1D, piecewise linear FEM with exact quadrature is exact at nodes.
    """
    bench = get_quadratic_benchmark()
    for n_elems in [4, 8, 16]:
        mesh = Mesh1D.uniform(0.0, 1.0, n_elems)
        u_fem = solve_poisson_1d(mesh, bench.f, quad_order=2)
        err = compute_errors(mesh, u_fem, bench.u, bench.u_prime, quad_order=2)

        # Nodal error should be around float machine epsilon
        assert err.linf_nodal_error < 1e-14, f"Nodal error too large: {err.linf_nodal_error}"
