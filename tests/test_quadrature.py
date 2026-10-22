"""Unit tests for Gauss-Legendre quadrature accuracy."""

import pytest
import numpy as np
from fem1d.quadrature import GaussLegendre1D


@pytest.mark.parametrize("order", [1, 2, 3, 4, 5])
def test_weights_sum_to_interval_length(order):
    # On reference [-1, 1], sum of weights must be 2.0
    _, weights = GaussLegendre1D.get_rule(order)
    assert np.sum(weights) == pytest.approx(2.0, abs=1e-14)

    # On physical interval [a, b], sum of weights must be b - a
    a, b = 0.3, 0.9
    _, w_phys = GaussLegendre1D.transform_to_interval(a, b, order=order)
    assert np.sum(w_phys) == pytest.approx(b - a, abs=1e-14)


def test_polynomial_integration_exactness():
    # 2-point Gauss integrates up to degree 3 polynomials exactly
    # int_{-1}^1 (x^3 - 2x^2 + 5x + 3) dx = [-2*(2/3) + 3*2] = -4/3 + 6 = 14/3
    pts, wts = GaussLegendre1D.get_rule(order=2)
    poly = lambda x: x**3 - 2 * x**2 + 5 * x + 3
    integral = np.sum(wts * poly(pts))
    assert integral == pytest.approx(14.0 / 3.0, abs=1e-13)

    # 3-point Gauss integrates up to degree 5 polynomials exactly
    pts3, wts3 = GaussLegendre1D.get_rule(order=3)
    poly5 = lambda x: 4 * x**5 - 3 * x**4 + 2 * x - 1
    # int_{-1}^1 = 0 - 3*(2/5) + 0 - 2 = -6/5 - 2 = -16/5 = -3.2
    integral3 = np.sum(wts3 * poly5(pts3))
    assert integral3 == pytest.approx(-3.2, abs=1e-13)
