"""Unit tests for 1D Mesh generation and properties."""

import pytest
import numpy as np
from fem1d.mesh import Mesh1D


def test_uniform_mesh():
    mesh = Mesh1D.uniform(0.0, 1.0, 4)
    assert mesh.num_elements == 4
    assert mesh.num_nodes == 5
    assert np.allclose(mesh.nodes, [0.0, 0.25, 0.5, 0.75, 1.0])
    assert np.allclose(mesh.element_lengths, 0.25)
    assert mesh.h_max == pytest.approx(0.25)
    assert mesh.h_min == pytest.approx(0.25)


def test_element_connectivity():
    mesh = Mesh1D.uniform(0.0, 2.0, 3)
    expected_elements = np.array([[0, 1], [1, 2], [2, 3]])
    assert np.array_equal(mesh.elements, expected_elements)
    left, right = mesh.get_element_coords(1)
    assert left == pytest.approx(2.0 / 3.0)
    assert right == pytest.approx(4.0 / 3.0)


def test_invalid_nodes():
    with pytest.raises(ValueError, match="strictly increasing"):
        Mesh1D.from_nodes(np.array([0.0, 0.5, 0.4, 1.0]))
