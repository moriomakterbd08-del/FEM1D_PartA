"""Verifies nodal exactness of 1D linear FEM and explains the Green's function property."""

from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from fem1d.mesh import Mesh1D
from fem1d.solver import solve_poisson_1d
from fem1d.manufactured import get_sin_benchmark, get_exponential_benchmark, get_quadratic_benchmark


def main():
    print("=" * 80)
    print("1D FEM NODAL EXACTNESS INVESTIGATION")
    print("=" * 80)
    print("Theorem:")
    print("In one space dimension, the Green's function G(x, y) for -u'' = delta(x-y)")
    print("with u(0) = u(1) = 0 is continuous and piecewise linear in y:")
    print("      G(x, y) = y(1 - x) for 0 <= y <= x")
    print("      G(x, y) = x(1 - y) for x <= y <= 1")
    print("When x = x_j is an interior mesh node, G(x_j, y) lies EXACTLY in the")
    print("finite element space V_h of continuous piecewise linear functions!")
    print("Therefore, by Galerkin orthogonality:")
    print("      u(x_j) - u_h(x_j) = a(u - u_h, G(x_j, .)) = 0")
    print("provided the load vector is integrated exactly.")
    print("=" * 80)

    benchmarks = [
        ("Quadratic: u(x) = x(1-x)", get_quadratic_benchmark(), 2),
        ("Sine: u(x) = sin(pi*x)", get_sin_benchmark(), 8),
        ("Polynomial-Exponential: u(x) = x(1-x)e^x", get_exponential_benchmark(), 8),
    ]

    for name, bench, q_order in benchmarks:
        print(f"\nBenchmark: {name}")
        print(f"{'N Elements':>12} {'Max Nodal Error |u(x_j) - u_h(x_j)|':>36}")
        print("-" * 52)
        for n_elems in [4, 8, 16, 32, 64]:
            mesh = Mesh1D.uniform(0.0, 1.0, n_elems)
            u_fem = solve_poisson_1d(mesh, bench.f, quad_order=q_order)
            u_exact = bench.u(mesh.nodes)
            max_nodal_err = np.max(np.abs(u_exact - u_fem))
            print(f"{n_elems:12d} {max_nodal_err:36.4e}")

    # Visualization of Nodal Exactness vs Inter-element Interpolation Error
    out_dir = Path("output/figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    n_demo = 6
    mesh = Mesh1D.uniform(0.0, 1.0, n_demo)
    bench = get_sin_benchmark()
    u_fem = solve_poisson_1d(mesh, bench.f, quad_order=8)

    # Dense sampling across entire domain
    x_dense = np.linspace(0.0, 1.0, 1000)
    u_exact_dense = bench.u(x_dense)

    # Evaluate piecewise linear u_h at all x_dense points
    u_h_dense = np.interp(x_dense, mesh.nodes, u_fem)
    pointwise_error = np.abs(u_exact_dense - u_h_dense)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Top plot: Solutions
    ax1.plot(x_dense, u_exact_dense, "k-", lw=2, label="Exact $u(x) = \\sin(\\pi x)$")
    ax1.plot(x_dense, u_h_dense, "r--", lw=1.5, label=f"FEM Solution $u_h$ ($N={n_demo}$)")
    ax1.plot(mesh.nodes, u_fem, "ro", ms=8, label="Mesh Nodes $(x_j, u_h(x_j))$")
    ax1.set_ylabel("$u(x)$", fontsize=12)
    ax1.set_title("Nodal Exactness Phenomenon in 1D FEM", fontsize=14, fontweight="bold")
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(fontsize=11)

    # Bottom plot: Pointwise error showing zeros at every node
    ax2.plot(x_dense, pointwise_error, "b-", lw=1.8, label="Pointwise Error $|u(x) - u_h(x)|$")
    ax2.plot(mesh.nodes, np.zeros_like(mesh.nodes), "ro", ms=8, label="Zero Error at Nodes")
    for xj in mesh.nodes:
        ax2.axvline(xj, color="gray", linestyle=":", alpha=0.4)
    ax2.set_xlabel("$x$", fontsize=12)
    ax2.set_ylabel("Error $|u(x) - u_h(x)|$", fontsize=12)
    ax2.set_title("Pointwise Error Across Domain: Vanishes at Nodes $x_j$", fontsize=12)
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(fontsize=11)

    fig.tight_layout()
    out_path = out_dir / "nodal_exactness.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"\nSaved nodal exactness illustration to: {out_path}")


if __name__ == "__main__":
    main()
