"""Runs a single 1D FEM simulation and plots exact vs numerical solution."""

from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from fem1d.mesh import Mesh1D
from fem1d.solver import solve_poisson_1d
from fem1d.errors import compute_errors
from fem1d.manufactured import get_sin_benchmark, get_exponential_benchmark, get_quadratic_benchmark


def main():
    parser = argparse.ArgumentParser(description="1D Poisson FEM Single Simulation")
    parser.add_argument(
        "--case",
        choices=["sin", "poly_exp", "quadratic"],
        default="sin",
        help="Benchmark problem case",
    )
    parser.add_argument(
        "--elements",
        "-n",
        type=int,
        default=16,
        help="Number of elements in mesh",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/figures",
        help="Directory to save generated figures",
    )
    args = parser.parse_args()

    # Select benchmark
    benchmarks = {
        "sin": get_sin_benchmark(),
        "poly_exp": get_exponential_benchmark(),
        "quadratic": get_quadratic_benchmark(),
    }
    bench = benchmarks[args.case]

    # Generate mesh and solve
    mesh = Mesh1D.uniform(0.0, 1.0, args.elements)
    u_fem = solve_poisson_1d(mesh, bench.f, quad_order=5)

    # Compute errors
    errors = compute_errors(mesh, u_fem, bench.u, bench.u_prime, quad_order=5)

    print("=" * 60)
    print(f"1D Poisson FEM Solution: {bench.description}")
    print(f"Mesh: {args.elements} elements, {mesh.num_nodes} nodes (h = {mesh.h_max:.4e})")
    print("-" * 60)
    print(f"L2 Error:            {errors.l2_error:.6e}")
    print(f"H1-Semi Error:       {errors.h1_semi_error:.6e}")
    print(f"Full H1 Error:       {errors.h1_error:.6e}")
    print(f"Max Nodal Error:     {errors.linf_nodal_error:.6e}")
    print("=" * 60)

    # Plot solution
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    x_dense = np.linspace(0.0, 1.0, 500)
    u_dense = bench.u(x_dense)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    ax1.plot(x_dense, u_dense, "k-", label="Exact Solution $u(x)$", lw=1.5)
    ax1.plot(mesh.nodes, u_fem, "ro--", label=f"FEM Solution $u_h$ ($N={args.elements}$)", ms=5)
    ax1.set_ylabel("$u(x)$")
    ax1.set_title(f"1D Poisson FEM: {bench.description} ($N={args.elements}$)")
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend()

    # Pointwise error on nodes
    u_exact_nodes = bench.u(mesh.nodes)
    nodal_err = np.abs(u_exact_nodes - u_fem)
    ax2.stem(mesh.nodes, nodal_err, linefmt="b-", markerfmt="bo", basefmt="k-")
    ax2.set_xlabel("$x$")
    ax2.set_ylabel("Nodal Absolute Error $|u - u_h|$")
    ax2.set_yscale("log")
    ax2.grid(True, linestyle=":", alpha=0.6)

    fig.tight_layout()
    out_path = out_dir / f"solution_{args.case}_N{args.elements}.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved figure to: {out_path}")


if __name__ == "__main__":
    main()
