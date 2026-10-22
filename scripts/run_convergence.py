"""Automated convergence study for 1D Poisson FEM verifying Céa's lemma and Aubin-Nitsche."""

from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from fem1d.mesh import Mesh1D
from fem1d.solver import solve_poisson_1d
from fem1d.errors import compute_errors
from fem1d.manufactured import get_sin_benchmark, get_exponential_benchmark


def run_benchmark_convergence(bench, n_elements_list):
    results = []
    for n in n_elements_list:
        mesh = Mesh1D.uniform(0.0, 1.0, n)
        u_fem = solve_poisson_1d(mesh, bench.f, quad_order=5)
        err = compute_errors(mesh, u_fem, bench.u, bench.u_prime, quad_order=5)
        results.append(err)
    return results


def print_table(bench_name, results):
    print(f"\nConvergence Table: {bench_name}")
    print("=" * 82)
    print(
        f"{'N':>5} {'h':>10} {'L2 Error':>14} {'Rate(L2)':>10} {'H1-Semi Err':>14} {'Rate(H1)':>10} {'Nodal Max':>12}"
    )
    print("-" * 82)

    for i, res in enumerate(results):
        if i == 0:
            rate_l2_str = "    -   "
            rate_h1_str = "    -   "
        else:
            prev = results[i - 1]
            rate_l2 = np.log(res.l2_error / prev.l2_error) / np.log(res.h_max / prev.h_max)
            rate_h1 = np.log(res.h1_semi_error / prev.h1_semi_error) / np.log(res.h_max / prev.h_max)
            rate_l2_str = f"{rate_l2:8.2f}"
            rate_h1_str = f"{rate_h1:8.2f}"

        n_elem = int(round(1.0 / res.h_max))
        print(
            f"{n_elem:5d} {res.h_max:10.4e} {res.l2_error:14.6e} {rate_l2_str} {res.h1_semi_error:14.6e} {rate_h1_str} {res.linf_nodal_error:12.4e}"
        )
    print("=" * 82)


def main():
    n_elements_list = [4, 8, 16, 32, 64, 128, 256]
    benchmarks = [
        ("Trigonometric: u(x) = sin(pi*x)", get_sin_benchmark()),
        ("Polynomial-Exponential: u(x) = x(1-x)e^x", get_exponential_benchmark()),
    ]

    out_dir = Path("output/figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, (label, bench) in zip(axes, benchmarks):
        results = run_benchmark_convergence(bench, n_elements_list)
        print_table(label, results)

        h_vals = np.array([r.h_max for r in results])
        l2_errs = np.array([r.l2_error for r in results])
        h1_errs = np.array([r.h1_semi_error for r in results])

        # Fitted slopes
        slope_l2 = np.polyfit(np.log(h_vals), np.log(l2_errs), 1)[0]
        slope_h1 = np.polyfit(np.log(h_vals), np.log(h1_errs), 1)[0]

        # Log-log plot
        ax.loglog(h_vals, l2_errs, "s-", color="navy", label=f"$L^2$ error (slope = {slope_l2:.2f})", lw=1.8)
        ax.loglog(h_vals, h1_errs, "o-", color="crimson", label=f"$H^1_0$ error (slope = {slope_h1:.2f})", lw=1.8)

        # Reference triangles
        ref_h = h_vals
        ref_l2 = l2_errs[0] * (ref_h / ref_h[0]) ** 2
        ref_h1 = h1_errs[0] * (ref_h / ref_h[0]) ** 1
        ax.loglog(ref_h, ref_l2, "k--", alpha=0.5, label="Theoretical $\\mathcal{O}(h^2)$")
        ax.loglog(ref_h, ref_h1, "k:", alpha=0.5, label="Theoretical $\\mathcal{O}(h)$")

        ax.set_xlabel("Mesh size $h$", fontsize=12)
        ax.set_ylabel("Error Norm", fontsize=12)
        ax.set_title(label, fontsize=12)
        ax.grid(True, which="both", linestyle=":", alpha=0.6)
        ax.legend(fontsize=10)

    fig.suptitle("1D Poisson FEM: Convergence Study (Log-Log)", fontsize=14, fontweight="bold")
    fig.tight_layout()

    out_path = out_dir / "convergence_study.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"\nSaved convergence figures to: {out_path}")


if __name__ == "__main__":
    main()
