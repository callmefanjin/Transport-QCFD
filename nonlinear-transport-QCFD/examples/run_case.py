"""Run the copied solver and compare phit with a finite-difference solution.

The five copied source modules are deliberately left unchanged. This driver
overrides case parameters for demonstrations and evaluates the existing Pauli
product using expanded gates or the equivalent direct statevector action.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import platform
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qiskit
import scipy
from qiskit.quantum_info import Statevector
from scipy.integrate import solve_ivp

import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import carleman
import evolution
import schrodingerisation
import transport
from pauli_action import evolve_pauli_product
from scientific_plot_style import configure_matplotlib, save_png_pdf, style_axes, style_legend




def read_source_parameters():
    tree = ast.parse((ROOT / "src" / "main.py").read_text(encoding="utf-8-sig"))
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main")
    parameters = {}
    for node in function.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                parameters[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                pass
    return parameters


def relative_error(value, reference):
    denominator = np.linalg.norm(reference)
    return float(np.linalg.norm(value - reference) / denominator) if denominator else None


def source_status():
    manifest = json.loads((ROOT / "source_manifest.json").read_text(encoding="utf-8"))
    checks = {name: hashlib.sha256((ROOT / "src" / name).read_bytes()).hexdigest() == expected
              for name, expected in manifest["sha256"].items()}
    return {"source_files_match_snapshot": checks}


def make_plots(output):
    configure_matplotlib()
    data = np.load(output / "results.npz")
    x = data["x"]
    curves = [
        ("phi_nonlinear", "Finite difference", "#222222", "o", "-"),
        ("phi_trotter", "phit", "#CC0077", "D", "--"),
    ]
    fig, ax = plt.subplots(figsize=(11.5, 7.6))
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.16, top=0.80)
    for key, label, color, marker, linestyle in curves:
        ax.plot(x, data[key].real, color=color, marker=marker, markersize=7,
                markerfacecolor="white", linewidth=2, linestyle=linestyle, label=label)
    style_axes(ax, xlabel=r"$x$", ylabel=r"$\phi$")
    style_legend(ax, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    save_png_pdf(fig, output / "solution_comparison")
    plt.close(fig)

    (output / "figure_captions.txt").write_text(
        "solution_comparison: real parts at the final time. Finite difference = "
        "the same periodic central spatial differences, with the full linear and "
        "quadratic source terms, integrated directly in time using Radau. "
        "phit = the original Carleman + WPT + Trotter solution recovered at p=0. "
        "Lines join the discrete grid points without smoothing.\n",
        encoding="utf-8")


def run(args):
    output = ROOT / "outputs" / args.name
    if args.plot_only:
        make_plots(output)
        return
    if output.exists():
        raise FileExistsError(f"Choose a new --name; existing results are preserved: {output}")
    source = read_source_parameters()
    config = dict(source)
    for key, value in (("start", args.start), ("end", args.end)):
        if value is not None:
            config[key] = value
    config["deltax"] = (config["end"] - config["start"]) / args.nx
    if args.mu is not None:
        config["initial_condition"] = "exp(-(x-mu)**2)"
        config["mu"] = args.mu
    for key, value in (("Np", args.np), ("r", args.r), ("T", args.time)):
        if value is not None:
            config[key] = value
    if args.nx < 3 or config["K"] < 1 or config["Np"] < 2 or config["Np"] % 2 or config["r"] < 1:
        raise ValueError("Use Nx >= 3, K >= 1, positive r, and an even Np >= 2.")
    output.mkdir(parents=True)
    (output / "parameters.json").write_text(json.dumps({"source": source, "run": config}, indent=2), encoding="utf-8")
    started = time.perf_counter()
    x = transport.build_grid(config["start"], config["end"], config["deltax"])
    if len(x) != args.nx:
        raise ValueError("Floating-point grid endpoint produced an unexpected point count.")
    phi0 = (np.exp(-(x - args.mu)**2) if args.mu is not None
            else transport.initial_condition(x))
    A = transport.linear_matrix(x, config["u"], config["D"], config["alpha"])
    F = transport.nonlinear_part(x, config["kappa"])
    Gamma = carleman.build_Gamma(A, F, config["K"])
    Y0 = carleman.build_Y(phi0, config["K"])
    H1, H2 = schrodingerisation.hermitian_split(Gamma)
    spectrum = np.linalg.eigvalsh(H1)
    if spectrum[-1] > 1e-10:
        raise ValueError("This runner uses p=0; a positive eigenvalue of H1 requires a different recovery strategy.")
    p, dp = schrodingerisation.build_p_grid(config["Lp"], config["Np"])
    psi0, eta = schrodingerisation.WPT_method(Y0, config["Lp"], config["Np"])
    W0 = np.outer(Y0, np.exp(-np.abs(p)))
    assert np.allclose(schrodingerisation.inverse_fft(psi0), W0)
    M = len(Y0)
    dim = 1 << (M - 1).bit_length()
    psi_t = np.zeros_like(psi0)
    tails = []
    norms = []
    raw_gates = []
    expanded_gates = []
    print(f"Run: Nx={len(x)}, K={config['K']}, M={M}, padded={dim}, Np={len(p)}, r={config['r']}, T={config['T']}", flush=True)
    for j, eta_j in enumerate(eta):
        mode_started = time.perf_counter()
        H = eta_j * H1 - H2
        H_pad = np.zeros((dim, dim), dtype=complex)
        H_pad[:M, :M] = H
        amplitude = np.linalg.norm(psi0[:, j])
        if amplitude == 0:
            tails.append(0.0); norms.append(0.0); raw_gates.append(0); expanded_gates.append(0)
            continue
        state0 = np.zeros(dim, dtype=complex)
        state0[:M] = psi0[:, j] / amplitude
        if args.pauli_action:
            state_t, gate_count = evolve_pauli_product(H_pad, state0, config["T"], config["r"])
            expanded_count = None
        else:
            qc = evolution.trotter_2nd_evolution(H_pad, config["T"], config["r"])
            # Each original gate contains one Pauli term; expansion preserves that
            # gate exactly and retains the original symmetric product ordering.
            expanded = qc.decompose()
            state_t = Statevector(state0).evolve(expanded).data
            gate_count, expanded_count = len(qc.data), len(expanded.data)
        psi_t[:, j] = amplitude * state_t[:M]
        tails.append(float(np.sum(np.abs(state_t[M:]) ** 2)))
        norms.append(float(abs(np.linalg.norm(state_t) - 1)))
        raw_gates.append(gate_count); expanded_gates.append(expanded_count)
        print(f"mode {j+1}/{len(eta)} eta={eta_j:.5g}: {time.perf_counter()-mode_started:.2f}s; padded probability={tails[-1]:.3e}", flush=True)
        np.savez_compressed(output / "mode_checkpoint.npz", completed_modes=j+1,
                            psi_t=psi_t,
                            padded_probability=np.array(tails))
    Wt = schrodingerisation.inverse_fft(psi_t)
    j0 = int(np.argmin(np.abs(p)))
    assert np.isclose(p[j0], 0)
    phi_trotter = Wt[:len(x), j0]
    # Independent nonlinear reference on the SAME spatial grid and equation.
    def nonlinear_rhs(t, phi):
        return A @ phi + config["kappa"] * phi**2
    def jacobian(t, phi):
        return A + np.diag(2 * config["kappa"] * phi)
    sol = solve_ivp(nonlinear_rhs, (0, config["T"]), phi0, method="Radau",
                    jac=jacobian, rtol=1e-10, atol=1e-12, t_eval=[config["T"]])
    if not sol.success:
        raise RuntimeError(sol.message)
    phi_nonlinear = sol.y[:, -1]
    np.savez_compressed(output / "results.npz", x=x, phi0=phi0, p=p, eta=eta,
                        phi_trotter=phi_trotter, phi_nonlinear=phi_nonlinear,
                        psi0=psi0, psi_t=psi_t, Wt=Wt,
                        padded_probability=np.array(tails))
    np.savetxt(output / "profiles.csv", np.column_stack((x, phi0, phi_nonlinear,
               phi_trotter.real, phi_trotter.imag)), delimiter=",", comments="",
               header="x,phi0,finite_difference,phit_real,phit_imag")
    summary = {
        "parameters": config, "Nx": len(x), "M": M, "padded_dimension": dim,
        "n_qubits": (dim-1).bit_length(), "deltap": dp,
        "H1_min_eigenvalue": float(spectrum[0]), "H1_max_eigenvalue": float(spectrum[-1]),
        "errors": {
            "phit_vs_finite_difference": relative_error(phi_trotter, phi_nonlinear),
        },
        "max_abs_imaginary_phi": float(np.max(np.abs(phi_trotter.imag))),
        "max_padded_probability": max(tails), "max_full_state_norm_error": max(norms),
        "raw_gate_counts": raw_gates, "expanded_gate_counts": expanded_gates,
        "runtime_seconds": time.perf_counter()-started,
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "qiskit": qiskit.__version__, "matplotlib": matplotlib.__version__},
        "source_status": source_status(),
        "method": ("Classical direct Pauli action with the same decomposition tolerance, term order, angles and repetitions as the copied Trotter circuit; no hardware run."
                   if args.pauli_action else
                   "Copied Trotter circuit, single-Pauli gates expanded with qc.decompose(), Qiskit Statevector simulation; no hardware run."),
        "execution_backend": "compiled Pauli action" if getattr(args, "compiled_pauli", False) else ("NumPy Pauli action" if args.pauli_action else "Qiskit statevector"),
        "limitations": "Single coarse-grid run, not a convergence study. p=0 recovery is approximate on the finite auxiliary grid. Plot shows real components; errors use full complex values. Padding leakage is recorded, not renormalized away.",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    make_plots(output)
    print(json.dumps(summary["errors"], indent=2), flush=True)
    print(f"Saved to {output}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=8)
    parser.add_argument("--start", type=float)
    parser.add_argument("--end", type=float)
    parser.add_argument("--mu", type=float, help="Use exp(-(x-mu)**2) as the initial profile.")
    parser.add_argument("--np", type=int)
    parser.add_argument("--r", type=int)
    parser.add_argument("--time", type=float)
    parser.add_argument("--name", default="first_run_20261016")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--pauli-action", action="store_true", help="Evaluate the same Pauli product directly on the statevector.")
    parser.add_argument("--compiled-pauli", action="store_true", help="Use the optional compiled worker (build it first).")
    args = parser.parse_args()
    if args.compiled_pauli:
        from compiled_pauli_action import evolve_pauli_product
        args.pauli_action = True
    run(args)
