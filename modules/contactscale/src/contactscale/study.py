"""Reproduce the synthetic contact-resolution experiment and all paper data."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter
from scipy.integrate import quad
from scipy.linalg import expm

from .analytical import contact_flux, contact_heat, contact_scale
from .solver import Slab, ThermalChain, two_slabs

LEFT = Slab(0.06, 16, 4e6)
RIGHT = Slab(0.03, 2, 2e6)
RESISTANCE = 2e-4
DELTA = 200.0


def _write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(output):
    """Write deterministic numerical CSV/JSON, a figure, and generated LaTeX data."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    args = (DELTA, LEFT.effusivity, RIGHT.effusivity, RESISTANCE)
    times = np.geomspace(1e-3, 10, 101)
    history, refinement = [], []
    exact_flux = contact_flux(times, *args)
    exact_heat = contact_heat(times, *args)
    max_balance = 0.0
    finest = None
    for n in (20, 40, 80, 160, 320):
        model, face = two_slabs(LEFT, RIGHT, n, RESISTANCE)
        initial = np.repeat([500.0, 300.0], n)
        result = model.evolve(initial, times)
        flux = np.array([model.flux(t)[face] for t in result])
        heat = (result[:, n:] - 300) @ model.mass[n:]
        balance = (result - initial) @ model.mass
        max_balance = max(max_balance, float(np.max(np.abs(balance))))
        q_at_one = model.flux(model.evolve(initial, [1])[0])[face]
        heat_at_one = (model.evolve(initial, [1])[0][n:] - 300) @ model.mass[n:]
        r_mesh = LEFT.length / (2 * n * LEFT.conductivity) + RIGHT.length / (
            2 * n * RIGHT.conductivity
        )
        refinement.append(
            {
                "cells_per_slab": n,
                "mesh_resistance_over_contact": r_mesh / RESISTANCE,
                "initial_flux_ratio": RESISTANCE / (RESISTANCE + r_mesh),
                "flux_error_at_1s_relative": abs(q_at_one / contact_flux(1, *args) - 1),
                "heat_error_at_1s_relative": abs(heat_at_one / contact_heat(1, *args) - 1),
                "max_flux_error_relative": float(np.max(np.abs(flux / exact_flux - 1))),
                "max_heat_error_relative": float(np.max(np.abs(heat / exact_heat - 1))),
            }
        )
        for i, time in enumerate(times):
            left_t, right_t = model.surface_temperatures(result[i], face)
            history.append(
                {
                    "cells_per_slab": n,
                    "time_s": time,
                    "flux_W_m2": flux[i],
                    "reference_flux_W_m2": exact_flux[i],
                    "heat_J_m2": heat[i],
                    "reference_heat_J_m2": exact_heat[i],
                    "left_surface_K": left_t,
                    "right_surface_K": right_t,
                }
            )
        if n == 320:
            finest = flux
    # Domain error control: same h, twice the distance to insulated boundaries.
    # Use n=160 -> 320 per slab so the reference comparison is affordable.
    base, face = two_slabs(LEFT, RIGHT, 160, RESISTANCE)
    doubled, face2 = two_slabs(
        Slab(2 * LEFT.length, 16, 4e6), Slab(2 * RIGHT.length, 2, 2e6), 320, RESISTANCE
    )
    base_flux = base.flux(base.evolve(np.repeat([500.0, 300.0], 160), [10])[0])[face]
    double_flux = doubled.flux(doubled.evolve(np.repeat([500.0, 300.0], 320), [10])[0])[face2]
    smooth = []
    for n in (10, 20, 40, 80):
        h = 1 / n
        model = ThermalChain(np.full(n, h), np.ones(n), np.ones(n))
        initial = np.diff(np.sin(np.pi * np.linspace(0, 1, n + 1))) / (np.pi * h)
        error = np.linalg.norm(
            model.evolve(initial, [0.1])[0] - initial * np.exp(-(np.pi**2) * 0.1)
        ) / np.sqrt(n)
        smooth.append({"cells": n, "rms_cell_average_error": error})
    # Independent dense exponential on a small heterogeneous chain.
    model, face = two_slabs(LEFT, RIGHT, 4, RESISTANCE)
    initial = np.repeat([500.0, 300.0], 4)
    k = np.diag(model.diagonal) - np.diag(model.conductance, 1) - np.diag(model.conductance, -1)
    difference = np.max(
        np.abs(model.evolve(initial, [0.1])[0] - expm(-k / model.mass[:, None] * 0.1) @ initial)
    )
    quadrature = quad(lambda u: 2 * u * contact_flux(u * u, *args), 0, 1, epsabs=1e-7)[0]
    metrics = {
        "schema_version": 1,
        "parameters": {
            "left_length_m": 0.06,
            "right_length_m": 0.03,
            "left_k_W_mK": 16,
            "right_k_W_mK": 2,
            "left_capacity_J_m3K": 4e6,
            "right_capacity_J_m3K": 2e6,
            "initial_left_K": 500,
            "initial_right_K": 300,
            "contact_resistance_m2K_W": RESISTANCE,
            "time_min_s": 0.001,
            "time_max_s": 10,
            "time_count": 101,
        },
        "crossover_time_s": contact_scale(LEFT.effusivity, RIGHT.effusivity, RESISTANCE),
        "reference_flux_at_1s_W_m2": float(contact_flux(1, *args)),
        "reference_heat_at_1s_J_m2": float(contact_heat(1, *args)),
        "perfect_contact_heat_at_1s_J_m2": float(contact_heat(1, *args[:-1], 0)),
        "perfect_contact_heat_excess_at_1s_relative": float(
            contact_heat(1, *args[:-1], 0) / contact_heat(1, *args) - 1
        ),
        "dense_exponential_error_K": float(difference),
        "heat_quadrature_error_J_m2": float(abs(quadrature - contact_heat(1, *args))),
        "max_energy_balance_residual_J_m2": max_balance,
        "domain_doubling_flux_relative_at_10s": float(
            abs(base_flux - double_flux) / abs(double_flux)
        ),
        "contact_flux_orders_at_1s": [
            float(np.log2(a["flux_error_at_1s_relative"] / b["flux_error_at_1s_relative"]))
            for a, b in zip(refinement[:-1], refinement[1:], strict=True)
        ],
        "smooth_orders": [
            float(np.log2(a["rms_cell_average_error"] / b["rms_cell_average_error"]))
            for a, b in zip(smooth[:-1], smooth[1:], strict=True)
        ],
        "refinement": refinement,
    }
    (output / "metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n"
    )
    _write_csv(output / "history.csv", history)
    _write_csv(output / "refinement.csv", refinement)
    _write_csv(output / "smooth.csv", smooth)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), layout="constrained")
    axes[0].loglog(times, exact_flux / 1e6, color="#163b65", label="Finite contact: analytical")
    axes[0].loglog(
        times,
        contact_flux(times, *args[:-1], 0) / 1e6,
        "--",
        color="#c26732",
        label="Perfect contact",
    )
    axes[0].loglog(times, finest / 1e6, ":", color="#1c8a73", label="Finite volume: 320/slab")
    axes[0].axvline(metrics["crossover_time_s"], color="0.5", lw=0.8)
    axes[0].set(xlabel="Time (s)", ylabel="Contact flux (MW/m²)")
    axes[0].legend(fontsize=8)
    ns = [r["cells_per_slab"] for r in refinement]
    axes[1].loglog(
        ns, [1 - r["initial_flux_ratio"] for r in refinement], "o-", label="Initial flux error"
    )
    axes[1].loglog(
        ns, [r["flux_error_at_1s_relative"] for r in refinement], "s-", label="Flux error at 1 s"
    )
    axes[1].loglog(
        ns, [r["heat_error_at_1s_relative"] for r in refinement], "^-", label="Heat error at 1 s"
    )
    axes[1].set(xlabel="Cells per slab", ylabel="Relative error")
    axes[1].set_xticks(ns, [str(n) for n in ns])
    axes[1].xaxis.set_minor_formatter(NullFormatter())
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.2, which="both")
    fig.savefig(output / "contact-study.png", dpi=200)
    plt.close(fig)
    # Publication table is always generated from the same measured arrays.
    table = [
        r"\begin{tabular}{rrrr}",
        r"\toprule",
        r"Cells/slab & $q_h(0)/q(0)$ & $E_q(1)$ & $E_Q(1)$ \\",
        r"\midrule",
    ]
    for row in refinement:
        table.append(
            f"{row['cells_per_slab']} & {row['initial_flux_ratio']:.5f} & "
            f"{row['flux_error_at_1s_relative']:.3e} & "
            f"{row['heat_error_at_1s_relative']:.3e} " + r"\\"
        )
    table.extend([r"\bottomrule", r"\end{tabular}"])
    (output / "table.tex").write_text("\n".join(table) + "\n", encoding="utf-8", newline="\n")
    macros = {
        "CrossTime": metrics["crossover_time_s"],
        "ReferenceFlux": metrics["reference_flux_at_1s_W_m2"],
        "ReferenceHeat": metrics["reference_heat_at_1s_J_m2"],
        "PerfectExcess": 100 * metrics["perfect_contact_heat_excess_at_1s_relative"],
        "DenseError": difference,
        "BalanceError": max_balance,
        "DomainError": metrics["domain_doubling_flux_relative_at_10s"],
        "FinalOrder": metrics["contact_flux_orders_at_1s"][-1],
        "SmoothOrder": metrics["smooth_orders"][-1],
    }
    (output / "numbers.tex").write_text(
        "\n".join(
            "\\newcommand{\\" + name + "}{" + f"{value:.5g}" + "}" for name, value in macros.items()
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    metrics = run(args.output)
    print(
        json.dumps(
            {k: v for k, v in metrics.items() if k not in ("parameters", "refinement")}, indent=2
        )
    )


if __name__ == "__main__":
    main()
