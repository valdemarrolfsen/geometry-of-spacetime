"""Regenerate every current figure, numerical table, and paper result macro."""

import csv
import importlib.metadata
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from spacetime.model import LocalizedPerturbation
from spacetime.rays import born_angular_shift, observer_frame, radial_emission_shift, shoot_to_source, trace_ray

FIGURES = ROOT / "paper/figures"
RESULTS = ROOT / "results"
GENERATED = ROOT / "paper/generated"
for directory in [FIGURES, RESULTS, GENERATED]:
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.labelsize": 10, "axes.titlesize": 11,
                     "svg.hashsalt": "geometry-of-spacetime-v0.5", "savefig.dpi": 170})
BLUE, ORANGE, GREEN, GRAY = "#2166ac", "#b35806", "#008675", "#89939f"


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    svg = FIGURES / f"{name}.svg"
    fig.savefig(svg, bbox_inches="tight", metadata={"Date": None})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    fig.savefig(FIGURES / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def csv_write(name, columns, rows):
    with (RESULTS / name).open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def main():
    model = LocalizedPerturbation()
    b = model.background
    eo = model.eta_observer
    reference = LocalizedPerturbation(epsilon=0)
    eta = np.linspace(0.35, eo, 600)
    chi = np.linspace(0, np.pi, 1201)
    mu = np.cos(chi)
    phi = model.potential(eo, mu)
    density = model.density_comoving(eo, mu)
    density_N = model.density_newtonian(eo, mu)
    clock_offset = np.array([quad(lambda e: b.A(e) * model.potential(e, m),
                                 0.4, eo, epsabs=1e-14)[0] for m in mu])
    # Integrated normal-observer clock offset, first order, relative to the
    # background cosmic-time interval. This is a specified gauge comparison.
    csv_write("profile.csv", ["chi_rad", "potential", "delta_comoving", "delta_Newtonian", "clock_offset_Amax_over_c"],
              zip(chi, phi, density, density_N, clock_offset))

    # Corresponding conformal and temporal-radial pictures of the same ray.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.1), constrained_layout=True)
    x = eo - eta
    axes[0].plot(x, eta, color=BLUE, lw=2.2)
    for et, label in [(eo, "Reception"), (1.2, "Later emission"), (0.65, "Earlier emission")]:
        xx, aa = eo-et, b.A(et)/b.A(eo)
        axes[0].axhline(et, color=GRAY, ls=":", lw=1)
        axes[0].plot(xx, et, "o", color=BLUE, ms=4)
        axes[0].annotate(label, (xx, et), xytext=(8, 6), textcoords="offset points", fontsize=8)
        theta = np.linspace(0, 1.55, 160)
        axes[1].plot(aa*np.cos(theta), aa*np.sin(theta), ":", color=GRAY)
        axes[1].plot(aa*np.cos(xx), aa*np.sin(xx), "o", color=BLUE, ms=4)
    axes[0].annotate("", xy=(0.95, eo-0.95), xytext=(0.7, eo-0.7), arrowprops={"arrowstyle": "->", "color": BLUE})
    axes[0].set(xlabel=r"Spatial coordinate $\chi$ (rad)", ylabel=r"Conformal time $\eta$", title="The past-directed null ray", xlim=(-.05, 1.8), ylim=(.3, 2))
    radius = b.A(eta)/b.A(eo)
    axes[1].plot(radius*np.cos(x), radius*np.sin(x), color=BLUE, lw=2.2)
    axes[1].plot([0, 1], [0, 0], "--", color=GRAY, lw=1)
    axes[1].text(1.02, 0, "Observer", fontsize=9)
    axes[1].text(0.01, -.085, "Earlier layers lie inward", fontsize=9)
    axes[1].set(xlabel=r"$(A/A_o)\cos\chi$", ylabel=r"$(A/A_o)\sin\chi$", title="The same ray in nested layers", xlim=(-.05, 1.24), ylim=(-.12, 1.05), aspect="equal")
    save(fig, "01-view-lines")

    # Spatial embedding deformation and potential. Magnification is explicit.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axes[0].plot(chi, phi/model.epsilon, color=BLUE, label=r"$\varphi/\epsilon$")
    axes[0].axhline(0, color=GRAY, lw=.8)
    axes[0].set(xlabel=r"Distance from mass $\chi$ (rad)", ylabel="Normalized potential", xlim=(0, np.pi), title="Localized well and global compensation")
    theta = np.linspace(-np.pi, np.pi, 1000)
    magnification = 1000
    deformation = -magnification*model.potential(eo, np.cos(theta))
    axes[1].plot(np.cos(theta), np.sin(theta), "--", color=GRAY, label="Unperturbed slice")
    axes[1].plot((1+deformation)*np.cos(theta), (1+deformation)*np.sin(theta), color=ORANGE, label="Spatial radius change ×1000")
    axes[1].plot(1.1, 0, "o", color=ORANGE, ms=4)
    axes[1].annotate("Overdensity", (1.1, 0), xytext=(.5, .48), arrowprops={"arrowstyle": "->", "color": ORANGE}, fontsize=9)
    axes[1].set(aspect="equal", title="Great-circle section of the spatial embedding", xlabel=r"$X_1/A_o$", ylabel=r"$X_2/A_o$")
    axes[1].legend(fontsize=8, loc="lower center", bbox_to_anchor=(.5, -.35), frameon=False)
    save(fig, "02-spatial-deformation")

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True)
    axes[0].plot(chi, 100*density, color=ORANGE, label="Comoving density contrast")
    axes[0].plot(chi, 100*density_N, color=BLUE, ls="--", label="Newtonian-gauge contrast")
    axes[0].axhline(0, color=GRAY, lw=.8)
    axes[0].set(xlim=(0, 1.3), xlabel=r"$\chi$ (rad)", ylabel="Density perturbation (%)", title="Small positive core; compensating surroundings")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(chi, 1e6*clock_offset, color=GREEN)
    axes[1].axhline(0, color=GRAY, lw=.8)
    axes[1].set(xlim=(0, 1.3), xlabel=r"$\chi$ (rad)", ylabel=r"Clock offset ($10^{-6} A_{\max}/c$)", title=r"Elapsed-clock change, $\eta=0.4$ to $1.8$")
    save(fig, "03-density-and-clocks")

    # Rays are in a totally geodesic S^2 containing the observer and lens.
    angles = [0.12, 0.25, 0.4]
    shootings = [shoot_to_source(model, baseline_angle=a) for a in angles]
    csv_write("matched-sources.csv", ["baseline_angle_rad", "angular_shift_rad", "additional_conformal_lookback", "source_residual", "normal_observer_redshift"],
              [(a, s["angular_shift_rad"], s["additional_conformal_lookback"], s["source_residual"], s["ray"]["redshift"][-1]) for a, s in zip(angles, shootings)])
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    colors = [BLUE, ORANGE, GREEN]
    for angle, shot, color in zip(angles, shootings, colors):
        ray = shot["ray"]
        n = ray["n"]
        radius_chi = np.arccos(np.clip(n[:, 0], -1, 1))
        azimuth = np.arctan2(n[:, 2], n[:, 1])
        axes[0].plot(radius_chi*np.cos(azimuth), radius_chi*np.sin(azimuth), color=color,
                     label=fr"$\alpha_0={angle:.2f}$ rad")
        axes[0].plot(radius_chi[-1]*np.cos(azimuth[-1]), radius_chi[-1]*np.sin(azimuth[-1]), "o", color=color, ms=4)
    axes[0].plot(0, 0, "*", color="black", ms=9)
    axes[0].plot(.65, 0, "s", color="black", ms=5)
    axes[0].text(.65, -.06, "Observer", ha="center", fontsize=8)
    axes[0].text(0, -.06, "Mass", ha="center", fontsize=8)
    axes[0].set(aspect="equal", xlabel=r"$\chi\cos\beta$ (rad)", ylabel=r"$\chi\sin\beta$ (rad)", title="Spatial projection of matched-source rays", ylim=(-.14, .46), xlim=(-.75, .82))
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(angles, [s["angular_shift_rad"]*1e6 for s in shootings], "o-", color=BLUE)
    axes[1].set(xlabel=r"Unperturbed observer angle $\alpha_0$ (rad)", ylabel="Observer-frame angular shift (µrad)", title="Deflection relative to the same source")
    save(fig, "04-matched-light-rays")

    # Keep the inward trajectory visible in a 3D section of the 4D picture.
    fig = plt.figure(figsize=(7.1, 5.0), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    theta, beta = np.meshgrid(np.linspace(0, 1.4, 18), np.linspace(-.13, .8, 14))
    for et in [.65, 1.2, eo]:
        radius = b.A(et)/b.A(eo)
        ax.plot_wireframe(radius*np.cos(theta), radius*np.sin(theta)*np.cos(beta),
                          radius*np.sin(theta)*np.sin(beta), color=GRAY, alpha=.23, linewidth=.45)
    for shot, color in zip(shootings, colors):
        ray = shot["ray"]
        coords = b.A(ray["eta"])[:, None] / b.A(eo) * ray["n"]
        ax.plot(*coords.T, color=color, lw=2)
        ax.scatter(*coords[-1], color=color, s=15)
    n0, _ = observer_frame()
    ax.scatter(*n0, color="black", marker="s", s=25)
    ax.text(*n0, "  Reception", fontsize=9)
    ax.scatter(0, 0, 0, color=GRAY, s=10)
    ax.set(xlabel=r"$X_1/A_o$", ylabel=r"$X_2/A_o$", zlabel=r"$X_3/A_o$", title="Past light paths cross successively earlier layers")
    ax.set_box_aspect((1, 1, .7))
    ax.view_init(elev=23, azim=35)
    save(fig, "05-nested-view-lines")

    delay = radial_emission_shift(model)
    amplitudes = [2.5e-5, 5e-5, 1e-4, 2e-4, 4e-4]
    convergence = []
    for amplitude in amplitudes:
        result = radial_emission_shift(LocalizedPerturbation(epsilon=amplitude))
        convergence.append([amplitude, result["additional_conformal_lookback"], result["born_conformal_lookback"],
                            abs(result["additional_conformal_lookback"] - result["born_conformal_lookback"])])
    convergence = np.array(convergence)
    csv_write("delay-convergence.csv", ["epsilon", "numerical_delay", "first_order_delay", "absolute_difference"], convergence)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), constrained_layout=True)
    axes[0].loglog(convergence[:, 0], convergence[:, 1], "o", color=BLUE, label="Ray integration")
    axes[0].loglog(convergence[:, 0], convergence[:, 2], "-", color=ORANGE, label="First-order path integral")
    axes[0].set(xlabel=r"Potential amplitude $\epsilon$", ylabel=r"Additional conformal lookback $\Delta s$", title="Timing change is linear at leading order")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].loglog(convergence[:, 0], convergence[:, 3], "o-", color=GREEN, label="Difference")
    axes[1].loglog(convergence[:, 0], convergence[0, 3]*(convergence[:, 0]/convergence[0, 0])**2, "--", color=GRAY, label=r"$\epsilon^2$ reference")
    axes[1].set(xlabel=r"$\epsilon$", ylabel="Absolute remainder", title="Retained higher-order terms are controlled")
    axes[1].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.set_xticks([2.5e-5, 1e-4, 4e-4],
                        [r"$2.5\times10^{-5}$", r"$10^{-4}$", r"$4\times10^{-4}$"])
        axis.xaxis.set_minor_formatter(NullFormatter())
    save(fig, "06-delay-convergence")

    core_radius = brentq(lambda cc: model.density_comoving(eo, np.cos(cc)), .05, .6)
    core_integral = quad(lambda cc: model.density_comoving(eo, np.cos(cc))*np.sin(cc)**2,
                         0, core_radius, epsabs=1e-14)[0]
    mass_geometric = 1.5*b.A(eo)*(b.h(eo)**2+1)*core_integral
    chosen = shootings[1]
    ray = chosen["ray"]
    n, p = ray["n"], ray["p"]
    sample_eta = np.linspace(.4, eo, 80)
    velocity_max = max(np.max(np.abs(-2*model.epsilon*(model.growth(et, True)+b.h(et)*model.growth(et))
                                    / (3*(b.h(et)**2+1))*model.profile(mu, 1)*np.sin(chi))) for et in sample_eta)
    diagnostics = {
        "model": "closed dust, Lambda=0, scalar growing mode, linear Einstein equations",
        "units": "c=G=A_max=1; angles and conformal time dimensionless",
        "parameters": {"epsilon": model.epsilon, "sigma_rad": model.sigma, "eta_observer": eo,
                       "chi_observer": .65, "baseline_source_distance_rad": 1.3},
        "potential_at_core": float(phi[0]),
        "max_abs_density_comoving": float(np.max(np.abs(density))),
        "max_abs_density_Newtonian": float(np.max(np.abs(density_N))),
        "max_peculiar_speed_over_c": float(velocity_max),
        "core_radius_rad": float(core_radius),
        "positive_core_mass_GM_over_c2Amax": float(mass_geometric),
        "core_lapse": float(model.lapse(eo, 1)),
        "core_radial_embedding_fraction_first_order": float(-phi[0]),
        "core_clock_offset_Amax_over_c": float(clock_offset[0]),
        "radial_timing": delay,
        "matched_source_025": {key: float(chosen[key]) for key in ["angular_shift_rad", "additional_conformal_lookback", "source_residual"]},
        "first_order_angular_shift_rad": float(born_angular_shift(model)),
        "max_sphere_constraint_error": float(np.max(np.abs(np.sum(n*n, axis=1)-1))),
        "max_tangency_constraint_error": float(np.max(np.abs(np.sum(n*p, axis=1)))),
        "delay_remainder_power": float(np.polyfit(np.log(convergence[:, 0]), np.log(convergence[:, 3]), 1)[0]),
        "software": {package: importlib.metadata.version(package) for package in ["numpy", "scipy", "sympy", "matplotlib"]},
    }
    (RESULTS / "diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    macros = {
        "CoreDensityPercent": f"{100*density[0]:.4f}",
        "CoreLapse": f"{model.lapse(eo, 1):.8f}",
        "CoreClockOffset": f"{clock_offset[0]*1e6:.3f}",
        "RadialDelay": f"{delay['additional_conformal_lookback']*1e6:.3f}",
        "BornDelay": f"{delay['born_conformal_lookback']*1e6:.3f}",
        "AngularShift": f"{chosen['angular_shift_rad']*1e6:.3f}",
        "BornAngularShift": f"{born_angular_shift(model)*1e6:.3f}",
        "CoreRadius": f"{core_radius:.4f}",
        "CoreMassScaled": f"{mass_geometric*1e6:.4f}",
        "RemainderPower": f"{diagnostics['delay_remainder_power']:.3f}",
        "MaxPeculiarSpeed": f"{velocity_max*1e6:.3f}",
    }
    (GENERATED / "numbers.tex").write_text("% Generated by scripts/reproduce.py; do not edit.\n" +
        "".join("\\newcommand{\\"+key+"}{"+value+"}\n" for key, value in macros.items()))
    print(json.dumps({"figures": 6, "core_density_percent": 100*density[0],
                      "angular_shift_microrad": chosen["angular_shift_rad"]*1e6,
                      "delay_remainder_power": diagnostics["delay_remainder_power"]}, indent=2))


if __name__ == "__main__":
    main()
