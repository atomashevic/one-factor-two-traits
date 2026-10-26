"""Figures 1 to 5 as vector PDFs; default output is build/figures/ for review.

Figures 1 and 2 use population CSV exports. Figures 3 to 5 use archived sample
CSV exports from make_data.R --samples-only. Select figures and an explicit output
directory to generate review previews without changing the manuscript.
"""
import argparse
import csv
import math
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.environ.get("FIG_OUT", os.path.join(HERE, "..", "build", "figures"))
parser = argparse.ArgumentParser(description="Plot archived manuscript figure exports.")
parser.add_argument("--figures", nargs="+", type=int, choices=range(1, 6), default=list(range(1, 6)))
parser.add_argument("--output-dir", default=DEFAULT_OUT)
parser.add_argument("--suffix", default="")
options = parser.parse_args()
OUT = os.path.abspath(options.output_dir)
os.makedirs(OUT, exist_ok=True)


def save(fig, number):
    destination = os.path.join(OUT, f"fig{number}{options.suffix}.pdf")
    fig.savefig(destination)
    plt.close(fig)
    print("Wrote", destination)


def number(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite {key}: {row}")
    return value


def points(rows, key):
    rows = sorted(rows, key=lambda r: number(r, "phi"))
    if not rows or len({r["phi"] for r in rows}) != len(rows):
        raise ValueError("Missing or duplicate curve coordinates")
    return [number(r, "phi") for r in rows], [number(r, key) for r in rows]


plt.rcParams.update({
    "font.family": "serif", "font.serif": ["cmr10"], "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True,
    "mathtext.fontset": "cm", "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.6, "lines.linewidth": 1.4,
    "pdf.fonttype": 42, "legend.frameon": False,
})
BLACK, DARK, MID, LIGHT = "0.0", "0.3", "0.5", "0.7"


def dec(x, nd=1):
    s = f"{x:.{nd}f}"
    return s[1:] if s.startswith("0.") else s


def fmt_axis(ax, xticks, yticks, ypct=False):
    ax.set_xticks(xticks)
    ax.set_xticklabels(["0" if t == 0 else dec(t) for t in xticks])
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{round(t * 100)}%" for t in yticks] if ypct else
                       ["0" if t == 0 else ("1" if t == 1 else dec(t)) for t in yticks])
    ax.grid(axis="y", color="0.9", linewidth=0.5)
    ax.set_axisbelow(True)


def read(name):
    with open(os.path.join(HERE, name)) as fh:
        return list(csv.DictReader(fh))



if 1 in options.figures:
    # ------------------------------------------------------------------ Figure 1
    rows = read("fig1.csv")
    phi = [float(r["phi"]) for r in rows]
    fig, ax = plt.subplots(figsize=(5.6, 3.1))
    ax.plot(phi, [float(r["uls_A"]) for r in rows], color=MID, ls="--", lw=1.2)
    ax.plot(phi, [float(r["uls_B"]) for r in rows], color=MID, ls="--", lw=1.2, label="ULS")
    ax.plot(phi, [float(r["ml_A"]) for r in rows], color=BLACK, label="ML")
    ax.plot(phi, [float(r["ml_B"]) for r in rows], color=BLACK)
    for x, lab, y in [(0.4089, r"ML: $\varphi^*$ = .409", 0.79), (0.125, r"ULS: 1/$p$ = .125", 0.06)]:
        ax.axvline(x, color=LIGHT, ls=":", lw=1)
        ax.text(x + 0.01, y, lab, fontsize=8, color=DARK)
    ax.text(0.03, 0.735, "ML, favoured trait", fontsize=8)
    ax.text(0.17, 0.17, "ML, other trait", fontsize=8)
    ax.text(0.16, 0.465, "ULS, composite", fontsize=8, color=DARK)
    ax.text(0.62, 0.56, "Composite, both estimators", fontsize=8)
    ax.set_xlim(0, 0.95); ax.set_ylim(0, 0.82)
    fmt_axis(ax, [0, 0.2, 0.4, 0.6, 0.8], [0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xlabel(r"Trait correlation $\varphi$"); ax.set_ylabel("Block loading")
    ax.legend(loc="lower right")
    fig.tight_layout(); save(fig, 1)


if 2 in options.figures:

    rows = read("fig2.csv")
    style = {"eq": dict(color=MID, ls="-", lw=1.4, label="Equal loadings .70"),
             "mir": dict(color=MID, ls="--", lw=1.4, label="Mirrored loadings .60 to .80"),
             "mild": dict(color=BLACK, ls="-", lw=1.0, label="Unequal blocks, .68 vs .72"),
             "strong": dict(color=BLACK, ls="-", lw=2.0, label="Unequal blocks, .60 vs .80")}
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    for d, st in style.items():
        pts = [(float(r["phi"]), float(r["share"])) for r in rows if r["design"] == d]
        ax.plot(*zip(*pts), **st)
    loc = [(float(r["phi"]), float(r["local"])) for r in rows if r["design"] == "mild" and r["local"] not in ("", "NA")]
    loc.append((0.2415, 0.7127))                     # the fold (Proposition 1)
    ax.plot(*zip(*loc), color=BLACK, ls=":", lw=1.2, label="Local minimum, .68 vs .72")
    ax.plot([0.2415], [0.7127], marker="o", mfc="white", mec=BLACK, ms=4)
    ax.text(0.25, 0.70, r"fold, $\varphi$ = .24", fontsize=8)
    for x, lab in [(0.4089, r"$\varphi^* = .409$"), (0.4411, ".441")]:
        ax.plot([x], [0.5], marker="o", color=DARK, ms=3.5)
    ax.text(0.36, 0.475, r"$\varphi^*$ = .409", fontsize=8, color=DARK, ha="right")
    ax.text(0.455, 0.515, ".441", fontsize=8, color=DARK)
    ax.set_xlim(0, 0.9); ax.set_ylim(0.45, 1.02)
    fmt_axis(ax, [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], [0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel(r"Trait correlation $\varphi$"); ax.set_ylabel("Share of the followed trait")
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.13, 1, 1)); save(fig, 2)



if 3 in options.figures:
    rows = read("fig3.csv")
    if len(rows) != 64:
        raise ValueError("Figure 3 requires 64 source rows")
    vstyle = {"ml": dict(color=BLACK, ls="-", marker="o", label="ML, continuous items"),
              "rob": dict(color=MID, ls=":", marker="s", label="WLSMV, robust"),
              "sc": dict(color=DARK, ls="--", marker="^", label="WLSMV, scaled"),
              "naive": dict(color=LIGHT, ls="-.", marker="v", label="WLSMV, unscaled")}
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.2), sharey=True)
    for ax, (idx, title) in zip(axes, [("cfi", r"CFI $\geq .95$"), ("rmsea", r"RMSEA $\leq .06$")]):
        trials = [("INCPAD", 0.653), ("ESCAPE", 0.865), ("SCOPE", 0.912)]
        for t, x in trials:
            ax.axvline(x, color=LIGHT, ls=":", lw=1)
        top = ax.secondary_xaxis("top")
        top.set_xticks([x for _, x in trials])
        top.set_xticklabels([t for t, _ in trials], rotation=90, fontsize=7, color=DARK)
        top.tick_params(length=2)
        for variant, st in vstyle.items():
            selected = [r for r in rows if r["metric"] == idx and r["variant"] == variant]
            if len(selected) != 8:
                raise ValueError("Figure 3 curve requires eight rows")
            for r in selected:
                if int(r["n"]) <= 0 or abs(int(r["n_pass"]) / int(r["n"]) - number(r, "p_accept")) > 1e-12:
                    raise ValueError("Figure 3 numerator/denominator mismatch")
            ax.plot(*points(selected, "p_accept"), ms=3.5, lw=1.2, **st)
        ax.text(0.03, 0.97, title, transform=ax.transAxes, va="top", fontsize=9)
        ax.set_xlim(0.48, 0.93); ax.set_ylim(-0.02, 1.03)
        fmt_axis(ax, [0.5, 0.6, 0.7, 0.8, 0.9], [0, 0.25, 0.5, 0.75, 1], ypct=True)
        ax.set_xlabel(r"Trait correlation $\varphi$")
    axes[0].set_ylabel("Samples meeting the cutoff")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.12, 1, 1)); save(fig, 3)

if 4 in options.figures:
    rows = read("fig4.csv")
    if len(rows) != 36:
        raise ValueError("Figure 4 requires 36 source rows")
    styles = {"ML": dict(color=BLACK, ls="-", marker="o"),
              "ULS": dict(color=MID, ls="--", marker="s"),
              "WLSMV": dict(color=MID, ls=":", marker="^")}
    fig, ax = plt.subplots(figsize=(5.6, 3.0))
    for field, label in [("ml_threshold", r"ML: $\varphi^*$"), ("uls_threshold", r"ULS: $1/p$")]:
        values = {number(r, field) for r in rows}
        if len(values) != 1:
            raise ValueError("Figure 4 threshold is not constant")
        x = values.pop()
        ax.axvline(x, color=LIGHT, ls=":", lw=1)
        ax.text(x + 0.008, 0.88, label + " = " + dec(x, 3), fontsize=8, color=DARK)
    for estimator, st in styles.items():
        selected = [r for r in rows if r["est"] == estimator]
        if len(selected) != 12:
            raise ValueError("Figure 4 curve requires twelve rows")
        if any(int(r["n_C_magnitude_valid"]) <= 0 for r in selected):
            raise ValueError("Figure 4 requires valid contrast observations")
        ax.plot(*points(selected, "mean_abs_C_magnitude"), ms=3.5, lw=1.3,
                label=estimator, **st)
    ax.set_xlim(-0.01, 0.81); ax.set_ylim(-0.02, 1.03)
    fmt_axis(ax, [0, .1, .2, .3, .4, .5, .6, .7, .8], [0, .25, .5, .75, 1])
    # Quarter ticks need two decimals to avoid rounding .25 to .2.
    ax.set_yticklabels(["0", ".25", ".50", ".75", "1"])
    ax.set_xlabel(r"Trait correlation $\varphi$")
    ax.set_ylabel(r"Mean absolute magnitude contrast $|C_{\mathrm{mag}}|$")
    ax.legend(loc="center right")
    fig.tight_layout(); save(fig, 4)

if 5 in options.figures:
    rows = read("fig5_left.csv")
    if len(rows) != 24:
        raise ValueError("Figure 5 left requires 24 source rows")
    styles = [("B_l68_l72", "Loadings .68 / .72", BLACK),
              ("C_9_7_l75", "PHQ-ADS, 7 items first", MID),
              ("B_l60_l80", "Loadings .60 / .80", LIGHT)]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.0), gridspec_kw={"width_ratios": [1.5, 1]})
    for condition, label, color in styles:
        selected = [r for r in rows if r["cond"] == condition]
        for r in selected:
            den = int(r["n_flag"]); num = int(r["n_higher_objective_weak"])
            if den <= 0 or abs(num / den - number(r, "p_default_on_weak_listed_first")) > 1e-12:
                raise ValueError("Figure 5 left numerator/denominator mismatch")
        a1.plot(*points(selected, "p_default_on_weak_listed_first"), color=color,
                marker="o", ms=3, lw=1.6, label=label)
        # Only Design B fold annotations have recorded source provenance.
        if condition.startswith("B_"):
            folds = {number(r, "phi_fold") for r in selected}
            if len(folds) != 1:
                raise ValueError("Ambiguous population fold")
            fold = folds.pop()
            a1.plot([fold, fold], [-0.06, -0.02], color=color, lw=2, clip_on=False)
    a1.set_xlim(-0.01, 0.51); a1.set_ylim(-0.02, 1.03)
    fmt_axis(a1, [0, .1, .2, .3, .4, .5], [0, .25, .5, .75, 1], ypct=True)
    a1.set_xlabel(r"Trait correlation $\varphi$")
    a1.set_title("Higher objective, weaker block first", fontsize=9, loc="left")
    a1.legend(loc="center right", bbox_to_anchor=(1.0, 0.45), fontsize=7.5)
    rows = read("fig5_right.csv")
    estimators = ["ML", "MLs", "ULS", "WLSMV"]
    if len(rows) != 4 or {r["est"] for r in rows} != set(estimators):
        raise ValueError("Figure 5 right requires four estimator rows")
    values = []
    for estimator in estimators:
        r = next(r for r in rows if r["est"] == estimator)
        den = int(r["n_selected"]); num = int(r["n_first"])
        if den <= 0 or den != num + int(r["n_second"]) or abs(num / den - number(r, "p_first")) > 1e-12:
            raise ValueError("Figure 5 right numerator/denominator mismatch")
        values.append(number(r, "p_first"))
    labels = ["ML", "ML on\nscores", "ULS", "WLSMV"]
    a2.bar(range(4), values, color=[BLACK, DARK, MID, LIGHT], width=0.62)
    for i, value in enumerate(values):
        a2.text(i, value + .02, f"{round(value * 100)}%", ha="center", fontsize=8)
    a2.axhline(.5, color=MID, ls="--", lw=.8)
    a2.set_xticks(range(4)); a2.set_xticklabels(labels, fontsize=8)
    a2.set_ylim(0, 1.03); a2.set_yticks([0, .25, .5, .75, 1])
    a2.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    a2.grid(axis="y", color="0.9", linewidth=.5); a2.set_axisbelow(True)
    a2.set_title("Larger magnitude on first block", fontsize=9, loc="left")
    fig.tight_layout(); save(fig, 5)
