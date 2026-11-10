import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import os
import re
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

DATA = str(DATA_FILE)
OUT = Path(OUTPUT_DIR / "Figure_8_Genes_Over_Time")
OUT.mkdir(parents=True, exist_ok=True)
INK, GRID = "#16324f", "#e3e8ee"
SINGLE, MULTI = "#0072B2", "#D55E00"
PERIODS = ["Before 2013", "2013–2016", "2017–2020", "2021–2025"]
EDGES = [2012.5, 2016.5, 2020.5]
plt.rcParams.update({"font.family": "Arial", "font.size": 20, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.2, "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})

d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
d["y"] = d.Year.astype(int)
d["g"] = pd.to_numeric(d["Number of Genes Identified"], errors="coerce")
d = d[d.g > 0].copy()
d["multi"] = d["Primary Organism"].str.contains(",")
d["period"] = pd.cut(d.y, [1990, 2012, 2016, 2020, 2030], labels=PERIODS)
TICKS = [1, 10, 100, 1000]


def trend(sub):
    m = sm.OLS(np.log10(sub.g), sm.add_constant(sub.y - sub.y.min())).fit()
    b, (lo, hi) = m.params.iloc[1], m.conf_int().iloc[1]
    pct = lambda v: (10 ** v - 1) * 100
    return pct(b), pct(lo), pct(hi), m.pvalues.iloc[1], m.rsquared


t_all, t_s = trend(d), trend(d[~d.multi])
rho = stats.spearmanr(d.y, d.g)

fig = plt.figure(figsize=(20, 15.5), facecolor="white")
gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 1], hspace=0.2, left=0.075, right=0.985, top=0.965, bottom=0.065)
ax_a, ax_b = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

# A
rng = np.random.default_rng(42)
for flag, col in ((False, SINGLE), (True, MULTI)):
    s = d[d.multi == flag]
    ax_a.scatter(s.y + rng.uniform(-0.3, 0.3, len(s)), s.g, s=120, color=col, alpha=0.72, edgecolors="white", linewidths=0.8, zorder=3)
med = d.groupby("y").g.median()
ax_a.plot(med.index, med.values, color=INK, lw=3.6, marker="o", ms=8, zorder=4)
for e in EDGES:
    ax_a.axvline(e, color="#aab6c3", ls="--", lw=1.8, zorder=1)
ax_a.set_yscale("log")
ax_a.set_ylim(1.5, 4200)
ax_a.set_xlim(d.y.min() - 0.8, d.y.max() + 0.8)
ax_a.set_xticks(range(2006, 2026, 2))
ax_a.yaxis.set_major_locator(FixedLocator(TICKS[1:] + [1000]))
ax_a.yaxis.set_major_formatter(FixedFormatter(["10", "100", "1000"]))
ax_a.yaxis.set_minor_locator(NullLocator())
ax_a.grid(axis="y", color=GRID, lw=1.1, zorder=0)
ax_a.tick_params(labelsize=20, length=6)
ax_a.set_ylabel("Genes identified per paper (log scale)", fontsize=22)
ax_a.set_xlabel("Publication year", fontsize=22)
for lab, x in zip(PERIODS, [(d.y.min() - 0.8 + 2012.5) / 2, 2014.5, 2018.5, 2023]):
    ax_a.text(x, 3200, lab, ha="center", va="center", fontsize=21, fontweight="bold")
ax_a.legend(handles=[Line2D([], [], marker="o", ls="", ms=14, color=SINGLE, label="One species study"), Line2D([], [], marker="o", ls="", ms=14, color=MULTI, label="Several species study"),
                     Line2D([], [], color=INK, lw=3.6, marker="o", ms=8, label="Median per year")], loc="lower left", ncol=3, frameon=False, fontsize=20, columnspacing=1.5, bbox_to_anchor=(0.0, -0.01))

# B
pos, data, colors, rows = [], [], [], []
for i, p in enumerate(PERIODS):
    for off, flag, col in ((-0.2, False, SINGLE), (0.2, True, MULTI)):
        v = d[(d.period == p) & (d.multi == flag)].g
        if len(v):
            pos.append(i + off)
            data.append(v.values)
            colors.append(col)
            ax_b.text(i + off, 3100, f"n = {len(v)}", ha="center", va="center", fontsize=19)
            ax_b.text(i + off + (0.165 if flag else -0.165), np.median(v), f"{np.median(v):g}", ha="left" if flag else "right", va="center", fontsize=20, fontweight="bold", color=col)
            if len(v) == 1:
                ax_b.scatter([i + off], v.values, s=150, color=col, zorder=4)
    a = d[d.period == p].g
    rows.append(dict(period=p, papers=len(a), mean=round(a.mean(), 1), median=a.median(), over_100=int((a > 100).sum()), over_100_percent=round((a > 100).mean() * 100, 1), maximum=int(a.max()),
                     single_n=int((d[d.period == p].multi == False).sum()), single_median=d[(d.period == p) & (~d.multi)].g.median(), multi_n=int(d[d.period == p].multi.sum()), multi_median=d[(d.period == p) & d.multi].g.median()))
bp = ax_b.boxplot(data, positions=pos, widths=0.3, patch_artist=True, showfliers=True, flierprops=dict(marker="o", markersize=7, markerfacecolor="none", markeredgecolor="#6b7785"),
                  medianprops=dict(color="white", lw=3.4), whiskerprops=dict(color=INK, lw=1.8), capprops=dict(color=INK, lw=1.8), boxprops=dict(linewidth=0))
for patch, col in zip(bp["boxes"], colors):
    patch.set_facecolor(col)
ax_b.set_yscale("log")
ax_b.set_ylim(1.5, 4200)
ax_b.yaxis.set_major_locator(FixedLocator(TICKS[1:]))
ax_b.yaxis.set_major_formatter(FixedFormatter(["10", "100", "1000"]))
ax_b.yaxis.set_minor_locator(NullLocator())
ax_b.set_xticks(range(4), PERIODS)
ax_b.set_xlim(-0.6, 3.6)
ax_b.grid(axis="y", color=GRID, lw=1.1, zorder=0)
ax_b.tick_params(labelsize=20, length=6)
ax_b.set_ylabel("Genes identified per paper (log scale)", fontsize=22)
ax_b.set_xlabel("Publication period", fontsize=22)

fig.text(0.006, 0.975, "A", fontsize=34, fontweight="bold", va="center")
fig.text(0.006, ax_b.get_position().y1 + 0.012, "B", fontsize=34, fontweight="bold", va="center")


def save(name, **kw):
    tmp = OUT / (name + ".tmp" + Path(name).suffix)
    fig.savefig(tmp, **kw)
    for _ in range(8):
        try:
            os.replace(tmp, OUT / name)
            return
        except OSError:
            time.sleep(5)
    raise OSError(f"{name} is locked; new version left as {tmp.name}")


save("Figure_8_Genes_Over_Time.png", dpi=300)
save("Figure_8_Genes_Over_Time.pdf")
save("Figure_8_Genes_Over_Time.svg")
t = pd.DataFrame(rows)
t.to_csv(OUT / "Figure_8_period_statistics.csv", index=False, encoding="utf-8-sig")
print(t.to_string())
print("all: %/yr", [round(x, 3) for x in t_all], "| single only", [round(x, 3) for x in t_s], "| spearman", round(rho[0], 3), rho[1], "| multi papers", int(d.multi.sum()), "median all", d.g.median(), "multi vs single", d[d.multi].g.median(), d[~d.multi].g.median())
print("over 100 by group", int((d[d.multi].g > 100).sum()), int((d[~d.multi].g > 100).sum()))
