import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import math
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import statsmodels.api as sm

DATA = str(DATA_FILE)
OUT = Path(OUTPUT_DIR / "Figure_3_Publication_Trends")
OUT.mkdir(parents=True, exist_ok=True)
INK, GRID = "#16324f", "#e3e8ee"
BLUE, RED = "#0072B2", "#D55E00"
plt.rcParams.update({"font.family": "Arial", "font.size": 20, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "text.color": INK, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.2, "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})

d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
d["y"] = d.Year.astype(int)
N = len(d)
years = np.arange(d.y.min(), d.y.max() + 1)
annual = d.y.value_counts().reindex(years, fill_value=0)

X = sm.add_constant(years - years.min())
pois = sm.GLM(annual.values, X, family=sm.families.Poisson()).fit()
alpha = max(((annual.values - pois.mu) ** 2 - pois.mu).sum() / (pois.mu ** 2).sum(), 1e-6)
nb = sm.GLM(annual.values, X, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
b, (lo, hi) = nb.params[1], nb.conf_int()[1]
grow = (math.exp(b) - 1) * 100
dbl, dbl_lo, dbl_hi = math.log(2) / b, math.log(2) / hi, math.log(2) / lo
pred = nb.get_prediction(X).summary_frame(alpha=0.05)

fig = plt.figure(figsize=(21, 7.6), facecolor="white")
gs = fig.add_gridspec(1, 2, wspace=0.14, left=0.052, right=0.995, top=0.93, bottom=0.125)
ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])

# A
ax_a.bar(years, annual.values, color=BLUE, width=0.74, zorder=2)
ax_a.fill_between(years, pred["mean_ci_lower"], pred["mean_ci_upper"], color=RED, alpha=0.18, lw=0, zorder=3)
ax_a.plot(years, pred["mean"], color=RED, lw=3.2, zorder=4)
for x in (annual.idxmax(),):
    ax_a.text(x, annual[x] + 1.5, str(annual[x]), ha="center", va="bottom", fontsize=20, fontweight="bold")
ax_a.set_ylabel("Papers per year", fontsize=22)
ax_a.set_ylim(0, annual.max() * 1.14)
ax_a.grid(axis="y", color=GRID, lw=1, zorder=0)
ax_a.text(0.03, 0.95, f"Annual growth {grow:.1f}% (95% CI {(math.exp(lo) - 1) * 100:.1f}–{(math.exp(hi) - 1) * 100:.1f}%)\nDoubling time {dbl:.1f} years (95% CI {dbl_lo:.1f}–{dbl_hi:.1f})",
          transform=ax_a.transAxes, va="top", fontsize=19, color=INK)
ax_a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=BLUE), Line2D([], [], color=RED, lw=3.2)], labels=["Observed", "Fitted trend, 95% CI"], loc="center left", bbox_to_anchor=(0.02, 0.56), frameon=False, fontsize=18)

# B
cum = annual.cumsum()
ax_b.fill_between(years, cum.values, color=BLUE, alpha=0.22, lw=0)
ax_b.plot(years, cum.values, color=BLUE, lw=3.4, marker="o", ms=8, mfc="white", mew=2.2)
share_recent = annual.loc[2022:2025].sum() / N * 100
ax_b.axvspan(2021.5, years.max() + 0.5, color="#f3e3d8", lw=0, zorder=0)
ax_b.text(2021.3, cum.max() * 0.74, f"{share_recent:.1f}% of all papers\nin 2022–2025", ha="right", va="center", fontsize=19)
ax_b.set_ylabel("Cumulative papers", fontsize=22)
ax_b.set_xlabel("Publication year", fontsize=22)
ax_a.set_xlabel("Publication year", fontsize=22)
ax_b.set_ylim(0, N * 1.14)
ax_b.grid(axis="y", color=GRID, lw=1)
for ax in (ax_a, ax_b):
    ax.set_xlim(years.min() - 0.7, years.max() + 0.7)
    ax.set_xticks(range(2006, 2026, 3))
    ax.tick_params(labelsize=19, length=6)

for ax, letter, xo in ((ax_a, "A", -0.10), (ax_b, "B", -0.10)):
    ax.text(xo, 1.0, letter, transform=ax.transAxes, fontsize=32, fontweight="bold", va="bottom", ha="left")


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


save("Figure_3_Publication_Trends_Collaboration.png", dpi=400)
save("Figure_3_Publication_Trends_Collaboration.pdf")
save("Figure_3_Publication_Trends_Collaboration.svg")
pd.DataFrame({"papers": annual, "cumulative": cum, "fit": pred["mean"].values, "fit_lo": pred["mean_ci_lower"].values, "fit_hi": pred["mean_ci_upper"].values}).round(2).to_csv(OUT / "Figure_3_annual_data.csv", encoding="utf-8-sig")
print("alpha", round(alpha, 3), "growth", round(grow, 1), round((math.exp(lo) - 1) * 100, 1), round((math.exp(hi) - 1) * 100, 1), "doubling", round(dbl, 2), round(dbl_lo, 2), round(dbl_hi, 2), "recent", round(share_recent, 1), "2025 share", round(annual[2025] / N * 100, 1))
