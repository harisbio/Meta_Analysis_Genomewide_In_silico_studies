import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import math
import os
import re
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

DATA = str(DATA_FILE)
OUT = Path(OUTPUT_DIR / "Figure_4_Team_Size")
OUT.mkdir(parents=True, exist_ok=True)
INK, GRID = "#16324f", "#e3e8ee"
RED = "#D55E00"
CLASSES = [("1–3 authors", 1, 3, "#7fb8e6"), ("4–6 authors", 4, 6, "#3a8fd0"), ("7–10 authors", 7, 10, "#0a5fa3"), (">10 authors", 11, 999, "#06305c")]
plt.rcParams.update({"font.family": "Arial", "font.size": 20, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.2, "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})

d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
d["y"] = d.Year.astype(int)
d["n"] = d.Authors.map(lambda s: len([a for a in re.split(r",|;", s) if a.strip()]))
N = len(d)
assert d.n.min() >= 1

X = sm.add_constant(d.y - d.y.min())
fit = sm.GLM(d.n, X, family=sm.families.Poisson()).fit(cov_type="HC1")
b, (lo, hi) = fit.params.iloc[1], fit.conf_int().iloc[1]
p = fit.pvalues.iloc[1]
rho, rho_p = stats.spearmanr(d.y, d.n)
ygrid = np.arange(d.y.min(), d.y.max() + 1)
pred = fit.get_prediction(sm.add_constant(ygrid - d.y.min())).summary_frame(alpha=0.05)

fig = plt.figure(figsize=(21, 8.2), facecolor="white")
gs = fig.add_gridspec(1, 2, wspace=0.14, left=0.052, right=0.995, top=0.93, bottom=0.115)
ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])

# A
counts = d.n.value_counts().reindex(range(1, d.n.max() + 1), fill_value=0)
col = {k: next(c for _, l, h, c in CLASSES if l <= k <= h) for k in counts.index}
ax_a.bar(counts.index, counts.values, color=[col[k] for k in counts.index], width=0.8, zorder=2)
for k, v in counts.items():
    if v:
        ax_a.text(k, v + counts.max() * 0.012, str(v), ha="center", va="bottom", fontsize=17)
ax_a.set_xlabel("Authors per paper", fontsize=22)
ax_a.set_ylabel("Number of papers", fontsize=22)
ax_a.set_xticks(range(1, d.n.max() + 1))
ax_a.set_xlim(0.3, d.n.max() + 0.7)
ax_a.set_ylim(0, counts.max() * 1.12)
ax_a.grid(axis="y", color=GRID, lw=1, zorder=0)
ax_a.tick_params(labelsize=18, length=6)
handles, table = [], []
for name, l, h, c in CLASSES:
    n = int(((d.n >= l) & (d.n <= h)).sum())
    handles.append(Patch(color=c, label=f"{name}: {n} ({n / N * 100:.1f}%)"))
    table.append(dict(team_class=name, papers=n, percent=round(n / N * 100, 1)))
ax_a.legend(handles=handles, loc="upper right", frameon=False, fontsize=19, handlelength=1.2, borderaxespad=0.2)
ax_a.text(0.985, 0.60, f"Mean {d.n.mean():.2f}, median {d.n.median():.0f}, range {d.n.min()}–{d.n.max()}", transform=ax_a.transAxes, ha="right", va="center", fontsize=19)

# B
rng = np.random.default_rng(11)
ax_b.scatter(d.y + rng.uniform(-0.28, 0.28, N), d.n + rng.uniform(-0.2, 0.2, N), s=62, color="#0a5fa3", alpha=0.55, edgecolors="none", zorder=2)
ax_b.fill_between(ygrid, pred["mean_ci_lower"], pred["mean_ci_upper"], color=RED, alpha=0.2, lw=0, zorder=3)
ax_b.plot(ygrid, pred["mean"], color=RED, lw=3.4, zorder=4)
ax_b.set_xlabel("Publication year", fontsize=22)
ax_b.set_ylabel("Authors per paper", fontsize=22)
ax_b.set_xlim(d.y.min() - 0.7, d.y.max() + 0.7)
ax_b.set_ylim(0, d.n.max() * 1.3)
ax_b.set_xticks(range(2006, 2026, 3))
ax_b.set_yticks(range(0, d.n.max() + 4, 3))
ax_b.grid(axis="y", color=GRID, lw=1, zorder=0)
ax_b.tick_params(labelsize=18, length=6)
pct = (math.exp(b) - 1) * 100
ptxt = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
ax_b.text(0.03, 0.97, f"Change per year {pct:+.1f}% (95% CI {(math.exp(lo) - 1) * 100:+.1f} to {(math.exp(hi) - 1) * 100:+.1f}%), {ptxt}", transform=ax_b.transAxes, va="top", fontsize=19)
ax_b.legend(handles=[Line2D([], [], marker="o", ls="", ms=10, color="#0a5fa3", alpha=0.65), Line2D([], [], color=RED, lw=3.4)], labels=["Individual papers", "Fitted trend, 95% CI"], loc="upper left", bbox_to_anchor=(0.01, 0.9), frameon=False, fontsize=19)

for ax, letter in ((ax_a, "A"), (ax_b, "B")):
    ax.text(-0.10, 1.0, letter, transform=ax.transAxes, fontsize=32, fontweight="bold", va="bottom", ha="left")


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


save("Figure_4_Team_Size.png", dpi=400)
save("Figure_4_Team_Size.pdf")
save("Figure_4_Team_Size.svg")
pd.DataFrame(table).to_csv(OUT / "Figure_4_team_classes.csv", index=False, encoding="utf-8-sig")
counts.rename("papers").rename_axis("authors_per_paper").to_csv(OUT / "Figure_4_distribution.csv", encoding="utf-8-sig")
print("mean", round(d.n.mean(), 2), "median", d.n.median(), "range", d.n.min(), d.n.max(), "classes", table)
print("trend %/yr", round(pct, 2), round((math.exp(lo) - 1) * 100, 2), round((math.exp(hi) - 1) * 100, 2), "p", p, "spearman", round(rho, 3), rho_p, "fit first/last", round(pred["mean"].iloc[0], 2), round(pred["mean"].iloc[-1], 2))
print("max rows", d[d.n >= 15][["Paper ID", "n"]].values.tolist(), "single author", d[d.n == 1]["Paper ID"].tolist())
