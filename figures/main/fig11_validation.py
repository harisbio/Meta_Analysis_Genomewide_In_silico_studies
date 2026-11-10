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
import matplotlib.patheffects as pe
from matplotlib.patches import Circle
import numpy as np
import pandas as pd

DATA = str(DATA_FILE)
OUT = Path(OUTPUT_DIR / "Figure_11_Validation")
OUT.mkdir(parents=True, exist_ok=True)
INK = "#16324f"
plt.rcParams.update({"font.family": "Arial", "font.size": 20, "text.color": INK, "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})
RINGS = [("Only in silico", "#7d8a97"), ("RNA-seq only", "#2fa985"), ("qPCR only", "#e8a93a"), ("qPCR + RNA-seq", "#1f7fc0")]   # inner to outer

d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
d["y"] = d.Year.astype(int)
N = len(d)
q, r = d["Experimental Validation (qPCR)"] == "Yes", d["Expression Analysis (RNA-seq)"] == "Yes"
d["k"] = np.select([q & r, q, r], ["qPCR + RNA-seq", "qPCR only", "RNA-seq only"], "Only in silico")
years = list(range(d.y.min(), d.y.max() + 1))
tab = pd.crosstab(d.y, d.k).reindex(index=years, columns=[n for n, _ in RINGS], fill_value=0)

fig, ax = plt.subplots(figsize=(17, 14.8), facecolor="white")
fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax.set_aspect("equal")
ax.axis("off")
GAP = 62
step = (360 - GAP) / (len(years) - 1)
start = 90 - GAP / 2
RAD = {n: 0.42 + i * 0.2 for i, (n, _) in enumerate(RINGS)}


def pt(a, r_):
    return r_ * math.cos(math.radians(a)), r_ * math.sin(math.radians(a))


ang = {y: start - i * step for i, y in enumerate(years)}
for n, col in RINGS:
    t = np.radians(np.linspace(start, start - (360 - GAP), 300))
    ax.plot(RAD[n] * np.cos(t), RAD[n] * np.sin(t), color=col, lw=2.4, alpha=0.4, zorder=1)
for y, a in ang.items():
    x0, y0 = pt(a, 0.33)
    x1, y1 = pt(a, 1.09)
    ax.plot([x0, x1], [y0, y1], color="#e3e8ee", lw=1.3, zorder=0)
    lx, ly = pt(a, 1.14)
    ax.text(lx, ly, str(y), ha="left" if lx > 0.03 else ("right" if lx < -0.03 else "center"), va="center", fontsize=24, fontweight="bold")
vmax = tab.values.max()
for n, col in RINGS:
    for y in years:
        v = tab.at[y, n]
        if not v:
            continue
        x, y_ = pt(ang[y], RAD[n])
        ax.scatter(x, y_, s=300 + 115 * v, color=col, edgecolors="white", linewidths=1.5, zorder=3)
        if v >= 2:
            ax.text(x, y_, str(v), ha="center", va="center", fontsize=25 if v >= 10 else 21, fontweight="bold", color="white", zorder=4,
                    path_effects=[pe.withStroke(linewidth=2.5, foreground=col)])
for n, col in RINGS:
    ax.text(0, RAD[n], f"{n} ({int(tab[n].sum())})", ha="center", va="center", fontsize=22, fontweight="bold", color="#16324f", zorder=5, path_effects=[pe.withStroke(linewidth=6, foreground="white")])
ax.set_xlim(-1.38, 1.38)
ax.set_ylim(-1.23, 1.17)


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


save("Figure_11_Validation.png", dpi=300)
save("Figure_11_Validation.pdf")
save("Figure_11_Validation.svg")
out = tab.copy()
out["total"] = out.sum(axis=1)
out.to_csv(OUT / "Figure_11_data.csv", encoding="utf-8-sig")
print({k: (int(v), round(v / N * 100, 1)) for k, v in d.k.value_counts().items()}, "qPCR", int(q.sum()), "RNA", int(r.sum()))
m = d["Expression (microarray/EST)"] == "Yes"
print("microarray in each type", pd.crosstab(d.k, m).to_dict())
per = pd.cut(d.y, [1990, 2012, 2016, 2020, 2030], labels=["<2013", "13-16", "17-20", "21-25"])
print(pd.crosstab(per, d.k, normalize="index").round(3).mul(100).to_string())
