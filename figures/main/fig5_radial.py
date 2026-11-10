import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import os
import textwrap
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mc
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

DATA = str(DATA_FILE)
OUT = Path(OUTPUT_DIR / "Figure_5_Journals_Methods")
OUT.mkdir(parents=True, exist_ok=True)
INK = "#16324f"
TYPES = {"In silico + RNA-seq + qPCR": "#0072B2", "In silico + RNA-seq": "#009E73", "In silico + qPCR": "#E69F00", "Only in silico": "#6b7785"}
plt.rcParams.update({"font.family": "Arial", "text.color": INK, "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})

d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
q, r = d["Experimental Validation (qPCR)"] == "Yes", d["Expression Analysis (RNA-seq)"] == "Yes"
d["Study_Type"] = np.select([q & r, r, q], list(TYPES)[:3], list(TYPES)[3])
vc = d.Journal.value_counts()
top = sorted(vc.index, key=lambda j: (-vc[j], j))[:12]
counts = vc.reindex(top)
comp = pd.crosstab(d.Journal, d.Study_Type).reindex(index=top, columns=list(TYPES), fill_value=0)
share = comp.div(comp.sum(axis=1), axis=0)

n = len(top)
theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
width = 2 * np.pi / n * 0.9
fig = plt.figure(figsize=(15.5, 16.5), facecolor="white")
ax = fig.add_axes([0.15, 0.17, 0.70, 0.70], polar=True)
ax.set_facecolor("white")
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi / 2)
ax.set_ylim(0, 1.0)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.spines["polar"].set_visible(False)

INNER, THICK, OUTER, OUTER_MAX = 0.22, 0.30, 0.545, 0.40
ramp = mc.LinearSegmentedColormap.from_list("bars", ["#a58bd6", "#3f1f7a"])
for i, j in enumerate(top):
    a = theta[i]
    start = INNER
    for cat, col in TYPES.items():
        h = THICK * share.loc[j, cat]
        if h > 0:
            ax.bar(a, h, width=width * 0.9, bottom=start, color=col, edgecolor="white", linewidth=1.4, zorder=3)
            start += h
    bar = OUTER_MAX * counts[j] / counts.max()
    ax.bar(a, bar, width=width * 0.9, bottom=OUTER, color=ramp((counts[j] - counts.min()) / (counts.max() - counts.min())), edgecolor="white", linewidth=1.4, zorder=4)
    sx, sy = np.sin(a), np.cos(a)   # position around the circle: x to the right, y up
    ha = "center" if abs(sx) < 0.2 else ("left" if sx > 0 else "right")
    va = "bottom" if sy > 0.75 else ("top" if sy < -0.75 else "center")
    ax.text(a, OUTER + bar + 0.03, textwrap.fill(f"{j} ({counts[j]})", 24 if j.startswith("Molecular") else 25), ha=ha, va=va, fontsize=23, fontweight="bold", rotation=0, zorder=6, linespacing=1.08, multialignment=ha)
fig.legend(handles=[Patch(color=c, label=k) for k, c in TYPES.items()], loc="lower center", ncol=2, frameon=False, fontsize=24, bbox_to_anchor=(0.5, 0.135), handlelength=1.2, columnspacing=1.6)


def save(name, **kw):
    tmp = OUT / (name + ".tmp" + Path(name).suffix)
    fig.savefig(tmp, bbox_inches="tight", pad_inches=0.25, **kw)
    for _ in range(8):
        try:
            os.replace(tmp, OUT / name)
            return
        except OSError:
            time.sleep(5)
    raise OSError(f"{name} is locked; new version left as {tmp.name}")


save("Figure_5_Journals_Radial.png", dpi=400)
save("Figure_5_Journals_Radial.pdf")
save("Figure_5_Journals_Radial.svg")
comp.assign(total=counts).to_csv(OUT / "Figure_5_radial_data.csv", encoding="utf-8-sig")
print(counts.to_dict())
