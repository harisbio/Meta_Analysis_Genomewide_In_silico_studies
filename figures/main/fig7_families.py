import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import collections
import os
import re
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pandas as pd

DATA = str(DATA_FILE)
OUT = Path(OUTPUT_DIR / "Figure_7_Gene_Families")
OUT.mkdir(parents=True, exist_ok=True)
STYLE = sys.argv[1] if len(sys.argv) > 1 else "final"
if STYLE != "final":
    OUT = OUT / "options"
    OUT.mkdir(exist_ok=True)
INK, GRID, BLUE = "#16324f", "#c3cfdb", "#0a5fa3"
MIN_PAPERS = 5
plt.rcParams.update({"font.family": "Arial", "font.size": 20, "text.color": INK, "axes.edgecolor": INK, "xtick.color": INK, "ytick.color": INK, "axes.labelcolor": INK, "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})

d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
d["y"] = d.Year.astype(int)
N = len(d)


def used(v):
    return v.strip().lower() not in {"", "no", "none", "false", "n/a", "na", "not used", "not specified", "not available"}


CLASSES = {"Domain": ["HMMER (with Pfam domains ID)", "BLAST (E-value cut-off used)", "SMART", "CDD"], "Motifs": ["MEME", "Others (Tools for Motif Analysis)"],
           "Phylogeny": ["MEGA", "IQ-TREE", "Phylogenetic Structure"], "Synteny": ["MCScanX", "Other (Tools for Synteny)"], "Validation": ["Experimental Validation (qPCR)", "Expression Analysis (RNA-seq)"]}
for k, cols in CLASSES.items():
    d[k] = d[cols].apply(lambda r: any(used(v) for v in r), axis=1)
MAPK_SET = {"MAPK", "MAPKK", "MAPKKK"}
SHORT = {"SPL (SBP-box)": "SPL"}


def families(v):
    out = [x.strip() for x in re.split(r",(?![^(]*\))", v) if x.strip()]
    out = ["MAPK cascade" if x in MAPK_SET else SHORT.get(x, x) for x in out]
    return list(dict.fromkeys(out))


d["fams"] = d["Gene Family"].map(families)
cnt = collections.Counter(f for x in d.fams for f in x)
fams = sorted([f for f, v in cnt.items() if v >= MIN_PAPERS], key=lambda f: (-cnt[f], f))
PAL = dict(zip(fams, ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#6b3fa0"]))
labels = list(CLASSES)
ang = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist() + [0]
scores = pd.DataFrame({f: d[d.fams.map(lambda x: f in x)][labels].mean() * 100 for f in fams}).T

fig = plt.figure(figsize=(19.5, 17.2) if STYLE == "final" else (23, 11.5), facecolor="white")
outer = fig.add_gridspec(1, 2, width_ratios=[1.75, 1], wspace=0.3 if STYLE in ("bubble", "heat") else 0.13, left=0.052, right=0.99, top=0.905, bottom=0.085)
if STYLE == "final":
    outer = fig.add_gridspec(2, 1, height_ratios=[3.05, 1], hspace=0.2, left=0.075, right=0.965, top=0.925, bottom=0.075)
ga = outer[0].subgridspec(2, 3, wspace=0.62, hspace=0.5)
years = np.arange(d.y.min(), d.y.max() + 1)
series = {}
for i, f in enumerate(fams):
    r, c = divmod(i, 3)
    ax = fig.add_subplot(ga[r, c], polar=True)
    v = scores.loc[f].tolist()
    col = PAL[f]
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.plot(ang, v + v[:1], color=col, lw=3.6)
    ax.fill(ang, v + v[:1], color=col, alpha=0.35)
    ax.scatter(ang[:-1], v, s=80, color=col, zorder=3, clip_on=False)
    ax.set_xticks(ang[:-1], labels, fontsize=19)
    ax.tick_params(axis="x", pad=14)
    ax.set_ylim(0, 100)
    ax.set_yticks([50, 100], ["50", "100%"], fontsize=15, color="#55677a")
    ax.set_rlabel_position(36)
    ax.grid(color=GRID, lw=1.2)
    ax.spines["polar"].set_color(GRID)
    ax.set_title(f"{f} (n = {cnt[f]})", fontsize=23, fontweight="bold", pad=32)
    series[f] = d[d.fams.map(lambda x: f in x)].y.value_counts().reindex(years, fill_value=0)

if STYLE == "lines":
    bx = fig.add_subplot(outer[1])
    for f in fams:
        s_ = series[f].cumsum()
        bx.step(np.append(years, years.max() + 1), np.append(s_.values, s_.values[-1]), where="post", color=PAL[f], lw=4, label=f"{f} ({cnt[f]})", zorder=3)
    bx.set_xlim(years.min(), years.max() + 1)
    bx.set_ylim(0, max(cnt[f] for f in fams) * 1.06)
    bx.set_xticks(range(2006, 2026, 4))
    bx.set_yticks(range(0, 13, 2))
    bx.tick_params(labelsize=20, length=6)
    bx.grid(axis="y", color="#e3e8ee", lw=1)
    for sp in ("top", "right"):
        bx.spines[sp].set_visible(False)
    bx.set_ylabel("Cumulative papers", fontsize=23)
    bx.set_xlabel("Publication year", fontsize=23)
    bx.legend(loc="upper left", frameon=False, fontsize=21, handlelength=1.3)
    bpos = bx.get_position().x0 - 0.045
elif STYLE == "bubble":
    bx = fig.add_subplot(outer[1])
    ys = {f: len(fams) - 1 - i for i, f in enumerate(fams)}
    for f in fams:
        bx.hlines(ys[f], years.min() - 0.5, years.max() + 0.5, color="#e3e8ee", lw=1.5, zorder=1)
        for y, v in series[f].items():
            if v:
                bx.scatter(y, ys[f], s=300 + 420 * (v - 1), color=PAL[f], edgecolors="white", linewidths=1.5, zorder=3)
                if v > 1:
                    bx.text(y, ys[f], str(v), ha="center", va="center", color="white", fontsize=20, fontweight="bold", zorder=4)
    bx.set_yticks(list(ys.values()), [f"{f} ({cnt[f]})" for f in ys], fontsize=22)
    bx.set_xlim(years.min() - 0.8, years.max() + 0.8)
    bx.set_ylim(-0.7, len(fams) - 0.3)
    bx.set_xticks(range(2006, 2026, 4))
    bx.tick_params(axis="x", labelsize=20, length=6)
    bx.tick_params(axis="y", length=0, pad=6)
    for sp in ("top", "right", "left"):
        bx.spines[sp].set_visible(False)
    bx.set_xlabel("Publication year", fontsize=23)
    bpos = bx.get_position().x0 - 0.135
elif STYLE in ("areas", "final"):
    gb = outer[1].subgridspec(1, 6, wspace=0.14) if STYLE == "final" else outer[1].subgridspec(3, 2, wspace=0.12, hspace=0.42)
    first = None
    ymax = max(cnt[f] for f in fams)
    for i, f in enumerate(fams):
        r, c = (0, i) if STYLE == "final" else divmod(i, 2)
        bx = fig.add_subplot(gb[r, c])
        first = first or bx
        s_ = series[f].cumsum()
        xs, vs = np.append(years, years.max() + 1), np.append(s_.values, s_.values[-1])
        bx.fill_between(xs, vs, step="post", color=PAL[f], alpha=0.4, lw=0)
        bx.step(xs, vs, where="post", color=PAL[f], lw=3.6)
        bx.set_xlim(years.min(), years.max() + 1)
        bx.set_ylim(0, ymax * 1.08)
        bx.set_yticks([0, 4, 8, 12])
        bx.set_xticks([2010, 2018, 2025] if STYLE == "final" else [2008, 2016, 2024])
        bx.tick_params(labelsize=18, length=5, labelleft=c == 0, labelbottom=STYLE == "final" or r == 2)
        bx.grid(axis="y", color="#e3e8ee", lw=1)
        for sp in ("top", "right"):
            bx.spines[sp].set_visible(False)
        bx.set_title(f, fontsize=21, fontweight="bold", pad=6)
    pos = outer[1].get_position(fig)
    if STYLE == "final":
        first.set_ylabel("Cumulative papers", fontsize=22)
        fig.text(pos.x0 + pos.width / 2, 0.012, "Publication year", ha="center", fontsize=22)
        fig.text(0.006, pos.y1 + 0.022, "B", fontsize=34, fontweight="bold", va="bottom")
        bpos = None
    else:
        fig.text(first.get_position().x0 - 0.04, 0.5, "Cumulative papers", rotation=90, va="center", fontsize=22)
        fig.text(pos.x0 + pos.width / 2, 0.012, "Publication year", ha="center", fontsize=22)
        bpos = first.get_position().x0 - 0.045
elif STYLE == "heat":
    bx = fig.add_subplot(outer[1])
    yrs = [y for y in years if y >= 2008]
    grid = np.array([[series[f][y] for y in yrs] for f in fams])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("h", ["#ffffff", "#7fb8e6", "#0a5fa3", "#06305c"])
    bx.imshow(grid, aspect="auto", cmap=cmap, vmin=0, vmax=grid.max())
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            if grid[i, j]:
                bx.text(j, i, str(grid[i, j]), ha="center", va="center", fontsize=18, fontweight="bold", color="white" if grid[i, j] >= 2 else "#16324f")
    bx.set_xticks(np.arange(len(yrs))[::3], [str(y) for y in yrs[::3]], fontsize=19)
    bx.set_yticks(range(len(fams)), [f"{f} ({cnt[f]})" for f in fams], fontsize=21)
    bx.set_xticks(np.arange(-0.5, len(yrs)), minor=True)
    bx.set_yticks(np.arange(-0.5, len(fams)), minor=True)
    bx.grid(which="minor", color="#d4dbe3", lw=1.2)
    bx.tick_params(which="both", length=0)
    for sp in bx.spines.values():
        sp.set_visible(False)
    bx.set_xlabel("Publication year", fontsize=23)
    bpos = bx.get_position().x0 - 0.135
fig.text(0.006, 0.955, "A", fontsize=34, fontweight="bold", va="bottom")
if bpos is not None:
    fig.text(bpos, 0.955, "B", fontsize=34, fontweight="bold", va="bottom")


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


NAME = "Figure_7_Gene_Families" if STYLE == "final" else f"Figure_7_panelB_{STYLE}"
save(NAME + ".png", dpi=300)
save(NAME + ".pdf")
save(NAME + ".svg")
scores.round(1).assign(papers=[cnt[f] for f in fams]).to_csv(OUT / "Figure_7_capability_scores.csv", encoding="utf-8-sig")
pd.DataFrame(cnt.most_common(), columns=["family", "papers"]).to_csv(OUT / "Figure_7_all_families.csv", index=False, encoding="utf-8-sig")
print("family names", len(cnt), "shown", len(fams), {f: cnt[f] for f in fams})
print(scores.round(0).astype(int).to_string())
print("single-paper families", sum(v == 1 for v in cnt.values()))
print({f: int(d[d.fams.map(lambda x: f in x)].y.min()) for f in fams})
sel = d[d.fams.map(lambda x: any(f in x for f in fams))]
print("papers in shown families", len(sel), "of those 2022-25", int((sel.y >= 2022).sum()))
