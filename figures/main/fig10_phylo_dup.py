import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import collections
import math
import os
import re
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

DATA = str(DATA_FILE)
OUT = Path(OUTPUT_DIR / "Figure_10_Phylogeny_Duplication")
OUT.mkdir(parents=True, exist_ok=True)
INK, GRID = "#16324f", "#dfe5ec"
plt.rcParams.update({"font.family": "Arial", "font.size": 20, "text.color": INK, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})
TYPE_COL = {"Group": "#e0b45f", "Subfamily": "#5f9fb0", "Clade": "#6fbaa2", "Subgroup": "#a48fd0", "Class": "#d98a8a", "Cluster": "#e09a74", "Type": "#8b9aaa", "Other": "#aeb7c2"}
PERIODS = ["Before 2013", "2013–2016", "2017–2020", "2021–2025"]
PER_COL = dict(zip(PERIODS, ["#cfe3f4", "#9cc6e8", "#5f9fd3", "#2f6fae"]))

d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
d["y"] = d.Year.astype(int)
d["period"] = pd.cut(d.y, [1990, 2012, 2016, 2020, 2030], labels=PERIODS)
d["genes"] = pd.to_numeric(d["Number of Genes Identified"], errors="coerce")
N = len(d)


def clean(v):
    return " ".join(str(v).replace("\x96", "-").replace("\xa0", " ").split())


def word_to_num(t):
    for w, n in {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11, "twelve": 12, "thirteen": 13, "fourteen": 14, "fifteen": 15, "sixteen": 16, "seventeen": 17, "eighteen": 18, "nineteen": 19, "twenty": 20}.items():
        t = re.sub(rf"\b{w}\b", str(n), t, flags=re.I)
    return t


def stype(t):
    l = clean(t).lower().replace("-", "")
    for key, name in (("subfamil", "Subfamily"), ("clade", "Clade"), ("cluster", "Cluster"), ("subgroup", "Subgroup"), ("group", "Group"), ("class", "Class"), ("type", "Type")):
        if key in l:
            return name
    return "Not specified" if "not specified" in l or not l else "Other"


def parts(t):
    l = word_to_num(clean(t).lower())
    if "not specified" in l or not l:
        return np.nan
    for pat in (r"(?:into|in|classified|divided)\s*(?:into)?\s+(\d+)\s+(?:distinct\s+)?(?:main\s+|major\s+)?(?:subfamil|clade|group|class|subgroup|type|cluster|famil)", r"^\s*(\d+)\s+(?:distinct\s+)?(?:main\s+|major\s+)?(?:sub-?\s?famil|sub-?clade|clade|group|class|sub-?group|type|cluster|famil)"):
        m = re.search(pat, l)
        if m:
            return float(m.group(1))
    for n in re.findall(r"\b(\d+)\b", l):
        if 2 <= int(n) <= 60:
            return float(n)
    return np.nan


d["stype"] = d["Phylogenetic Structure"].map(stype)
d["parts"] = d["Phylogenetic Structure"].map(parts)
p = d[(d.stype != "Not specified") & d.parts.notna()].copy()
order = list(p.stype.value_counts().index)
FN = lambda v: np.log(np.maximum(v, 1.2))

fig = plt.figure(figsize=(20, 24), facecolor="white")
ax = fig.add_axes([0.14, 0.525, 0.72, 0.445], polar=True)
n = len(order)
ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
wid = 2 * np.pi / n
rng = np.random.default_rng(42)
pmax = float(p.parts.max())
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
ticks = [2, 5, 10, 20, 40]
for i, t in enumerate(order):
    s = p[p.stype == t]
    col = TYPE_COL.get(t, "#aeb7c2")
    ax.bar(ang[i], FN(pmax) * 1.05, width=wid * 0.84, color=col, alpha=0.14, zorder=0)
    sizes = np.clip(np.sqrt(s.genes.fillna(s.genes.median())) * 22, 110, 1100)
    ax.scatter(ang[i] + rng.uniform(-wid * 0.36, wid * 0.36, len(s)), FN(s.parts.values + rng.uniform(-0.18, 0.18, len(s))), s=sizes, color=col, edgecolors="white", linewidths=1.0, alpha=0.9, zorder=3)
    x, y = math.sin(ang[i]), math.cos(ang[i])
    ax.text(ang[i], FN(pmax) * 1.13, f"{t} ({len(s)})", ha="center" if abs(x) < 0.25 else ("left" if x > 0 else "right"), va="center", fontsize=25, fontweight="bold", clip_on=False)
ax.set_yticks([FN(v) for v in ticks], [str(v) for v in ticks], fontsize=20, color="#5f7186")
ax.set_rlabel_position(np.degrees(wid / 2))
ax.set_ylim(0, FN(pmax) * 1.05)
ax.set_xticks([])
ax.spines["polar"].set_visible(False)
ax.grid(color=GRID, lw=1.1)

# B
RAW = {"fragment": "Segmental"}
UNNAMED = {"duplication", "paralogous", "homoeologous copies", "type not named"}
SETS = ["Segmental", "Tandem", "WGD", "Dispersed", "Proximal", "Transposed"]


def dup(v):
    out = set()
    for x in v.split(","):
        x = x.strip()
        lo = x.lower()
        if lo in UNNAMED or lo in ("not specified", ""):
            continue
        out.add(RAW.get(lo, x))
    return frozenset(out)


col = "Gene Duplication Events (TANDEM/SEGMENTAL)"
d["dup"] = d[col].map(dup)
named = d[d.dup.map(len) > 0]
unnamed = int(d[col].map(lambda v: any(x.strip().lower() in UNNAMED for x in v.split(",")) and not dup(v)).sum())
not_spec = int((d[col] == "Not specified").sum())
combo = named.dup.value_counts()
MIN = 2
main = [c for c in combo.index if combo[c] >= MIN]
rest = [c for c in combo.index if combo[c] < MIN]
d["combo"] = d.dup.map(lambda c: c if c in main else ("rest" if len(c) else None))
labels = main + (["rest"] if rest else [])
ga = fig.add_gridspec(2, 2, left=0.17, right=0.985, top=0.475, bottom=0.035, width_ratios=[0.3, 1], height_ratios=[1.45, 1], wspace=0.03, hspace=0.04)
bar = fig.add_subplot(ga[0, 1])
mat = fig.add_subplot(ga[1, 1], sharex=bar)
side = fig.add_subplot(ga[1, 0], sharey=mat)
X = np.arange(len(labels))
bottom = np.zeros(len(labels))
tab = []
for per in PERIODS:
    v = [int(((d.combo.map(lambda c: c == l)) & (d.period == per)).sum()) for l in labels]
    bar.bar(X, v, bottom=bottom, color=PER_COL[per], width=0.74, edgecolor="white", linewidth=1, label=per, zorder=2)
    bottom += v
    tab.append(v)
for x, t in zip(X, bottom):
    bar.text(x, t + 1.2, f"{int(t)}", ha="center", va="bottom", fontsize=22, fontweight="bold")
bar.set_ylabel("Number of papers", fontsize=23)
bar.set_ylim(0, bottom.max() * 1.15)
bar.grid(axis="y", color=GRID, lw=1.1, zorder=0)
bar.tick_params(axis="x", length=0, labelbottom=False)
bar.tick_params(axis="y", labelsize=20, length=6)
for s_ in ("top", "right", "bottom"):
    bar.spines[s_].set_visible(False)
bar.legend(loc="upper right", frameon=False, fontsize=22, handlelength=1.2, title="Publication period", title_fontsize=22)
ys = np.arange(len(SETS))[::-1]
for j, y in enumerate(ys):
    if j % 2 == 0:
        mat.axhspan(y - 0.5, y + 0.5, color="#f4f6f9", zorder=0)
for x, l in zip(X, labels):
    if l == "rest":
        mat.text(x, ys.mean(), "Other\nrare\ncombinations", ha="center", va="center", fontsize=19, color="#5f7186")
        continue
    on = [ys[SETS.index(s_)] for s_ in SETS if s_ in l]
    mat.scatter([x] * len(ys), ys, s=170, color="#d5dbe2", zorder=1)
    mat.scatter([x] * len(on), on, s=260, color="#2f6fae", zorder=3)
    if len(on) > 1:
        mat.plot([x, x], [min(on), max(on)], color="#2f6fae", lw=4.5, zorder=2)
mat.set_yticks(ys, SETS, fontsize=23)
mat.set_ylim(-0.6, len(SETS) - 0.4)
mat.set_xlim(-0.6, len(labels) - 0.4)
mat.set_xticks([])
mat.tick_params(axis="y", length=0, labelleft=False)
for s_ in mat.spines.values():
    s_.set_visible(False)
tot = [int(d.dup.map(lambda c: s_ in c).sum()) for s_ in SETS]
side.barh(ys, tot, color="#8fbbe0", height=0.62, zorder=2)
for y, v in zip(ys, tot):
    side.text(v + max(tot) * 0.03, y, str(v), va="center", ha="right", fontsize=22, fontweight="bold")
side.invert_xaxis()
side.set_xlim(max(tot) * 1.2, 0)
side.yaxis.tick_right()
side.set_yticks(ys, SETS, fontsize=23)
side.set_ylim(-0.6, len(SETS) - 0.4)
side.tick_params(axis="y", length=0, pad=8)
side.set_xticks([])
for s_ in side.spines.values():
    s_.set_visible(False)
mat.yaxis.set_visible(False)
for lab, s_ in zip(side.get_yticklabels(), SETS):
    lab.set_ha("left")
for sp_ in ("left",):
    pass
mid = side.get_position()
for a_ in (bar, mat):
    q = a_.get_position()
    a_.set_position([0.345, q.y0, 0.64, q.height])
side.set_position([0.03, mid.y0, 0.215, mid.height])
fig.text(0.008, 0.975, "A", fontsize=36, fontweight="bold", va="center")
fig.text(0.008, 0.485, "B", fontsize=36, fontweight="bold", va="center")


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


save("Figure_10_Phylogeny_Duplication.png", dpi=280)
save("Figure_10_Phylogeny_Duplication.pdf")
save("Figure_10_Phylogeny_Duplication.svg")
rows = [dict(combination=" + ".join(s_ for s_ in SETS if s_ in l) if l != "rest" else "Other rare combinations", papers=int(bottom[i]), **{per: tab[k][i] for k, per in enumerate(PERIODS)}) for i, l in enumerate(labels)]
pd.DataFrame(rows).to_csv(OUT / "Figure_10B_duplication_combinations.csv", index=False, encoding="utf-8-sig")
p.groupby("stype").agg(papers=("parts", "size"), median_partitions=("parts", "median"), max_partitions=("parts", "max")).sort_values("papers", ascending=False).to_csv(OUT / "Figure_10A_structure_terms.csv", encoding="utf-8-sig")
print("structure", p.stype.value_counts().to_dict(), "total", len(p), "median", p.parts.median(), "max", pmax, "excluded", N - len(p), d[~d.index.isin(p.index)]["Phylogenetic Structure"].tolist()[:12])
print("not specified", not_spec, "unnamed", unnamed, "with named type", len(named))
print(pd.DataFrame(rows).to_string())
print("sets", dict(zip(SETS, tot)))
