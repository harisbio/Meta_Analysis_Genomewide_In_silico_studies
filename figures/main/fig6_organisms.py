import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import collections
import math
import os
import re
import time
from pathlib import Path as FP

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import matplotlib.patheffects as pe
import pandas as pd

DATA = str(DATA_FILE)
OUT = FP(OUTPUT_DIR / "Figure_6_Organism_Network")
OUT.mkdir(parents=True, exist_ok=True)
INK, BLUE, GREEN, EDGE = "#16324f", "#0a5fa3", "#00875f", "#5b7f9e"
MIN_PRIMARY, N_REF = 8, 8
plt.rcParams.update({"font.family": "Arial", "text.color": INK, "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})
COMMON = {"apple": "Malus domestica", "chickpea": "Cicer arietinum", "hordeum vulgare": "Hordeum vulgare", "hordium vulgare": "Hordeum vulgare", "malus x domestica": "Malus domestica", "malus × domestica": "Malus domestica",
          "fortunella hindisi": "Fortunella hindsii", "gossypium arboretum": "Gossypium arboreum", "physcomitrium patens": "Physcomitrella patens"}
SKIP = {"", "none", "not available", "not specified", "nan", "others"}


def species(text):
    """Return one clean name per species in a cell."""
    text = re.sub(r"[\xa0\u200b]", " ", str(text))
    text = re.sub(r"\([^)]*\)", " ", text).replace(" and ", ",")
    out = []
    for t in re.split(r"[,;/]", text):
        t = " ".join(t.split()).strip(" .")
        if t.lower() in SKIP:
            continue
        low = t.lower()
        if low in COMMON:
            out.append(COMMON[low])
            continue
        m = re.match(r"^([A-Z][a-zA-Z\-]+)\s+(?:[x×]\s+)?([a-z][a-z\-]+)", t)
        name = f"{m.group(1)} {m.group(2)}" if m and m.group(2) not in ("spp", "sp", "hybrid", "cv", "var") else t
        name = COMMON.get(name.lower(), name)
        out.append(name)
    return list(dict.fromkeys(out))


d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
N = len(d)
P, S = collections.Counter(), collections.Counter()
E = collections.Counter()
pairs_any = 0
for pr, se in zip(d["Primary Organism"], d["Secondary Organism"]):
    ps, ss = species(pr), [s for s in species(se)]
    ss = [s for s in ss if s not in ps]
    pairs_any += bool(ss)
    P.update(ps)
    S.update(ss)
    for a in ps:
        for b in ss:
            E[(a, b)] += 1
outer = sorted([a for a, v in P.items() if v >= MIN_PRIMARY], key=lambda a: (-P[a], a))
inner = sorted(S, key=lambda a: (-S[a], a))[:N_REF]


def abbreviate(name):
    parts = name.split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else name


fig, ax = plt.subplots(figsize=(19, 17.5), facecolor="white")
fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
ax.set_aspect("equal")
ax.axis("off")
R_OUT, R_IN = 1.0, 0.5


def place(names, radius, start=90.0):
    return {n: (radius * math.cos(math.radians(start - 360 * i / len(names))), radius * math.sin(math.radians(start - 360 * i / len(names)))) for i, n in enumerate(names)}


po = place(outer, R_OUT)
pi_ = place(inner, R_IN, start=112.5)
wmax = max(E[(a, b)] for a in outer for b in inner)
for a in outer:
    for b in inner:
        w = E[(a, b)]
        if not w:
            continue
        x0, y0 = po[a]
        x1, y1 = pi_[b]
        path = Path([(x0, y0), (x0 * 0.62, y0 * 0.62), (x1 * 1.28, y1 * 1.28), (x1, y1)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        ax.add_patch(PathPatch(path, fill=False, edgecolor=EDGE, lw=1.0 + 7 * (w / wmax), alpha=0.3 + 0.55 * (w / wmax), zorder=1, capstyle="round"))


def size(v, vmax):
    return 170 + 2300 * (v / vmax)


pmax, smax = max(P[a] for a in outer), max(S[a] for a in inner)
for a in outer:
    ax.scatter(*po[a], s=size(P[a], pmax), color=BLUE, edgecolors="white", linewidths=2, zorder=3)
    x, y = po[a]
    r = R_OUT + 0.035 + 0.0020 * math.sqrt(size(P[a], pmax))
    lx, ly = x / R_OUT * r, y / R_OUT * r
    ha = "center" if abs(x) < 0.12 else ("left" if x > 0 else "right")
    va = "bottom" if y > 0.93 else ("top" if y < -0.93 else "center")
    ax.text(lx, ly, f"{a} ({P[a]})", ha=ha, va=va, fontsize=24, fontstyle="italic", fontweight="bold", zorder=5)
for b in inner:
    ax.scatter(*pi_[b], s=size(S[b], smax), color=GREEN, edgecolors="white", linewidths=2, zorder=3)
    x, y = pi_[b]
    r = R_IN + 0.03 + 0.0019 * math.sqrt(size(S[b], smax))
    lx, ly = x / R_IN * r, y / R_IN * r
    ha = "left" if x > 0 else "right"
    va = "bottom" if y > 0.2 else ("top" if y < -0.2 else "center")
    ax.text(lx, ly, f"{abbreviate(b)} ({S[b]})", ha=ha, va=va, fontsize=22, fontstyle="italic", fontweight="bold", color="#00563c", zorder=6,
            path_effects=[pe.withStroke(linewidth=7, foreground="white")])
ax.set_xlim(-1.72, 1.72)
ax.set_ylim(-1.36, 1.2)
ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=19, color=BLUE, label="Study organisms"), Line2D([], [], marker="o", ls="", ms=19, color=GREEN, label="Reference organisms"),
                   Line2D([], [], color=EDGE, lw=4, label="Compared in the same paper")], loc="lower center", ncol=3, frameon=False, fontsize=23, bbox_to_anchor=(0.5, 0.0), columnspacing=1.6)


def save(name, **kw):
    tmp = OUT / (name + ".tmp" + FP(name).suffix)
    fig.savefig(tmp, bbox_inches="tight", pad_inches=0.2, **kw)
    for _ in range(8):
        try:
            os.replace(tmp, OUT / name)
            return
        except OSError:
            time.sleep(5)
    raise OSError(f"{name} is locked; new version left as {tmp.name}")


save("Figure_6_Organism_Network.png", dpi=400)
save("Figure_6_Organism_Network.pdf")
save("Figure_6_Organism_Network.svg")
pd.DataFrame(P.most_common(), columns=["study_organism", "papers"]).assign(percent=lambda x: (x.papers / N * 100).round(1)).to_csv(OUT / "Figure_6_study_organisms.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(S.most_common(), columns=["reference_organism", "papers"]).assign(percent=lambda x: (x.papers / N * 100).round(1)).to_csv(OUT / "Figure_6_reference_organisms.csv", index=False, encoding="utf-8-sig")
pd.DataFrame([dict(study=a, reference=b, papers=v) for (a, b), v in E.most_common()]).to_csv(OUT / "Figure_6_links.csv", index=False, encoding="utf-8-sig")
print("distinct study organisms", len(P), "shown", len(outer), "references", len(S), "papers with a reference", pairs_any)
print("top study", [(a, P[a], round(P[a] / N * 100, 1)) for a in outer[:12]])
print("top ref", [(a, S[a], round(S[a] / N * 100, 1)) for a in inner])
print("odd study names", [a for a in P if not re.fullmatch(r"[A-Z][a-z\-]+ [a-z\-]+", a)][:25], "odd ref", [a for a in S if not re.fullmatch(r"[A-Z][a-z\-]+ [a-z\-]+", a)][:25])
