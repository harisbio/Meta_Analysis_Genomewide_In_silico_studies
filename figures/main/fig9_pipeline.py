import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import collections
import os
import re
import time
from pathlib import Path as FP

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mc
import matplotlib.patheffects
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
import numpy as np
import pandas as pd

DATA = str(DATA_FILE)
OUT = FP(OUTPUT_DIR / "Figure_9_Pipeline_Tools")
OUT.mkdir(parents=True, exist_ok=True)
INK = "#16324f"
plt.rcParams.update({"font.family": "Arial", "font.size": 20, "text.color": INK, "axes.edgecolor": INK, "xtick.color": INK, "ytick.color": INK, "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})
PERIODS = ["Before 2013", "2013–2016", "2017–2020", "2021–2025"]

d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
N = len(d)
d["y"] = d.Year.astype(int)
d["period"] = pd.cut(d.y, [1990, 2012, 2016, 2020, 2030], labels=PERIODS)


def source(v):
    """First matching database; the order is the resources' overall frequency."""
    l = v.lower()
    for key, name in [("phytozome", "Phytozome"), ("ensembl", "Ensembl Plants"), ("gramene", "Ensembl Plants"), ("ncbi", "NCBI"), ("genbank", "NCBI"), ("refseq", "NCBI"),
                      ("cotton", "Cotton databases"), ("sol genomics", "Sol Genomics Network"), ("sgn", "Sol Genomics Network"), ("gdr", "GDR"), ("rosaceae", "GDR"),
                      ("brad", "BRAD"), ("brassica", "BRAD"), ("cruciferae", "BRAD"), ("rice genome annotation", "RGAP"), ("rgap", "RGAP"), ("msu", "RGAP"), ("tigr", "RGAP")]:
        if key in l:
            return name
    return "Not stated" if l.strip() in ("", "not available") else "Other databases"


def first_species(v):
    t = re.sub(r"[\xa0\u200b]", " ", v).split(",")[0]
    m = re.match(r"\s*([A-Z][a-z\-]+)\s+(?:[x×]\s+)?([a-z][a-z\-]+)", t)
    return f"{m.group(1)} {m.group(2)}" if m else t.strip()


q, r_ = d["Experimental Validation (qPCR)"] == "Yes", d["Expression Analysis (RNA-seq)"] == "Yes"
DESIGN = {"In silico + RNA-seq + qPCR": "#6c9fce", "In silico + qPCR": "#e5bd73", "In silico + RNA-seq": "#7cc2aa", "Only in silico": "#a9b3bd"}
d["design"] = np.select([q & r_, q, r_], list(DESIGN)[:3], list(DESIGN)[3])
d["source"] = d["Genome Source"].map(source)
d["org"] = d["Primary Organism"].map(first_species)
top = d.org.value_counts()
keep = list(top.index[:8])
d["org2"] = d.org.where(d.org.isin(keep), "Other organisms")
SRC_ORDER = [s for s in d.source.value_counts().index if s not in ("Other databases", "Not stated")] + [x for x in ("Other databases", "Not stated") if x in set(d.source)]
SRC_RANK = {k: i for i, k in enumerate([s_ for s_ in d.source.value_counts().index if s_ != "Other databases"] + ["Other databases"])}
ORG_ORDER = sorted(keep, key=lambda o: d[d.org2 == o].source.map(SRC_RANK).mean()) + ["Other organisms"]
DES_ORDER = list(DESIGN)
SRC_COL = dict(zip(SRC_ORDER, ["#7ba7d1", "#7dbda8", "#a896cc", "#dfa582", "#ddb0c8", "#e2c989", "#a3cfe6", "#c8b08c", "#d5dae0", "#d5dae0"]))

fig = plt.figure(figsize=(20, 24), facecolor="white")
gs = fig.add_gridspec(2, 1, height_ratios=[1.16, 1], hspace=0.09, left=0.015, right=0.985, top=0.975, bottom=0.02)
ax = fig.add_subplot(gs[0])
ax.axis("off")
XS = [0.215, 0.50, 0.775]
W, GAP = 0.016, 9


def stack(order, counts, gap=GAP):
    out = {}
    total = sum(counts[k] for k in order) + gap * (len(order) - 1)
    y = total
    for k in order:
        out[k] = [y - counts[k], y]
        y -= counts[k] + gap
    return out, total


def ribbon(x0, y0a, y0b, x1, y1a, y1b, col, alpha):
    mid = (x0 + x1) / 2
    verts = [(x0, y0b), (mid, y0b), (mid, y1b), (x1, y1b), (x1, y1a), (mid, y1a), (mid, y0a), (x0, y0a), (x0, y0b)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=col, edgecolor="white", linewidth=0.6, alpha=alpha, zorder=1))


cs, co, cd = d.source.value_counts(), d.org2.value_counts(), d.design.value_counts()
ns, tot_s = stack(SRC_ORDER, cs)
no, tot_o = stack(ORG_ORDER, co, gap=23)
nd, tot_d = stack(DES_ORDER, cd, gap=20)
top_h = max(tot_s, tot_o, tot_d)
for nodes_, tot_ in ((ns, tot_s),):
    for k in nodes_:
        nodes_[k] = [nodes_[k][0] + (top_h - tot_) / 2, nodes_[k][1] + (top_h - tot_) / 2]
offset_o, offset_d = (top_h - tot_o) / 2, (top_h - tot_d) / 2
no = {k: [v[0] + offset_o, v[1] + offset_o] for k, v in no.items()}
nd = {k: [v[0] + offset_d, v[1] + offset_d] for k, v in nd.items()}
ptr_s, ptr_o_in, ptr_o_out, ptr_d = {k: v[1] for k, v in ns.items()}, {k: v[1] for k, v in no.items()}, {k: v[1] for k, v in no.items()}, {k: v[1] for k, v in nd.items()}
for s_ in SRC_ORDER:
    for o in ORG_ORDER:
        n = int(((d.source == s_) & (d.org2 == o)).sum())
        if n:
            ribbon(XS[0] + W, ptr_s[s_] - n, ptr_s[s_], XS[1], ptr_o_in[o] - n, ptr_o_in[o], SRC_COL[s_], 0.6)
            ptr_s[s_] -= n
            ptr_o_in[o] -= n
for o in ORG_ORDER:
    for k in DES_ORDER:
        n = int(((d.org2 == o) & (d.design == k)).sum())
        if n:
            ribbon(XS[1] + W, ptr_o_out[o] - n, ptr_o_out[o], XS[2], ptr_d[k] - n, ptr_d[k], DESIGN[k], 0.5)
            ptr_o_out[o] -= n
            ptr_d[k] -= n
for k, (a, b) in ns.items():
    ax.add_patch(Rectangle((XS[0], a), W, b - a, facecolor=SRC_COL[k], edgecolor="none", zorder=2))
    ax.text(XS[0] - 0.008, (a + b) / 2, f"{k} ({cs[k]})", ha="right", va="center", fontsize=22, fontweight="bold")
for k, (a, b) in no.items():
    ax.add_patch(Rectangle((XS[1], a), W, b - a, facecolor="#8a97a5", edgecolor="none", zorder=2))
    ax.text(XS[1] + W / 2, b + 3, f"{k} ({co[k]})", ha="center", va="bottom", fontsize=21, fontweight="bold", fontstyle="normal" if k.startswith("Other") else "italic", zorder=4,
            path_effects=[matplotlib.patheffects.withStroke(linewidth=7, foreground="white")])
for k, (a, b) in nd.items():
    ax.add_patch(Rectangle((XS[2], a), W, b - a, facecolor=DESIGN[k], edgecolor="none", zorder=2))
    ax.text(XS[2] + W + 0.008, (a + b) / 2, f"{k.replace(' + ', ' +\n', 1) if len(k) > 16 else k}\n({cd[k]})", ha="left", va="center", fontsize=22, fontweight="bold", linespacing=1.08)
ax.set_xlim(0, 1)
ax.set_ylim(-4, top_h + 36)
for x, lab in zip([XS[0] - 0.06, XS[1] + W / 2, XS[2] + 0.07], ["Genome source", "Main study organism", "Study design"]):
    ax.text(x, top_h + 31, lab, ha="center", va="center", fontsize=26, fontweight="bold")

# B
col_or = lambda c: d[c].map(lambda v: v.strip().lower() not in ("", "no", "not used", "none"))


def has(col, word):
    return d[col].str.contains(rf"\b{word}\b", case=False, regex=True)


TOOLS = [("Identification", "BLAST", col_or("BLAST (E-value cut-off used)")), ("Identification", "HMMER", col_or("HMMER (with Pfam domains ID)")),
         ("Domain", "Pfam", has("Other (Tools for Domain Analysis)", "Pfam")), ("Domain", "SMART", d.SMART == "Yes"), ("Domain", "CDD", d.CDD == "Yes"),
         ("Motif", "MEME", d.MEME == "Yes"), ("Phylogeny", "MEGA", d.MEGA == "Yes"), ("Phylogeny", "IQ-TREE", d["IQ-TREE"] == "Yes"),
         ("Synteny", "MCScanX", d.MCScanX == "Yes"), ("Synteny", "PGDD", has("Other (Tools for Synteny)", "Plant Genome Duplication Database")),
         ("Visualization", "TBtools", d.TBTools == "Yes"), ("Visualization", "GSDS", has("Other (Tools for Visualization)", "GSDS")),
         ("Other", "ExPASy", has("Other Tools Used In Methodology", "ExPASy")), ("Other", "PlantCARE", has("Other Tools Used In Methodology", "PlantCARE")),
         ("Other", "ClustalW", has("Other Tools Used In Methodology", "ClustalW")), ("Other", "WoLF PSORT", has("Other Tools Used In Methodology", "WoLF PSORT"))]
rows = []
for group, name, mask in TOOLS:
    row = dict(step=group, tool=name, papers=int(mask.sum()), overall=mask.mean() * 100)
    for p in PERIODS:
        row[p] = mask[d.period == p].mean() * 100
    rows.append(row)
t = pd.DataFrame(rows)
SOFT = "#8fbbe0"
STEP_COL = {k: SOFT for k in ("Identification", "Domain", "Motif", "Phylogeny", "Synteny", "Visualization", "Other")}
bx = fig.add_subplot(gs[1])
nt = len(t)
ys_ = np.arange(nt)[::-1]
X_NAME, X_P, X_BAR, BAR_W = 0.0, [1.15, 2.15, 3.15, 4.15], 5.05, 2.4
for i, row in t.iterrows():
    y = ys_[i]
    col = STEP_COL[row.step]
    if i % 2 == 0:
        bx.add_patch(Rectangle((-1.75, y - 0.5), 9.6, 1, facecolor="#f6f8fa", edgecolor="none", zorder=0))
    bx.text(X_NAME + 0.42, y, row.tool, ha="right", va="center", fontsize=23, fontweight="bold")
    for x, p in zip(X_P, PERIODS):
        v = row[p]
        bx.scatter(x, y, s=max(v, 0.3) * 19, color=col, edgecolors="none", zorder=2)
        big = v >= 28
        bx.text(x, y, f"{v:.0f}", ha="center" if big else "left", va="center", fontsize=20, fontweight="bold", color=INK, zorder=3) if big else bx.text(x + 0.05 + 0.0052 * max(v, 0.3) ** 0.5 * 7, y, f"{v:.0f}", ha="left", va="center", fontsize=20, color=INK, zorder=3)
    bx.add_patch(Rectangle((X_BAR, y - 0.3), BAR_W * row.overall / 100, 0.6, facecolor=col, edgecolor="none", zorder=2))
    bx.text(X_BAR + BAR_W * row.overall / 100 + 0.05, y, f"{row.overall:.0f}% ({row.papers})", ha="left", va="center", fontsize=21, fontweight="bold")
prev, start = None, 0
for i, g in enumerate(list(t.step) + [None]):
    if g != prev:
        if prev is not None:
            yt, yb = ys_[start] + 0.38, ys_[i - 1] - 0.38
            bx.add_patch(Rectangle((-0.12 - 0.09 + -0.0, yb), 0.0, 0.0))
            bx.plot([-1.72, -1.72], [yb, yt], color="#c3ccd6", lw=5, solid_capstyle="round", zorder=3, clip_on=False)
            bx.text(-1.62, (yt + yb) / 2, prev, ha="left", va="center", fontsize=21, fontweight="bold", color=INK)
        prev, start = g, i
for x, p in zip(X_P, PERIODS):
    bx.text(x, nt - 0.02, f"{p}\n(n = {(d.period == p).sum()})", ha="center", va="bottom", fontsize=21, fontweight="bold")
bx.text(X_BAR + BAR_W / 2, nt - 0.02, "All papers", ha="center", va="bottom", fontsize=21, fontweight="bold")
bx.set_xlim(-1.8, X_BAR + BAR_W + 0.55)
bx.set_ylim(-0.6, nt + 0.75)
bx.axis("off")
fig.text(0.006, 0.982, "A", fontsize=34, fontweight="bold", va="center")
fig.text(0.006, bx.get_position().y1 + 0.004, "B", fontsize=34, fontweight="bold", va="center")


def save(name, **kw):
    tmp = OUT / (name + ".tmp" + FP(name).suffix)
    fig.savefig(tmp, **kw)
    for _ in range(8):
        try:
            os.replace(tmp, OUT / name)
            return
        except OSError:
            time.sleep(5)
    raise OSError(f"{name} is locked; new version left as {tmp.name}")


save("Figure_9_Pipeline_Tools.png", dpi=300)
save("Figure_9_Pipeline_Tools.pdf")
save("Figure_9_Pipeline_Tools.svg")
t.round(1).to_csv(OUT / "Figure_9_tool_use.csv", index=False, encoding="utf-8-sig")
d.groupby(["source", "org2", "design"], observed=True).size().rename("papers").reset_index().to_csv(OUT / "Figure_9_flows.csv", index=False, encoding="utf-8-sig")
print(cs.to_dict())
print(co.to_dict())
print(cd.to_dict())
print(t.round(0).to_string())
oth = d[d.source == "Other databases"]["Genome Source"].value_counts()
print(oth.head(25).to_dict())
b = d["BLAST (E-value cut-off used)"]
print(b.value_counts().to_dict())
