import math, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import global_analysis as ga  # noqa: E402
from global_analysis import C, OTHER, INK, INK2, GRID, SURF, PERIODS, LABEL, TAB, save, plt  # noqa: E402
from family_names import family_norm, key as family_key  # noqa: E402

PER = [p[2] for p in PERIODS]


def load():
    d = pd.read_csv(ga.GLOBAL_DATA_FILE, dtype=str, keep_default_na=False)
    a = d[d.in_analysis == "True"].copy()
    a["yr"] = pd.to_numeric(a.year, errors="coerce")
    a = a[a.yr.between(2000, 2025)].copy(); a["yr"] = a.yr.astype(int); a["period"] = a.yr.map(ga.period)
    for f in ga.LISTS:
        a[f] = a[f].map(lambda s: [x for x in s.split("; ") if x])
    a["cites"] = pd.to_numeric(a.cited_by, errors="coerce"); a["fwci_n"] = pd.to_numeric(a.fwci, errors="coerce")
    a["n_sp"] = pd.to_numeric(a.n_species_norm, errors="coerce"); a["genes"] = pd.to_numeric(a.genes_total, errors="coerce")
    fm = a[["key", "yr", "period", "gene_family", "family_class", "corr_country"]].assign(f=a.gene_family.str.split(r",\s*")).explode("f")
    fm = fm[fm.f.fillna("").str.strip().ne("") & fm.f.ne("Not stated")]
    disp = fm.assign(k=fm.f.map(family_key)).groupby("k").f.agg(lambda s: s.value_counts().index[0]).to_dict()
    fm["group"] = [family_norm(f, c, disp)[1] for f, c in zip(fm.f, fm.family_class)]
    fm = fm.drop_duplicates(["key", "group"])
    sp = a[["key", "yr", "period", "species_norm"]].assign(s=a.species_norm.str.split(", ")).explode("s")
    sp = sp[sp.s.notna() & ~sp.s.isin(["Not stated", ""])]
    return a, fm, sp


def heat(ax, m, rows, cols, fmt="{:.0f}", vmax=None, cmap="Blues", cbar_label="", fig=None):
    vmax = vmax or max(m.max(), 1)
    im = ax.imshow(m, cmap=cmap, aspect="auto", vmin=0, vmax=vmax)
    ax.set_xticks(range(len(cols)), cols); ax.set_yticks(range(len(rows)), rows)
    ax.grid(False); ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if m[i, j] > 0:
                ax.text(j, i, fmt.format(m[i, j]), ha="center", va="center", fontsize=6.8,
                        color=SURF if m[i, j] > 0.55 * vmax else INK)
    if fig is not None:
        cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label(cbar_label, fontsize=8); cb.outline.set_visible(False)
    return im


def main():
    a, fm, sp = load()
    N = len(a)

    # 13 Bradford curve
    j = pd.read_csv(TAB / "04_journals_bradford.csv")
    j["rank"] = np.arange(1, len(j) + 1); j["cum_pct"] = j.cum_papers / j.papers.sum() * 100
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.plot(j["rank"], j.cum_pct, color=C[0], lw=2); ax.set_xscale("log")
    ax.xaxis.set_major_formatter(ga.matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    for z, col in zip((1, 2), (C[1], C[2])):
        r = j[j.bradford_zone == z]["rank"].max()
        ax.axvline(r, color=INK2, lw=0.8, ls="--")
        ax.text(r * 1.08, 8 + 10 * z, f"zone {z} ends: {r} journals", fontsize=7.5, color=INK2)
    for y in (33.3, 66.7): ax.axhline(y, color=GRID, lw=0.8)
    ax.set_xlabel("Journals ranked by number of papers (log scale)"); ax.set_ylabel("Cumulative share of papers (%)")
    ax.set_ylim(0, 102); ax.set_title(f"Bradford curve: {len(j)} journals, three zones of equal papers", loc="left", fontsize=10)
    save(fig, "13_bradford_curve")

    # 14 collaboration network (top 20 countries by international co-papers)
    pairs = pd.read_csv(TAB / "03c_country_pairs.csv")
    strength = pd.concat([pairs.groupby("country_a").co_papers.sum(), pairs.groupby("country_b").co_papers.sum()]).groupby(level=0).sum()
    top = list(strength.sort_values(ascending=False).head(20).index)
    pe = pairs[pairs.country_a.isin(top) & pairs.country_b.isin(top) & (pairs.co_papers >= 5)]
    papers = a.corr_country.value_counts()
    order = sorted(top, key=lambda c: -strength[c])
    ang = {c: 2 * math.pi * i / len(order) + math.pi / 2 for i, c in enumerate(order)}
    pos = {c: (math.cos(t), math.sin(t)) for c, t in ang.items()}
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    mx = pe.co_papers.max()
    for r in pe.sort_values("co_papers").itertuples():
        (x1, y1), (x2, y2) = pos[r.country_a], pos[r.country_b]
        ax.plot([x1, x2], [y1, y2], color=C[0], alpha=0.15 + 0.7 * r.co_papers / mx, lw=0.5 + 6 * r.co_papers / mx,
                solid_capstyle="round", zorder=1)
    for c, (x, y) in pos.items():
        n = papers.get(c, 0)
        ax.scatter(x, y, s=30 + 10 * math.sqrt(n), color=C[1] if c == "China" else C[0], edgecolor=SURF, linewidth=1.5, zorder=3)
        t = ang[c]; ha = "left" if math.cos(t) > 0.1 else "right" if math.cos(t) < -0.1 else "center"
        rr = 1.13 + 0.0016 * math.sqrt(n)                 # push the label past large dots
        ax.text(rr * math.cos(t), rr * math.sin(t), f"{c}\n{int(strength[c]):,}", ha=ha, va="center", fontsize=7.5, color=INK)
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.45, 1.6); ax.axis("off"); ax.set_aspect("equal")
    ax.set_title("International co-authorship among the 20 most connected countries\n"
                 "(line width = co-authored papers, ≥ 5; number = all international co-papers; dot size = papers)", loc="left", fontsize=9.5)
    pe.to_csv(TAB / "14_network_edges.csv", index=False)
    save(fig, "14_collaboration_network")

    # 15 citations by publication year (two measures -> two panels)
    cy = pd.read_csv(TAB / "12b_citations_by_year.csv").set_index("yr").loc[2005:2025]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))
    axs[0].bar(cy.index, cy.median_citations, color=C[0], width=0.78, edgecolor=SURF)
    axs[0].set_ylabel("Median citations per paper"); axs[0].grid(axis="x", visible=False)
    axs[0].set_title("Citations (raw: older papers have had longer)", loc="left", fontsize=10)
    axs[1].plot(cy.index, cy.median_fwci, color=C[0], lw=2, marker="o", ms=3)
    axs[1].axhline(1, color=INK2, lw=0.8, ls="--"); axs[1].text(2005, 1.05, "world average", fontsize=7, color=INK2)
    axs[1].set_ylabel("Median FWCI"); axs[1].set_ylim(0, max(cy.median_fwci) * 1.15)
    axs[1].set_title("Field-weighted citation impact (age-adjusted)", loc="left", fontsize=10)
    for x in axs: x.xaxis.set_major_locator(ga.matplotlib.ticker.MaxNLocator(integer=True))
    save(fig, "15_citations_by_year")

    # 16 top 20 families: share of papers per period
    top20 = fm.group.value_counts().head(20).index
    t = fm[fm.group.isin(top20)].groupby(["group", "period"]).size().unstack(fill_value=0).reindex(columns=PER, fill_value=0)
    tot = a.groupby("period").size().reindex(PER)
    tp = (t / tot * 100).loc[top20]
    tp.round(2).to_csv(TAB / "16_top_families_pct_by_period.csv")
    fig, ax = plt.subplots(figsize=(6, 7))
    heat(ax, tp.values, list(tp.index), PER, fmt="{:.1f}", cbar_label="% of papers in the period", fig=fig)
    ax.set_title("Top 20 gene families: share of papers by period (%)", loc="left", fontsize=10)
    save(fig, "16_top_families_by_period")

    # 17 species x family matrix (top 15 x top 15)
    m = sp[["key", "s"]].merge(fm[["key", "group"]], on="key")
    ts = sp.s.value_counts().head(15).index; tf = fm.group.value_counts().head(15).index
    mm = m[m.s.isin(ts) & m.group.isin(tf)].groupby(["s", "group"]).key.nunique().unstack(fill_value=0).reindex(index=ts, columns=tf, fill_value=0)
    mm.to_csv(TAB / "17_species_by_family.csv")
    fig, ax = plt.subplots(figsize=(9, 6))
    heat(ax, mm.values.astype(float), [f"$\\it{{{s.replace(' ', '~')}}}$" for s in mm.index], list(mm.columns), cbar_label="Papers", fig=fig)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_title("Papers per species × gene family (top 15 each; blank = not yet studied)", loc="left", fontsize=10)
    save(fig, "17_species_by_family")

    # 18 species diversification over time
    top10 = sp.s.value_counts().head(10).index
    ky = a.set_index("key").yr
    by_paper = sp.groupby("key").s.apply(set)
    share_top = by_paper.map(lambda s: bool(s & set(top10))).groupby(ky).mean().loc[2005:] * 100
    distinct = sp.groupby("yr").s.nunique().loc[2005:]
    pd.DataFrame({"pct_papers_on_top10_species": share_top.round(1), "distinct_species_in_year": distinct}).to_csv(TAB / "18_species_diversification.csv")
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))
    axs[0].plot(share_top.index, share_top.values, color=C[0], lw=2)
    axs[0].set_ylim(0, 100); axs[0].set_ylabel("% of papers")
    axs[0].set_title("Papers on any of the 10 most studied species", loc="left", fontsize=10)
    axs[1].bar(distinct.index, distinct.values, color=C[0], width=0.78, edgecolor=SURF); axs[1].grid(axis="x", visible=False)
    axs[1].set_ylabel("Distinct species"); axs[1].set_title("Distinct species studied per year", loc="left", fontsize=10)
    for x in axs: x.xaxis.set_major_locator(ga.matplotlib.ticker.MaxNLocator(integer=True))
    save(fig, "18_species_diversification")

    # 19 family size per family (single-species papers, top 20 families)
    first = fm.drop_duplicates("key").set_index("key").group
    g1 = a[(a.n_sp == 1) & (a.genes > 0)].assign(group=lambda x: x.key.map(first)).dropna(subset=["group"])
    tf20 = [f for f in fm.group.value_counts().index if (g1.group == f).sum() >= 15][:20]
    data = [g1[g1.group == f].genes.values for f in tf20]
    order = np.argsort([np.median(x) for x in data])
    g1[g1.group.isin(tf20)].groupby("group").genes.describe().round(1).to_csv(TAB / "19_family_size_by_family.csv")
    fig, ax = plt.subplots(figsize=(7, 6))
    bp = ax.boxplot([data[i] for i in order], vert=False, widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color=INK, lw=1.2), boxprops=dict(facecolor=C[0], edgecolor=C[0], alpha=0.85),
                    whiskerprops=dict(color=INK2), capprops=dict(color=INK2))
    ax.set_yticks(range(1, len(order) + 1), [f"{tf20[i]} (n={len(data[i])})" for i in order]); ax.set_xscale("log")
    ax.xaxis.set_major_formatter(ga.matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.grid(axis="y", visible=False); ax.set_xlabel("Members identified in one species (log scale; box = middle 50%)")
    ax.set_title("Family size by gene family (single-species papers)", loc="left", fontsize=10)
    save(fig, "19_family_size_by_family")

    # 20 multi-species studies over time
    ms = a.dropna(subset=["n_sp"]).groupby("yr").agg(papers=("key", "size"), multi_pct=("n_sp", lambda x: (x > 1).mean() * 100),
                                                     mean_species=("n_sp", "mean")).loc[2005:].round(2)
    ms.to_csv(TAB / "20_multispecies_by_year.csv")
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.plot(ms.index, ms.multi_pct, color=C[0], lw=2); ax.set_ylim(0, max(ms.multi_pct) * 1.3)
    ax.set_ylabel("% of papers with ≥ 2 species")
    ax.set_title("Comparative (multi-species) studies per year", loc="left", fontsize=10)
    ax.xaxis.set_major_locator(ga.matplotlib.ticker.MaxNLocator(integer=True))
    save(fig, "20_multispecies_by_year")

    # 21 country x family specialisation (row % of each country's papers)
    tc = a[a.corr_country != ""].corr_country.value_counts().head(10).index; tfam = fm.group.value_counts().head(12).index
    cf = fm[fm.corr_country.isin(tc) & fm.group.isin(tfam)].groupby(["corr_country", "group"]).size().unstack(fill_value=0)
    cf = cf.reindex(index=tc, columns=tfam, fill_value=0)
    cfp = cf.div(a.corr_country.value_counts()[tc], axis=0) * 100
    cfp.round(2).to_csv(TAB / "21_country_by_family_pct.csv")
    fig, ax = plt.subplots(figsize=(9, 4.6))
    heat(ax, cfp.values, [f"{c} ({a.corr_country.value_counts()[c]:,})" for c in tc], list(tfam), fmt="{:.1f}",
         cbar_label="% of the country's papers", fig=fig)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_title("Which families each country studies (% of its papers; top 10 countries × top 12 families)", loc="left", fontsize=10)
    save(fig, "21_country_by_family")

    # 22 stress conditions per year (abiotic 4 + biotic combined)
    yrs = a[a.yr >= 2008]
    cond = {"Drought": ["DROUGHT"], "Salt": ["SALT"], "Cold": ["COLD"], "Heat": ["HEAT"],
            "Biotic (fungal, bacterial, viral, insect)": ["FUNGAL", "BACTERIAL", "VIRAL", "INSECT/NEMATODE"]}
    cy2 = pd.DataFrame({k: yrs.groupby("yr").conditions.apply(lambda s: s.map(lambda l: any(c in l for c in v)).mean() * 100)
                        for k, v in cond.items()}).round(1)
    cy2.to_csv(TAB / "22_conditions_by_year.csv")
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    for i, k in enumerate(cond):
        ax.plot(cy2.index, cy2[k], color=C[i], lw=2, label=k)
    ax.set_ylabel("% of papers"); ax.set_ylim(0, cy2.values.max() * 1.25)
    ax.legend(ncol=3, fontsize=7.5, loc="upper left")
    ax.set_title("Stress conditions studied, per year", loc="left", fontsize=10)
    ax.xaxis.set_major_locator(ga.matplotlib.ticker.MaxNLocator(integer=True))
    save(fig, "22_conditions_by_year")

    # 23 functional validation per year
    val = {"Any validation": None, "Overexpression": "OVEREXPRESSION", "Subcellular localization (exp.)": "SUBCELLULAR_LOCALIZATION",
           "Silencing (VIGS/RNAi)": "SILENCING", "Knockout/mutant": "KNOCKOUT"}
    vy = pd.DataFrame({k: yrs.groupby("yr").validation.apply(lambda s: (s.map(bool) if v is None else s.map(lambda l: v in l)).mean() * 100)
                       for k, v in val.items()}).round(1)
    vy.to_csv(TAB / "23_validation_by_year.csv")
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    for i, k in enumerate(val):
        ax.plot(vy.index, vy[k], color=C[i], lw=2.4 if i == 0 else 1.6, label=k)
    ax.set_ylabel("% of papers"); ax.set_ylim(0, vy.values.max() * 1.25)
    ax.legend(ncol=3, fontsize=7.5, loc="upper left")
    ax.set_title("Experimental validation named in the abstract, per year", loc="left", fontsize=10)
    ax.xaxis.set_major_locator(ga.matplotlib.ticker.MaxNLocator(integer=True))
    save(fig, "23_validation_by_year")
    print("figures 13-23 written")


if __name__ == "__main__":
    main()
