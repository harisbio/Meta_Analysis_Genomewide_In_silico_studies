import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import json, re, sys, time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from family_names import family_norm, key as family_key  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
HERE = OUTPUT_DIR / "Supplementary"
TAB = HERE / "tables"; FIG = HERE / "figures"
TAB.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)

# palette: reference categorical order (validated), one sequential hue, neutral ink
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER = "#a3a29c"; INK = "#0b0b0b"; INK2 = "#52514e"; GRID = "#e4e3df"; SURF = "#ffffff"
SEQ = "Blues"
PERIODS = [(2000, 2009, "2000–2009"), (2010, 2014, "2010–2014"), (2015, 2019, "2015–2019"), (2020, 2025, "2020–2025")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.dpi": 300, "legend.frameon": False})
plt.rcParams.update({"font.size": 10.5, "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10})
from matplotlib.axes import Axes as _Axes
from matplotlib.ticker import MaxNLocator as _MNL
_orig_title = _Axes.set_title
def _no_title(self, *a, **k):
    self._had_title = True
    return self.title
_Axes.set_title = _no_title
def _finish(fig):
    t = [a for a in fig.axes if getattr(a, "_had_title", False)]
    if len(t) >= 2:
        for i, a in enumerate(t):
            _orig_title(a, "ABCDEFGH"[i], loc="left", fontsize=14, fontweight="bold", pad=8)
    for a in fig.axes:
        for lab in list(a.get_xticklabels()) + list(a.get_yticklabels()):
            if lab.get_fontsize() < 9:
                lab.set_fontsize(9)
        for tx in a.texts:
            if tx.get_fontsize() < 7.5:
                tx.set_fontsize(7.5)
        lo, hi = a.get_xlim()
        if 1990 < lo < 2030 and 2000 < hi < 2035:
            a.set_xticks([t for t in _MNL(integer=True).tick_values(lo, hi) if lo <= t <= 2025])
            a.set_xlim(lo, hi)
LISTS = ["analyses", "expression", "validation", "conditions"]
LABEL = {"PHYLOGENY": "Phylogeny", "GENE_STRUCTURE": "Gene structure", "MOTIF": "Conserved motifs", "DOMAIN": "Conserved domains",
         "PHYSICOCHEMICAL": "Physicochemical properties", "SUBCELLULAR_PREDICTION": "Subcellular prediction",
         "CHROMOSOME_LOCATION": "Chromosome location", "DUPLICATION": "Gene duplication", "SYNTENY": "Synteny/collinearity",
         "KAKS": "Ka/Ks", "CIS_ELEMENTS": "Promoter cis-elements", "PROTEIN_STRUCTURE": "Protein structure", "PPI": "Protein interaction network",
         "MIRNA": "miRNA targeting", "GO_KEGG": "GO/KEGG annotation", "COEXPRESSION": "Co-expression", "SNP_HAPLOTYPE": "SNP/haplotype",
         "EVOLUTION_DIVERGENCE": "Divergence/evolution", "RNA_SEQ": "RNA-seq", "QRT_PCR": "qRT-PCR", "MICROARRAY": "Microarray",
         "PUBLIC_DATABASE": "Public expression database", "OVEREXPRESSION": "Overexpression", "KNOCKOUT": "Knockout/mutant",
         "SILENCING": "Silencing (VIGS/RNAi)", "Y2H": "Interaction assay (Y2H…)", "SUBCELLULAR_LOCALIZATION": "Subcellular localization (exp.)",
         "TRANSACTIVATION": "Transactivation/binding", "HETEROLOGOUS": "Heterologous expression", "ENZYME_ASSAY": "Enzyme assay",
         "PROMOTER_REPORTER": "Promoter reporter", "DROUGHT": "Drought", "SALT": "Salt", "COLD": "Cold", "HEAT": "Heat",
         "OSMOTIC": "Osmotic", "WATERLOGGING": "Waterlogging", "HEAVY_METAL": "Heavy metal", "NUTRIENT": "Nutrient",
         "LIGHT": "Light", "OXIDATIVE": "Oxidative", "ABA": "ABA", "HORMONE_OTHER": "Other hormones", "FUNGAL": "Fungal",
         "BACTERIAL": "Bacterial", "VIRAL": "Viral", "INSECT/NEMATODE": "Insect/nematode", "TISSUES": "Tissues/organs",
         "DEVELOPMENT": "Development", "OTHER": "Other"}


def save(fig, name):
    _finish(fig)
    for ext in ("png", "svg"):
        for attempt in range(5):          # Windows sometimes holds the old file open briefly (preview/antivirus)
            try:
                fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight"); break
            except OSError:
                time.sleep(2)
    plt.close(fig)


def period(y):
    for a, b, lab in PERIODS:
        if a <= y <= b: return lab
    return None


def hbar(ax, labels, values, color=C[0], fmt="{:,}"):
    y = np.arange(len(labels))[::-1]
    ax.barh(y, values, color=color, height=0.72, edgecolor=SURF, linewidth=1)
    ax.set_yticks(y, labels); ax.grid(axis="y", visible=False); ax.tick_params(axis="y", length=0)
    for yi, v in zip(y, values):
        ax.text(v, yi, " " + fmt.format(v), va="center", fontsize=7.5, color=INK2)
    ax.set_xlim(0, max(values) * 1.15)


def main():
    d = pd.read_csv(GLOBAL_DATA_FILE, dtype=str, keep_default_na=False)
    a = d[d.in_analysis == "True"].copy()
    a["yr"] = pd.to_numeric(a.year, errors="coerce")
    a = a[a.yr.between(2000, 2025)].copy(); a["yr"] = a.yr.astype(int); a["period"] = a.yr.map(period)
    for f in LISTS:
        a[f] = a[f].map(lambda s: [x for x in s.split("; ") if x])
    a["cites"] = pd.to_numeric(a.cited_by, errors="coerce"); a["fwci_n"] = pd.to_numeric(a.fwci, errors="coerce")
    a["n_sp"] = pd.to_numeric(a.n_species_norm, errors="coerce")
    a["genes"] = pd.to_numeric(a.genes_total, errors="coerce")
    a["intl"] = a.international == "True"
    N = len(a); S = {"papers": N, "years": [int(a.yr.min()), int(a.yr.max())]}
    print("papers in analysis:", N)

    # 1 publications per year
    py = a.groupby("yr").size().reindex(range(2000, 2026), fill_value=0)
    py.rename("papers").to_frame().assign(cumulative=py.cumsum()).to_csv(TAB / "01_papers_per_year.csv")
    g = np.polyfit(py.index[py > 0] - 2000, np.log(py[py > 0]), 1)[0]
    S["doubling_time_years"] = round(np.log(2) / g, 1)
    S["papers_2020_2025"] = int(py.loc[2020:].sum()); S["share_2020_2025"] = round(py.loc[2020:].sum() / N * 100, 1)
    S["by_source"] = a.source.value_counts().to_dict()
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.bar(py.index, py.values, color=C[0], width=0.78, edgecolor=SURF, linewidth=0.8)
    ax.set_ylabel("Papers"); ax.grid(axis="x", visible=False); ax.set_xticks(range(2000, 2026, 5))
    for yr in (2010, 2015, 2020, 2025):
        ax.text(yr, py[yr], f"{py[yr]:,}", ha="center", va="bottom", fontsize=7.5, color=INK2)
    ax.set_title("Genome-wide gene family identification papers per year", loc="left", fontsize=10, color=INK)
    save(fig, "01_papers_per_year")

    # 2 countries
    ac = a[a.corr_country != ""]
    ct = ac.corr_country.value_counts()
    cdf = ct.rename("papers").to_frame().assign(share_pct=lambda x: (x.papers / len(ac) * 100).round(1))
    cdf.to_csv(TAB / "02_countries.csv")
    S["countries_n"] = int(ct.size); S["country_known"] = int(len(ac))
    S["top_countries"] = {k: [int(v), round(v / len(ac) * 100, 1)] for k, v in ct.head(10).items()}
    fig, ax = plt.subplots(figsize=(6, 4.2))
    top = ct.head(15); hbar(ax, list(top.index), list(top.values))
    ax.set_xscale("symlog", linthresh=100); ax.set_xlim(0, top.max() * 2.2)
    ax.set_xlabel("Papers (corresponding author; log scale above 100)")
    ax.set_title("Top 15 countries", loc="left", fontsize=10)
    save(fig, "02_countries_top15")
    # country share over time: top 4 + other
    top4 = list(ct.head(4).index)
    sh = (ac.assign(c=ac.corr_country.where(ac.corr_country.isin(top4), "Other"))
          .groupby(["yr", "c"]).size().unstack(fill_value=0))
    sh = sh.loc[2008:]; shp = sh.div(sh.sum(1), axis=0) * 100
    shp.round(1).to_csv(TAB / "02b_country_share_by_year.csv")
    fig, ax = plt.subplots(figsize=(7, 3.2))
    for i, c in enumerate(top4 + ["Other"]):
        col = C[i] if c != "Other" else OTHER
        ax.plot(shp.index, shp[c], color=col, lw=2, label=c)
    ax.text(2025.3, shp["China"].iloc[-1], f'China {shp["China"].iloc[-1]:.0f}%', color=INK, fontsize=7.5, va="center")
    ax.set_ylabel("Share of papers (%)"); ax.set_xlim(2008, 2027.5); ax.set_ylim(0, 100)
    ax.legend(loc="upper left", ncol=5, fontsize=7.5)
    ax.set_title("Country share of papers per year (2008–2025)", loc="left", fontsize=10)
    save(fig, "02b_country_share_by_year")

    # 2c world map (corresponding-author country; log colour scale because China holds ~80%)
    import plotly.graph_objects as go, pycountry
    mp = ac.groupby("corr_country_code").agg(country=("corr_country", "first"), papers=("key", "size")).reset_index()
    mp["iso3"] = [getattr(pycountry.countries.get(alpha_2=c), "alpha_3", None) for c in mp.corr_country_code]
    mp = mp.dropna(subset=["iso3"]); mp.to_csv(TAB / "02c_map_countries.csv", index=False)
    ticks = [1, 10, 100, 1000, 7000]
    fig_m = go.Figure(go.Choropleth(
        locations=mp.iso3, z=np.log10(mp.papers), text=mp.country, customdata=mp.papers,
        hovertemplate="%{text}: %{customdata:,} papers<extra></extra>",
        colorscale=[[0, "#dbe9f8"], [0.5, "#5b9be6"], [1, "#0d3f80"]], zmin=0, zmax=np.log10(8000),
        marker_line_color="#ffffff", marker_line_width=0.4,
        colorbar=dict(title=dict(text="Papers", font=dict(size=12)), tickvals=np.log10(ticks), ticktext=[f"{t:,}" for t in ticks],
                      thickness=12, len=0.6, outlinewidth=0)))
    fig_m.update_geos(projection_type="natural earth", showframe=False, showcoastlines=False, showcountries=True,
                      countrycolor="#ffffff", showland=True, landcolor="#eeeeec", showocean=False, lataxis_range=[-58, 85])
    fig_m.update_layout(margin=dict(l=0, r=0, t=0, b=0), width=1100, height=560, paper_bgcolor=SURF,
                        font=dict(family="DejaVu Sans, Arial", color=INK))
    fig_m.write_image(FIG / "02c_world_map.png", scale=3); fig_m.write_image(FIG / "02c_world_map.svg")
    

    # 3 collaboration
    col = ac.groupby("yr").agg(papers=("key", "size"), international=("intl", "sum"))
    col["pct_international"] = (col.international / col.papers * 100).round(1)
    col.to_csv(TAB / "03_collaboration_by_year.csv")
    S["international_pct"] = round(ac.intl.mean() * 100, 1)
    bycountry = ac.groupby("corr_country").agg(papers=("key", "size"), pct_international=("intl", "mean"))
    bycountry["pct_international"] = (bycountry.pct_international * 100).round(1)
    bycountry.sort_values("papers", ascending=False).to_csv(TAB / "03b_collaboration_by_country.csv")
    pairs = {}
    for s in ac.all_countries:
        cs = sorted({x.strip() for x in s.split(";") if x.strip()})
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)): pairs[(cs[i], cs[j])] = pairs.get((cs[i], cs[j]), 0) + 1
    pd.Series(pairs).sort_values(ascending=False).rename("co_papers").rename_axis(["country_a", "country_b"]).to_csv(TAB / "03c_country_pairs.csv")
    S["top_pairs"] = [[f"{p[0]}–{p[1]}", int(v)] for p, v in sorted(pairs.items(), key=lambda x: -x[1])[:5]]
    fig, ax = plt.subplots(figsize=(7, 2.8))
    cc = col.loc[2005:]
    ax.plot(cc.index, cc.pct_international, color=C[0], lw=2)
    ax.set_ylabel("International papers (%)"); ax.set_ylim(0, max(cc.pct_international) * 1.25)
    ax.set_title("Share of papers with authors from more than one country", loc="left", fontsize=10)
    save(fig, "03_international_by_year")

    # 4 journals + Bradford
    a["jn"] = a.journal.str.strip().str.lower().str.replace(r"\s*\(.*?\)\s*", " ", regex=True).str.replace(r"[^a-z0-9& ]", "", regex=True).str.strip()
    disp = a.groupby("jn").journal.agg(lambda s: s.value_counts().index[0])
    jt = a[a.jn != ""].jn.value_counts()
    jdf = pd.DataFrame({"journal": disp[jt.index].values, "papers": jt.values})
    jdf["cum_papers"] = jdf.papers.cumsum(); third = jdf.papers.sum() / 3
    jdf["bradford_zone"] = np.where(jdf.cum_papers <= third, 1, np.where(jdf.cum_papers <= 2 * third, 2, 3))
    jdf.to_csv(TAB / "04_journals_bradford.csv", index=False)
    zones = jdf.groupby("bradford_zone").agg(journals=("journal", "size"), papers=("papers", "sum"))
    zones.to_csv(TAB / "04b_bradford_zones.csv")
    S["journals_n"] = int(len(jdf)); S["bradford_zones"] = zones.journals.astype(int).tolist()
    S["top_journals"] = jdf.head(5)[["journal", "papers"]].values.tolist()
    fig, ax = plt.subplots(figsize=(6, 4.2))
    tj = jdf.head(15); hbar(ax, [j[:46] for j in tj.journal], list(tj.papers))
    ax.set_xlabel("Papers"); ax.set_title("Top 15 journals", loc="left", fontsize=10)
    save(fig, "04_journals_top15")

    # 5 species
    sp = a[["key", "yr", "species_norm"]].assign(s=a.species_norm.str.split(", ")).explode("s")
    sp = sp[~sp.s.isin(["Not stated", ""]) & sp.s.notna()]
    st = sp.s.value_counts()
    st.rename("papers").to_csv(TAB / "05_species.csv")
    sp["genus"] = sp.s.str.split().str[0].str.strip("[]")
    sp.genus.value_counts().rename("papers").to_csv(TAB / "05b_genera.csv")
    S["species_distinct"] = int(st.size); S["genera_distinct"] = int(sp.genus.nunique())
    S["top_species"] = [[k, int(v)] for k, v in st.head(10).items()]
    ns = a.dropna(subset=["n_sp"])
    S["single_species_pct"] = round((ns.n_sp == 1).mean() * 100, 1)
    nsp = ns.groupby("period").n_sp.agg(papers="size", mean="mean", median="median",
                                         multi_pct=lambda x: (x > 1).mean() * 100).round(2)
    nsp.to_csv(TAB / "05c_species_per_paper_by_period.csv")
    newsp = sp.sort_values("yr").drop_duplicates("s").groupby("yr").size().reindex(range(2000, 2026), fill_value=0)
    newsp.rename("new_species").to_frame().assign(cumulative=newsp.cumsum()).to_csv(TAB / "05d_new_species_per_year.csv")
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.4), gridspec_kw={"width_ratios": [1.1, 1]})
    ts = st.head(20); hbar(axs[0], [f"$\\it{{{s.replace(' ', '~')}}}$" if not s.startswith("[") else s for s in ts.index], list(ts.values))
    axs[0].set_xlabel("Papers"); axs[0].set_title("Top 20 species", loc="left", fontsize=10)
    cum = newsp.cumsum()
    axs[1].plot(cum.index, cum.values, color=C[0], lw=2)
    axs[1].set_ylabel("Distinct species studied (cumulative)")
    axs[1].set_title("Species coverage over time", loc="left", fontsize=10)
    axs[1].text(2025, cum.iloc[-1], f" {cum.iloc[-1]:,}", va="center", fontsize=8, color=INK2)
    save(fig, "05_species")

    # 6 gene families + class
    fm = a[["key", "yr", "period", "gene_family", "family_class"]].assign(f=a.gene_family.str.split(r",\s*")).explode("f")
    fm = fm[fm.f.fillna("").str.strip().ne("") & fm.f.ne("Not stated")]
    fm["k"] = fm.f.map(family_key)
    display = fm.groupby("k").f.agg(lambda s: s.value_counts().index[0]).to_dict()   # most common spelling per key
    fg = [family_norm(f, c, display) for f, c in zip(fm.f, fm.family_class)]
    fm["family"] = [x[0] for x in fg]; fm["group"] = [x[1] for x in fg]
    fm[["key", "f", "family", "group"]].rename(columns={"f": "name_in_abstract"}).to_csv(TAB / "06x_family_name_mapping.csv", index=False)
    fam = fm.drop_duplicates(["key", "family"]).family.value_counts()
    fam.rename("papers").rename_axis("family").to_csv(TAB / "06a_families_detailed.csv")
    ft = fm.drop_duplicates(["key", "group"]).group.value_counts()
    fdf = pd.DataFrame({"family": ft.index, "papers": ft.values})
    fdf.to_csv(TAB / "06_gene_families.csv", index=False)
    S["family_names_raw"] = int(fm.f.str.strip().nunique()); S["families_detailed"] = int(len(fam))
    S["families_distinct"] = int(len(fdf)); S["top_families"] = fdf.head(10).values.tolist()
    S["families_studied_once"] = int((ft == 1).sum())
    first_fam = {k: g for k, g in fm.drop_duplicates("key")[["key", "group"]].values}
    cls = a[a.family_class != ""].groupby(["period", "family_class"]).size().unstack(fill_value=0)
    clsp = (cls.div(cls.sum(1), axis=0) * 100).round(1)
    cls.to_csv(TAB / "06b_family_class_by_period.csv"); clsp.to_csv(TAB / "06c_family_class_pct_by_period.csv")
    S["family_class"] = a.family_class.value_counts().to_dict()
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.8), gridspec_kw={"width_ratios": [1, 1.1]})
    tf = fdf.head(25); hbar(axs[0], list(tf.family), list(tf.papers))
    axs[0].set_xlabel("Papers"); axs[0].set_title("Top 25 gene families (subfamilies merged)", loc="left", fontsize=10)
    order = ["TF", "enzyme", "transporter/channel", "kinase/receptor", "structural/other protein", "mixed", "unclear"]
    order = [o for o in order if o in clsp.columns]
    left = np.zeros(len(clsp))
    for i, o in enumerate(order):
        colr = C[i] if o not in ("unclear", "mixed") else (OTHER if o == "unclear" else "#cfcec8")
        axs[1].barh(range(len(clsp))[::-1], clsp[o], left=left, color=colr, edgecolor=SURF, linewidth=1.2, label=o, height=0.6)
        left += clsp[o].values
    axs[1].set_yticks(range(len(clsp))[::-1], clsp.index); axs[1].set_xlim(0, 100); axs[1].grid(axis="y", visible=False)
    axs[1].set_xlabel("Share of papers (%)"); axs[1].legend(ncol=2, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.14))
    axs[1].set_title("Gene family class by period", loc="left", fontsize=10)
    save(fig, "06_gene_families")

    # 7-8 methods, expression, validation, conditions over periods (heatmaps of % of papers)
    def share_table(field, fname, title, fig_name):
        codes = sorted({c for l in a[field] for c in l})
        rows = {}
        for p in [x[2] for x in PERIODS]:
            sub = a[a.period == p]
            rows[p] = {c: sub[field].map(lambda l: c in l).mean() * 100 for c in codes}
        t = pd.DataFrame(rows).round(1)
        t["all"] = [round(a[field].map(lambda l: c in l).mean() * 100, 1) for c in t.index]
        t = t.sort_values("all", ascending=False)
        t.insert(0, "label", [LABEL.get(c, c) for c in t.index])
        t.to_csv(TAB / fname)
        m = t[[x[2] for x in PERIODS]].values
        fig, ax = plt.subplots(figsize=(6, 0.32 * len(t) + 1.2))
        im = ax.imshow(m, cmap=SEQ, aspect="auto", vmin=0, vmax=max(60, m.max()))
        ax.set_xticks(range(4), [x[2] for x in PERIODS]); ax.set_yticks(range(len(t)), t.label)
        ax.grid(False); ax.tick_params(length=0)
        for s in ax.spines.values(): s.set_visible(False)
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                ax.text(j, i, f"{m[i, j]:.0f}", ha="center", va="center", fontsize=7.5,
                        color=SURF if m[i, j] > 0.55 * max(60, m.max()) else INK)
        ax.set_title(title, loc="left", fontsize=10)
        cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("% of papers", fontsize=8); cb.outline.set_visible(False)
        save(fig, fig_name)
        return t
    ta = share_table("analyses", "07_analyses_by_period.csv", "In-silico analyses named in the abstract (% of papers)", "07_analyses_by_period")
    te = share_table("expression", "08_expression_by_period.csv", "Expression methods named in the abstract (% of papers)", "08_expression_by_period")
    tv = share_table("validation", "09_validation_by_period.csv", "Functional validation named in the abstract (% of papers)", "09_validation_by_period")
    tc = share_table("conditions", "10_conditions_by_period.csv", "Conditions / contexts studied (% of papers)", "10_conditions_by_period")
    S["analyses_pct"] = ta["all"].head(8).to_dict(); S["expression_pct"] = te["all"].to_dict()
    S["validation_any_pct"] = round(a.validation.map(bool).mean() * 100, 1)
    S["validation_any_by_period"] = a.groupby("period").validation.apply(lambda s: round(s.map(bool).mean() * 100, 1)).to_dict()
    S["conditions_pct"] = tc["all"].head(10).to_dict()
    # expression methods per year (line)
    ey = a.groupby("yr").apply(lambda s: pd.Series({c: s.expression.map(lambda l: c in l).mean() * 100
                                                    for c in ["QRT_PCR", "RNA_SEQ", "MICROARRAY"]})).loc[2005:]
    ey.round(1).to_csv(TAB / "08b_expression_by_year.csv")
    fig, ax = plt.subplots(figsize=(7, 3))
    for i, c in enumerate(["QRT_PCR", "RNA_SEQ", "MICROARRAY"]):
        ax.plot(ey.index, ey[c], color=C[i], lw=2, label=LABEL[c])
        ax.text(2025.3, ey[c].iloc[-1], LABEL[c], fontsize=7.5, va="center", color=INK)
    ax.set_xlim(2005, 2028); ax.set_ylabel("% of papers"); ax.legend(loc="upper left", fontsize=7.5, ncol=3)
    ax.set_title("Expression methods named in abstracts, per year", loc="left", fontsize=10)
    save(fig, "08b_expression_by_year")

    # 11 results: gene counts, groups, duplication, selection
    gn = a.dropna(subset=["genes"]); gn = gn[gn.genes > 0]
    gs = gn.groupby("period").genes.describe(percentiles=[.25, .5, .75]).round(1)
    gs.to_csv(TAB / "11_genes_total_by_period.csv")
    gcl = gn.groupby("family_class").genes.describe(percentiles=[.25, .5, .75]).round(1); gcl.to_csv(TAB / "11b_genes_total_by_class.csv")
    gfam = gn[gn.n_sp == 1].assign(group=lambda x: x.key.map(first_fam))
    gf = gfam[gfam.group.isin(fdf.family.head(25))].groupby("group").genes.agg(
        papers="size", median="median", q1=lambda x: x.quantile(.25), q3=lambda x: x.quantile(.75)).round(1)
    gf.sort_values("median", ascending=False).to_csv(TAB / "11c_genes_single_species_top_families.csv")
    S["genes_total_median"] = float(gn.genes.median()); S["genes_stated_pct"] = round(len(gn) / N * 100, 1)
    dup = a[a.duplication_norm.ne("Not stated") & a.duplication_norm.ne("")]
    dk = dup.duplication_norm.str.split(", ").explode().value_counts()
    dk.rename("papers").to_frame().assign(pct_of_reporting=lambda x: (x.papers / len(dup) * 100).round(1)).to_csv(TAB / "11d_duplication_types.csv")
    first = dup.duplication_norm.str.split(", ").str[0].value_counts()
    first.rename("papers_main_type").to_csv(TAB / "11e_duplication_main_type.csv")
    sel = a[a.selection.isin(["Purifying", "Positive", "Both"])].selection.value_counts()
    sel.rename("papers").to_csv(TAB / "11f_selection.csv")
    grp = a.groups.str.extract(r"^(\d+)\s+(\w+)"); grp.columns = ["n", "word"]; grp = grp.dropna(); grp["n"] = grp.n.astype(int)
    grp.n.describe().round(1).to_csv(TAB / "11g_groups_number.csv")
    grp.word.str.lower().value_counts().head(10).rename("papers").to_csv(TAB / "11h_groups_word.csv")
    S["duplication_reported_pct"] = round(len(dup) / N * 100, 1); S["duplication_types"] = dk.to_dict()
    S["duplication_main"] = first.to_dict(); S["selection"] = sel.to_dict()
    S["groups_median"] = float(grp.n.median())
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
    bins = np.logspace(0, np.log10(max(gn.genes.max(), 10)), 40)
    axs[0].hist(gn.genes, bins=bins, color=C[0], edgecolor=SURF, linewidth=0.6); axs[0].set_xscale("log")
    axs[0].axvline(gn.genes.median(), color=INK, lw=1, ls="--")
    axs[0].text(gn.genes.median(), axs[0].get_ylim()[1] * 0.92, f" median {gn.genes.median():.0f}", fontsize=7.5, color=INK)
    axs[0].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    axs[0].set_xlabel("Family members identified (log scale)"); axs[0].set_ylabel("Papers")
    axs[0].set_title("Size of identified gene families", loc="left", fontsize=10)
    dd = dk.reindex(["Segmental", "Tandem", "WGD", "Dispersed", "Proximal", "Transposed"]).dropna()
    hbar(axs[1], list(dd.index), list(dd.astype(int).values))
    axs[1].set_xlabel(f"Papers (of {len(dup):,} reporting a type)"); axs[1].set_title("Duplication types reported", loc="left", fontsize=10)
    save(fig, "11_gene_counts_duplication")

    # 12 citations
    cit = a.dropna(subset=["cites"])
    S["citations_total"] = int(cit.cites.sum()); S["citations_median"] = float(cit.cites.median())
    S["fwci_median"] = round(float(cit.fwci_n.median()), 2)
    S["top10pct_share"] = round((a.top10pct == "True").mean() * 100, 1)
    cc = cit[cit.corr_country != ""].groupby("corr_country").agg(papers=("key", "size"), citations=("cites", "sum"),
                                                                 median_citations=("cites", "median"), median_fwci=("fwci_n", "median"),
                                                                 pct_top10=("top10pct", lambda x: (x == "True").mean() * 100))
    cc = cc[cc.papers >= 30].round(2).sort_values("median_fwci", ascending=False); cc.to_csv(TAB / "12_citations_by_country.csv")
    cy = cit.groupby("yr").agg(papers=("key", "size"), median_citations=("cites", "median"), median_fwci=("fwci_n", "median")).round(2)
    cy.to_csv(TAB / "12b_citations_by_year.csv")
    cit.nlargest(20, "cites")[["key", "year", "title", "journal", "corr_country", "cites"]].to_csv(TAB / "12c_most_cited.csv", index=False)
    S["most_cited"] = cit.nlargest(3, "cites")[["year", "title", "cites"]].values.tolist()
    fig, ax = plt.subplots(figsize=(6.5, 0.27 * len(cc) + 1.2))
    y = np.arange(len(cc))[::-1]
    ax.scatter(cc.median_fwci, y, s=20 + 6 * np.sqrt(cc.papers), color=C[0], edgecolor=SURF, linewidth=1.5, zorder=3)
    ax.axvline(1, color=INK2, lw=0.8, ls="--"); ax.text(1.03, -0.6, "world average (FWCI = 1)", fontsize=7, color=INK2, va="center")
    ax.set_xlim(0.8, cc.median_fwci.max() * 1.08); ax.set_ylim(-1, len(cc))
    ax.set_yticks(y, [f"{c} ({n:,})" for c, n in zip(cc.index, cc.papers)]); ax.grid(axis="y", visible=False); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Median field-weighted citation impact"); ax.set_title("Citation impact by country (≥ 30 papers; dot size = papers)", loc="left", fontsize=10)
    save(fig, "12_citation_impact_by_country")

    (HERE / "summary.json").write_text(json.dumps(S, indent=1, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps({k: S[k] for k in ["papers", "doubling_time_years", "share_2020_2025", "countries_n", "international_pct",
                                        "journals_n", "bradford_zones", "species_distinct", "families_distinct",
                                        "validation_any_pct", "genes_total_median", "fwci_median"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
