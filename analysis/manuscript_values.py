import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR 
import collections
import re
import sys

import numpy as np
import pandas as pd

DATA = str(DATA_FILE)
d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
N = len(d)
d["y"] = d.Year.astype(int)
d["g"] = pd.to_numeric(d["Number of Genes Identified"])
PERIODS = ["Before 2013", "2013-2016", "2017-2020", "2021-2025"]
d["period"] = pd.cut(d.y, [1990, 2012, 2016, 2020, 2030], labels=PERIODS)


def pc(n, base=N, dp=1):
    return f"{n / base * 100:.{dp}f}%"


def np_(n, base=N):
    return f"{pc(n, base)} (n = {n})"


V = {"N": N}
# countries
c = d["Country of Origin"].value_counts()
V["countries"] = c
V["n_countries"] = len(c)
V["top3"] = int(c[["China", "India", "Pakistan"]].sum())
V["country_one"] = sorted(c[c == 1].index)
V["country_two"] = sorted(c[c == 2].index)
# years
a = d.y.value_counts()
V["annual"] = a
V["y0610"] = int(a.reindex(range(2006, 2011), fill_value=0).sum())
V["y1213"] = int(a.reindex([2012, 2013], fill_value=0).sum())
V["y2019"] = int(a[2019])
V["y2225"] = int(a.reindex(range(2022, 2026), fill_value=0).sum())
V["y2025"] = int(a[2025])
V["first_year"] = int(d.y.min())
# teams
au = d.Authors.map(lambda s: len([x for x in re.split(r",|;", s) if x.strip()]))
V["team_mean"], V["team_median"], V["team_max"], V["team_min"] = au.mean(), au.median(), int(au.max()), int(au.min())
V["team"] = [int(((au >= l) & (au <= h)).sum()) for l, h in ((1, 3), (4, 6), (7, 10), (11, 999))]
# journals
j = d.Journal.value_counts()
V["journals"] = j
V["n_journals"] = len(j)
V["top5j"] = int(j[["International Journal of Molecular Sciences", "BMC Genomics", "Plants", "BMC Plant Biology", "Frontiers in Plant Science"]].sum())
named = ["International Journal of Molecular Sciences", "BMC Genomics", "Plants", "BMC Plant Biology", "Frontiers in Plant Science", "Scientific Reports", "PLOS ONE", "Genes"]
V["rest_papers"] = N - int(j[named].sum())
V["rest_journals"] = len(j) - len(named)
V["rest_max"] = int(j.drop(named).max())
V["top5_is_top5"] = list(j.index[:5])
# organisms: each species once per paper
sp = lambda v: [x.strip() for x in v.split(",") if x.strip() and x.strip() != "Not Available"]
d["P"] = d["Primary Organism"].map(sp)
d["S"] = [[x for x in sp(s) if x not in p] for s, p in zip(d["Secondary Organism"], d.P)]
P = collections.Counter(x for r in d.P for x in r)
S = collections.Counter(x for r in d.S for x in r)
V["P"], V["S"] = P, S
V["n_primary"] = len(P)
V["with_ref"] = int((d.S.map(len) > 0).sum())
top5sp = [k for k, _ in sorted(P.items(), key=lambda kv: (-kv[1], kv[0]))[:5]]
V["top5sp"] = top5sp
V["top5sp_union"] = int(d.P.map(lambda r: any(x in top5sp for x in r)).sum())
V["wheat_rice_union"] = int(d.P.map(lambda r: "Triticum aestivum" in r or "Oryza sativa" in r).sum())
V["P_rank"] = sorted(P.items(), key=lambda kv: (-kv[1], kv[0]))[:16]
# families
d["fams"] = d["Gene Family"].map(lambda v: [x.strip() for x in re.split(r",(?![^(]*\))", v) if x.strip()])
F = collections.Counter(x for r in d.fams for x in r)
V["F"] = F
V["n_family_names"] = len(F)
V["n_family_entries"] = d["Gene Family"].nunique()
V["family_single"] = sum(v == 1 for v in F.values())
MAPK = {"MAPK", "MAPKK", "MAPKKK"}
SHORT = {"SPL (SBP-box)": "SPL"}
d["fam7"] = d.fams.map(lambda r: list(dict.fromkeys("MAPK cascade" if x in MAPK else SHORT.get(x, x) for x in r)))
F7 = collections.Counter(x for r in d.fam7 for x in r)
V["F7"] = F7
six = sorted([f for f, v in F7.items() if v >= 5], key=lambda f: (-F7[f], f))
V["six"] = six
V["fam_first"] = {f: int(d[d.fam7.map(lambda r: f in r)].y.min()) for f in list(six) + ["bZIP", "SWEET", "SOD", "GST"]}
trio = ["DOF", "MAPK cascade", "WRKY"]
V["trio_before2020"] = int(d[d.fam7.map(lambda r: any(f in r for f in trio)) & (d.y < 2020)].shape[0])
V["bzip_spl"] = int(d[d.fam7.map(lambda r: "bZIP" in r or "SPL" in r)].shape[0])
V["bzip_spl_first"] = int(d[d.fam7.map(lambda r: "bZIP" in r or "SPL" in r)].y.min())
sel = d[d.fam7.map(lambda r: any(f in r for f in six))]
V["six_papers"] = len(sel)
V["six_2225"] = int((sel.y >= 2022).sum())
V["six_before"] = int((sel.y < 2022).sum())
# genes
V["genes"] = {p: dict(n=int((d.period == p).sum()), mean=d[d.period == p].g.mean(), median=d[d.period == p].g.median(), over100=int((d[d.period == p].g > 100).sum()), max=int(d[d.period == p].g.max())) for p in PERIODS}
V["genes_all_median"] = d.g.median()
V["genes_le30"] = int((d.g <= 30).sum())


# sources (Figure 9 rule)
def source(v):
    l = v.lower()
    for key, name in [("phytozome", "Phytozome"), ("ensembl", "Ensembl Plants"), ("gramene", "Ensembl Plants"), ("ncbi", "NCBI"), ("genbank", "NCBI"), ("refseq", "NCBI"), ("cotton", "Cotton databases"), ("sol genomics", "Sol Genomics Network"), ("sgn", "Sol Genomics Network"),
                      ("gdr", "GDR"), ("rosaceae", "GDR"), ("brad", "BRAD"), ("brassica", "BRAD"), ("cruciferae", "BRAD"), ("rice genome annotation", "RGAP"), ("rgap", "RGAP"), ("msu", "RGAP"), ("tigr", "RGAP")]:
        if key in l:
            return name
    return "Other databases"


d["source"] = d["Genome Source"].map(source)
V["sources"] = d.source.value_counts()
V["tair_any"] = int(d["Genome Source"].str.contains("TAIR", case=False).sum())
V["org_src"] = {o: d[d.P.map(lambda r: o in r)].source.value_counts().head(3).to_dict() for o in ["Triticum aestivum", "Oryza sativa", "Solanum lycopersicum", "Crocus sativus", "Gossypium hirsutum", "Medicago truncatula"]}
q, r = d["Experimental Validation (qPCR)"] == "Yes", d["Expression Analysis (RNA-seq)"] == "Yes"
V["pipe_both_by_org"] = {o: (int((d.P.map(lambda x: o in x) & q & r).sum()), int(d.P.map(lambda x: o in x).sum())) for o in ["Triticum aestivum", "Oryza sativa", "Solanum lycopersicum"]}
# tools
b = d["BLAST (E-value cut-off used)"]
V["blast"] = int((b != "Not Used").sum())
V["blast_ns"] = int((b == "BLAST (not specified)").sum())
V["blast_stated"] = V["blast"] - V["blast_ns"]
V["blast_vals"] = b.value_counts()
h = d["HMMER (with Pfam domains ID)"]
V["hmmer"] = int((~h.str.startswith("No")).sum())
V["hmmer_id"] = int(h.str.contains(r"PF\d").sum())
for k, col in {"smart": "SMART", "cdd": "CDD", "meme": "MEME", "mega": "MEGA", "iq": "IQ-TREE", "mcscanx": "MCScanX", "tbtools": "TBTools"}.items():
    V[k] = int((d[col] == "Yes").sum())
ALLT = ["Other (Tools for Domain Analysis)", "Others (Tools for Motif Analysis)", "Other (Tools for Synteny)", "Other (Tools for Visualization)", "Other Tools Used In Methodology"]
has = lambda word, cols=ALLT: d[cols].apply(lambda row: any(re.search(rf"\b{word}\b", v, re.I) for v in row), axis=1)
V["pfam"] = int(has("Pfam", ["Other (Tools for Domain Analysis)"]).sum())
V["weblogo"] = int(has("WebLogo").sum())
V["pgdd"] = int(has("Plant Genome Duplication Database", ["Other (Tools for Synteny)"]).sum())
V["gsds"] = int(has("GSDS", ["Other (Tools for Visualization)"]).sum())
V["expasy"] = int(has("ExPASy", ["Other Tools Used In Methodology"]).sum())
V["protparam_any"] = int(has("ProtParam").sum())
V["wolf"] = int(has("WoLF PSORT", ["Other Tools Used In Methodology"]).sum())
V["clustalw"] = int(has("ClustalW", ["Other Tools Used In Methodology"]).sum())
V["clustalx"] = int(has("ClustalX", ["Other Tools Used In Methodology"]).sum())
# phylogeny (Figure 10 rule)


def stype(t):
    l = " ".join(str(t).replace("\xa0", " ").split()).lower().replace("-", "")
    for key, name in (("subfamil", "Subfamily"), ("clade", "Clade"), ("cluster", "Cluster"), ("subgroup", "Subgroup"), ("group", "Group"), ("class", "Class"), ("type", "Type")):
        if key in l:
            return name
    return "Not specified" if "not specified" in l or not l else "Other"


st = d["Phylogenetic Structure"].map(stype)
V["phylo"] = st.value_counts()
V["phylo_n"] = int((st != "Not specified").sum())
# duplication
col = d["Gene Duplication Events (TANDEM/SEGMENTAL)"]
V["dup"] = col.value_counts()
V["dup_ns"] = int((col == "Not specified").sum())
V["dup_unnamed"] = int((col == "Type not named").sum())
V["dup_named"] = N - V["dup_ns"] - V["dup_unnamed"]
V["dup_st"] = int((col == "Segmental, Tandem").sum())
V["dup_s"] = int((col == "Segmental").sum())
V["dup_t"] = int((col == "Tandem").sum())
V["dup_other"] = V["dup_named"] - V["dup_st"] - V["dup_s"] - V["dup_t"]
V["dup_seg_any"] = int(col.str.contains("Segmental").sum())
V["dup_tan_any"] = int(col.str.contains("Tandem").sum())
V["dup_wgd_any"] = int(col.str.contains("WGD").sum())
V["dup_st_2325"] = int(((col == "Segmental, Tandem") & (d.y >= 2023)).sum())
V["dup_st_by_year"] = d[col == "Segmental, Tandem"].y.value_counts().sort_index().to_dict()
V["dup_st_maxg"] = int(d[col == "Segmental, Tandem"].g.max())
V["dup_st_over300"] = int((d[col == "Segmental, Tandem"].g > 300).sum())
nm = ~col.isin(["Not specified", "Type not named"])
V["dup_named_le30"] = int((d[nm].g <= 30).sum())
# validation
V["qpcr"], V["rnaseq"] = int(q.sum()), int(r.sum())
V["both"], V["qonly"], V["ronly"], V["silico"] = int((q & r).sum()), int((q & ~r).sum()), int((~q & r).sum()), int((~q & ~r).sum())
V["both_2025"] = int(((q & r) & (d.y == 2025)).sum())
V["both_by_year"] = d[q & r].y.value_counts().sort_index().to_dict()
V["microarray"] = int((d["Expression (microarray/EST)"] == "Yes").sum())
V["silico_micro"] = int(((~q & ~r) & (d["Expression (microarray/EST)"] == "Yes")).sum())

if __name__ == "__main__":
    for k, v in V.items():
        if isinstance(v, pd.Series):
            v = v.to_dict()
        print(k, "=", v)
