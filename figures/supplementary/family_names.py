import re

# key -> canonical family. Keys are the normalised form from key(): upper case, no spaces/hyphens/slashes, no "family"/"genes".
SYN = {
    "HSF": "HSF", "GRF": "GRF", "GROWTHREGULATINGFACTOR": "GRF", "MADSBOX": "MADS-box", "MADS": "MADS-box", "AP2ERF": "AP2/ERF", "AP2EREBP": "AP2/ERF", "AUXIAA": "Aux/IAA",
    "HDZIP": "HD-Zip", "NFY": "NF-Y", "1433": "14-3-3", "SNRK2": "SnRK2", "SNRK": "SnRK",
    "NBSLRR": "NLR", "NBLRR": "NLR", "NLR": "NLR", "NBS": "NLR", "NBARC": "NLR", "CNL": "NLR", "TIRNBSLRR": "NLR",
    "CYP450": "CYP450", "P450": "CYP450", "CYP": "CYP450", "CYTOCHROMEP450": "CYP450",
    "SPL": "SPL (SBP-box)", "SBP": "SPL (SBP-box)", "SBPBOX": "SPL (SBP-box)",
    "HAKKUPKT": "HAK/KUP/KT", "KTHAKKUP": "HAK/KUP/KT", "KUPHAKKT": "HAK/KUP/KT", "HAK": "HAK/KUP/KT", "KUP": "HAK/KUP/KT",
    "C2H2": "C2H2-ZFP", "C2H2ZFP": "C2H2-ZFP", "C2H2ZF": "C2H2-ZFP", "C2H2ZINCFINGER": "C2H2-ZFP",
    "JMJC": "JmjC", "JMJ": "JmjC", "SHSP": "HSP20", "HSP20": "HSP20", "SMALLHEATSHOCKPROTEIN": "HSP20",
    "ZFHD": "ZF-HD", "ZHD": "ZF-HD", "BES1": "BES1/BZR1", "BZR": "BES1/BZR1", "BES1BZR1": "BES1/BZR1", "BZR1": "BES1/BZR1",
    "PYL": "PYL", "PYRPYL": "PYL", "PYRPYLRCAR": "PYL", "CDPK": "CDPK", "CPK": "CDPK", "MAPK": "MAPK", "MPK": "MAPK",
    "MAPKK": "MAPKK", "MKK": "MAPKK", "MAPKKK": "MAPKKK", "MAP3K": "MAPKKK", "EIN3EIL": "EIN3/EIL", "EIL": "EIN3/EIL",
    "FAR1": "FAR1/FHY3", "FAR1FHY3": "FAR1/FHY3", "WAK": "WAK/WAKL", "WAKL": "WAK/WAKL", "WAKWAKL": "WAK/WAKL",
    "RING": "RING", "RINGFINGER": "RING", "C3HC4TYPERINGFINGER": "RING", "ARGONAUTE": "AGO", "AGO": "AGO",
    "NPF": "NPF (NRT1/PTR)", "NRT1": "NPF (NRT1/PTR)", "NRT1PTR": "NPF (NRT1/PTR)", "AMT1": "AMT", "LHCB": "LHC",
    "YUCCA": "YUCCA", "YUC": "YUCCA", "LTP": "nsLTP", "NSLTP": "nsLTP", "SDG": "SET domain", "SET": "SET domain",
    "SETDOMAIN": "SET domain", "GSK": "GSK3", "GSK3": "GSK3", "SKP1": "SKP1", "SKP1LIKE": "SKP1", "ALKB": "ALKBH",
    "ROP": "Rop", "PHD": "PHD finger", "PHDFINGER": "PHD finger", "CESA": "CesA/Csl", "CESACSL": "CesA/Csl", "CSL": "CesA/Csl",
    "EXPANSIN": "Expansin", "EXP": "Expansin", "LACCASE": "Laccase", "LAC": "Laccase", "AQUAPORIN": "Aquaporin (MIP)",
    "AQP": "Aquaporin (MIP)", "MIP": "Aquaporin (MIP)", "CHITINASE": "Chitinase", "INVERTASE": "Invertase", "INV": "Invertase",
    "ANNEXIN": "Annexin", "ANN": "Annexin", "OLEOSIN": "Oleosin", "DEHYDRIN": "Dehydrin", "DHN": "Dehydrin",
    "CYSTATIN": "Cystatin", "MYOSIN": "Myosin", "KINESIN": "Kinesin", "TUBULIN": "Tubulin", "CALEOSIN": "Caleosin",
    "PHYTOCYANIN": "Phytocyanin", "METACASPASE": "Metacaspase", "FORMIN": "Formin", "LECTIN": "Lectin",
    "CLASSIIIPEROXIDASE": "Class III peroxidase", "PRX": "Class III peroxidase", "POD": "Class III peroxidase",
    "PEROXIDASE": "Class III peroxidase", "TRIHELIX": "Trihelix", "HOMEOBOX": "Homeobox", "CYCLIN": "Cyclin",
    "ETHYLENERECEPTOR": "Ethylene receptor", "SUBTILASE": "Subtilase (SBT)", "SBT": "Subtilase (SBT)",
    "PTYPEATPASE": "P-type ATPase", "DNAJ": "DnaJ (HSP40)", "HSP40": "DnaJ (HSP40)", "FBOX": "F-box", "UBOX": "U-box",
    "HSP": "HSP (unspecified)", "PR10": "PR-10", "PR1": "PR-1", "GLYI": "Glyoxalase I", "GLYII": "Glyoxalase II",
    "LIPOXYGENASE": "LOX", "LOX": "LOX", "SUCROSESYNTHASE": "SUS", "SUS": "SUS", "PHENYLALANINEAMMONIALYASE": "PAL",
    "TERPENESYNTHASE": "TPS", "TPS": "TPS", "GLUTATHIONESTRANSFERASE": "GST", "SUPEROXIDEDISMUTASE": "SOD",
    "HEATSHOCKTRANSCRIPTIONFACTOR": "HSF", "HEATSHOCKFACTOR": "HSF", "WRKYTRANSCRIPTIONFACTOR": "WRKY",
    "LATEEMBRYOGENESISABUNDANT": "LEA", "UDPGLYCOSYLTRANSFERASE": "UGT", "UDPGLUCOSYLTRANSFERASE": "UGT",
    "POLYGALACTURONASE": "PG", "DIRIGENT": "DIR", "PECTINMETHYLESTERASE": "PME", "CATALASE": "CAT",
}
# canonical family -> parent group (only where the family is a named subfamily or subgroup of a larger family)
GROUP = {
    "R2R3-MYB": "MYB", "MYB-related": "MYB", "MYB-CC": "MYB", "R2R3MYB": "MYB", "1R-MYB": "MYB", "3R-MYB": "MYB",
    "ERF": "AP2/ERF", "DREB": "AP2/ERF", "CBF": "AP2/ERF", "CBF/DREB1": "AP2/ERF", "AP2": "AP2/ERF", "RAV": "AP2/ERF",
    "HD-Zip I": "HD-Zip", "HD-Zip II": "HD-Zip", "HD-Zip III": "HD-Zip", "HD-Zip IV": "HD-Zip",
    "MIKC-type MADS-box": "MADS-box", "MIKCc-type MADS-box": "MADS-box", "NF-YA": "NF-Y", "NF-YB": "NF-Y", "NF-YC": "NF-Y",
    "JAZ": "TIFY", "ABCB": "ABC", "ABCG": "ABC", "ABCC": "ABC", "ABCA": "ABC", "ABCI": "ABC",
    "GA2ox": "GA oxidase", "GA20ox": "GA oxidase", "GA3ox": "GA oxidase", "GAox": "GA oxidase",
    "PHT1": "PHT", "NRT2": "NRT", "Q-type C2H2 ZFP": "C2H2-ZFP", "PP2C-A": "PP2C", "TIR1/AFB": "TIR1/AFB",
    "LRR-RLK": "RLK", "LecRLK": "RLK", "CrRLK1L": "RLK", "LysM-RLK": "RLK", "WAK/WAKL": "RLK",
}
GROUP_KEYS = {}   # filled below: key form of every GROUP entry


def key(s):
    s = re.sub(r"(?<=[A-Z0-9])s$", "", str(s).strip())            # WRKYs -> WRKY; keeps GRAS, TPS
    k = re.sub(r"[\s\-_/().,']", "", s.upper())
    k = re.sub(r"(GENES?|GENEFAMILY|FAMILY|SUPERFAMILY|PROTEINS?|TRANSCRIPTIONFACTORS?)$", "", k) or k
    return k


for g_name, parent in GROUP.items():
    GROUP_KEYS[key(g_name)] = (g_name, parent)


def family_norm(name, family_class="", display=None):
    """Canonical family and its parent group for one family name from an abstract."""
    k = key(name)
    if k == "ARF":
        fam = "ARF" if family_class == "TF" else "Arf GTPase"
        return fam, fam
    if k in GROUP_KEYS:
        return GROUP_KEYS[k]
    fam = SYN.get(k) or (display or {}).get(k) or name.strip()
    return fam, GROUP.get(fam, fam)
