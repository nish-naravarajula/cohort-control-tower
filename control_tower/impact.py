import os, sys, pathlib
import pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "data" / "tcga-brca"))
import reference_analyses as ra

# which files each analysis reads, and which columns in them (protein columns are protein names)
USES = {
    "A1": {"clinical.tsv": ["breast_carcinoma_estrogen_receptor_status"], "protein_rppa.tsv.gz": ["ERALPHA"]},
    "A2": {"clinical.tsv": ["PAM50Call_RNAseq"]},
}


def run_analyses(release):
    """Run A1/A2 with the organizers' code. An analysis whose input file is missing returns None."""
    r = pathlib.Path(release)
    out = {"A1": None, "A2": None}
    if not (r / "clinical.tsv").exists():
        return out
    clinical = pd.read_csv(r / "clinical.tsv", sep="\t", low_memory=False)
    out["A2"] = ra.a2(clinical)
    if (r / "protein_rppa.tsv.gz").exists():
        protein = pd.read_csv(r / "protein_rppa.tsv.gz", sep="\t", index_col=0)
        out["A1"] = ra.a1(clinical, protein)
    return out


def events(report):
    out = []
    for file, res in report.items():
        if isinstance(res, str):
            if res.startswith("only"):
                out.append({"kind": "file_" + res.replace(" ", "_"), "file": file, "samples": [], "columns": None})
            continue
        for s in res["removed"]:
            out.append({"kind": "removed", "file": file, "samples": [s], "columns": None})
        for new, src in res["duplicates"]:
            out.append({"kind": "duplicate", "file": file, "samples": [new], "source": src, "columns": None})
        for new, old in res["renamed"]:
            out.append({"kind": "renamed", "file": file, "samples": [new, old], "source": old, "columns": None})
        for s in res["added_new"]:
            out.append({"kind": "added", "file": file, "samples": [s], "columns": None})
        for a, b in res["swaps"]:
            out.append({"kind": "swap", "file": file, "samples": [a, b], "columns": None})
        for s, src in res["moved"]:
            out.append({"kind": "moved", "file": file, "samples": [s], "source": src, "columns": None})
        for s, cols in res["edited"].items():
            out.append({"kind": "edited", "file": file, "samples": [s], "columns": cols})
        for s in res["repeated_ids"]:
            out.append({"kind": "repeated_id", "file": file, "samples": [s], "columns": None})
        if res["columns_removed"]:
            out.append({"kind": "columns_removed", "file": file, "samples": [], "columns": res["columns_removed"]})
        if res["columns_added"]:
            out.append({"kind": "columns_added", "file": file, "samples": [], "columns": res["columns_added"]})
    return out


def hits(name, ev, res_A, res_B):
    reads = USES[name]
    if ev["file"] not in reads:
        return False                                    # (a) analysis doesn't read this file
    cols = reads[ev["file"]]
    if not ev["samples"]:                               # whole-file or whole-column change
        return ev["columns"] is None or any(c in cols for c in ev["columns"])
    used = set()
    for res in (res_A[name], res_B[name]):
        if res:
            used |= set(res["sample_ids"])
    if not any(s in used for s in ev["samples"]):
        return False                                    # (b) sample not in the analysis
    if ev["columns"] is None:
        return True                                     # whole sample moved/copied/removed
    return any(c in cols for c in ev["columns"])        # (c) edited column is one it reads


def add_impact(evs, res_A, res_B):
    for ev in evs:
        ev["affects"] = [n for n in USES if hits(n, ev, res_A, res_B)]
        ev["not_affected"] = [n for n in USES if n not in ev["affects"]]
    return evs
