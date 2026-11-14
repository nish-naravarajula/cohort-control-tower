import os
from detect import FILES, file_hash

# for each change type: finding kind, severity, repair action, limitation, repair justification
RULES = {
    "swap":      ("sample_swap", "high", "restore_from_previous_release",
                  "A swap is inferred from profiles trading places; the release does not say why.",
                  "Restore these profiles as they were in the previous release until the swap is explained."),
    "moved":     ("sample_swap", "high", "restore_from_previous_release",
                  "The sample now carries another sample's profile; the release does not say why.",
                  "Restore these profiles from the previous release until the change is explained."),
    "duplicate": ("duplicate", "high", "exclude",
                  "Identical profiles strongly suggest copies, not new measurements.",
                  "Exclude the duplicated profiles from analyses until the provider confirms them."),
    "removed":   ("membership_change", "medium", "flag_for_review",
                  "The tool sees the samples are gone but not why they were removed.",
                  "Ask the data provider why these samples lost their data."),
    "added":     ("membership_change", "medium", "flag_for_review",
                  "New samples with new values cannot be checked against the previous release.",
                  "Confirm with the data provider that these are genuine new samples."),
    "edited":    ("label_conflict", "high", "flag_for_review",
                  "The tool sees that the value changed but cannot tell which release has the correct label.",
                  "A person must confirm which label is correct before the analysis is used."),
    "renamed":   ("membership_change", "high", "relabel",
                  "The profile is identical to a sample whose ID disappeared; the tool assumes the ID was changed.",
                  "Confirm the ID change with the data provider, or relabel the sample back to its previous ID."),
    "repeated_id": ("duplicate", "high", "exclude",
                  "The same sample ID appears more than once in the file; the tool compared only its first row.",
                  "Keep one row per sample ID after the provider confirms which row is correct."),
    "columns_removed": ("missing_data", "medium", "flag_for_review",
                  "Whole columns disappeared; the tool compared only the columns present in both releases.",
                  "Ask the data provider why these columns were removed."),
    "columns_added": ("other", "low", "flag_for_review",
                  "New columns have no previous values to compare against.",
                  "Confirm the new columns with the data provider."),
    "file_only_in_A": ("missing_data", "medium", "flag_for_review",
                  "The file could not be compared because it exists only in release A.",
                  "Ask the data provider why this file is missing."),
    "file_only_in_B": ("membership_change", "low", "flag_for_review",
                  "The file could not be compared because it exists only in release B.",
                  "Confirm the new file with the data provider."),
}
EXPECTED_MISSING = {"rna_hiseq.tsv.gz"}   # judge releases ship without RNA on purpose


def group_key(ev):
    # one finding per (change type, file) - and for edits, per set of changed columns
    cols = tuple(sorted(ev["columns"])) if ev["kind"] == "edited" else ()
    return (ev["kind"], ev["file"], cols)


def describe(kind, file, group, cols):
    n = len(group)
    if kind == "swap":
        pairs = ", ".join(f"{e['samples'][0]} <-> {e['samples'][1]}" for e in group)
        return f"{n} pair(s) of samples traded their profiles in {file}: {pairs}."
    if kind == "duplicate": return f"{n} new sample ID(s) in {file} carry profiles identical to existing samples."
    if kind == "removed":   return f"{n} sample(s) present in release A are missing from {file} in release B."
    if kind == "added":     return f"{n} new sample(s) in {file} have values not seen in release A."
    if kind == "moved":     return f"{n} sample(s) in {file} now carry another sample's profile."
    if kind == "edited":    return f"{n} sample(s) changed {', '.join(cols)} in {file}."
    if kind == "renamed":   return f"{n} sample(s) in {file} reappear under a new ID with an identical profile."
    if kind == "repeated_id": return f"{n} sample ID(s) appear more than once in {file}."
    if kind in ("columns_removed", "columns_added"):
        c = group[0]["columns"]
        what = "removed from" if kind == "columns_removed" else "added to"
        return f"{len(c)} column(s) {what} {file}: {', '.join(c[:10])}{' ...' if len(c) > 10 else ''}."
    return f"{file} is {kind.replace('file_', '').replace('_', ' ')}."


def evidence(ev, name_A, name_B):
    k, f = ev["kind"], ev["file"]
    if k == "removed":
        return [{"kind": "sample", "id": ev["samples"][0], "release": name_A, "file": f,
                 "note": f"present in {name_A}, absent in {name_B}"}]
    if k in ("duplicate", "moved"):
        return [{"kind": "sample", "id": ev["samples"][0], "release": name_B, "file": f,
                 "note": f"identical to {ev['source']} in {name_A}"}]
    if k == "edited":
        note = "; ".join(f"{c}: '{a}' -> '{b}'" for c, (a, b) in ev["columns"].items())
        return [{"kind": "sample", "id": ev["samples"][0], "release": name_B, "file": f, "note": note}]
    if k == "renamed":
        new, old = ev["samples"]
        return [{"kind": "sample", "id": new, "release": name_B, "file": f,
                 "note": f"identical to {old} in {name_A}; {old} is absent from {name_B}"},
                {"kind": "sample", "id": old, "release": name_A, "file": f, "note": "previous ID"}]
    if k == "repeated_id":
        return [{"kind": "sample", "id": ev["samples"][0], "release": name_B, "file": f,
                 "note": "this ID has more than one row"}]
    if k in ("columns_removed", "columns_added"):
        return [{"kind": "file", "id": f, "release": name_B, "note": "columns: " + ", ".join(ev["columns"])}]
    if k.startswith("file_"):
        return [{"kind": "file", "id": f}]
    return [{"kind": "sample", "id": s, "release": name_B, "file": f} for s in ev["samples"]]


def summarize(res):
    if res is None:
        return None
    return {k: v for k, v in res.items() if k not in ("sample_ids", "analysis", "inputs")}


def release_files(path):
    return [{"file": n, "sha256": file_hash(f"{path}/{n}"), "hash_of": "content"}
            for n in FILES if os.path.exists(f"{path}/{n}")]


def build(dir_A, dir_B, evs, res_A, res_B, team="Cohort Control Tower team"):
    name_A = os.path.basename(os.path.normpath(dir_A))
    name_B = os.path.basename(os.path.normpath(dir_B))

    # RNA missing from a judge release is expected: report it as a note, not a finding
    notes = [f"{e['file']} is not in release B and was not compared (judge releases ship without RNA)."
             for e in evs if e["kind"] == "file_only_in_A" and e["file"] in EXPECTED_MISSING]
    evs = [e for e in evs if not (e["kind"] == "file_only_in_A" and e["file"] in EXPECTED_MISSING)]

    groups = {}
    for ev in evs:
        groups.setdefault(group_key(ev), []).append(ev)

    findings, repairs = [], []
    all_analyses = list(res_A)
    for i, ((kind, file, cols), group) in enumerate(groups.items(), 1):
        fkind, sev, action, limitation, why = RULES[kind]
        affects = sorted({a for e in group for a in e["affects"]})
        if not affects:
            sev = "low"
        refs = [r for e in group for r in evidence(e, name_A, name_B)]
        if kind in ("swap", "moved") and "A1" in affects and file == "protein_rppa.tsv.gz":
            refs.append({"kind": "protein", "id": "ERALPHA"})
        fid = f"F{i}"
        findings.append({"id": fid, "kind": fkind, "summary": describe(kind, file, group, cols),
                         "severity": sev, "confidence": "high", "evidence": refs,
                         "affects_analyses": affects,
                         "not_affected_analyses": [a for a in all_analyses if a not in affects],
                         "limitation": limitation, "needs_human_review": True})
        repairs.append({"id": f"R{i}", "action": action,
                        "targets": [r for r in refs if r["kind"] in ("sample", "file")],
                        "justification": why, "addresses": [fid]})

    analyses = []
    for name in all_analyses:
        a, b = summarize(res_A[name]), summarize(res_B[name])
        if a is None or b is None:            # an input file is missing, so it could not be rerun
            analyses.append({"id": name, "status": "not_rerun", **({"before": a} if a else {})})
            continue
        analyses.append({"id": name, "status": "unchanged" if a == b else "changed",
                         "before": a, "after": b})

    # Ready = nothing touches A1/A2.  Conditional = problems found, each with a proposed repair.
    hit = any(f["affects_analyses"] for f in findings)
    if any(a["status"] == "not_rerun" for a in analyses):
        verdict = "Insufficient Evidence"      # an analysis could not even be rerun on release B
    elif hit:
        verdict = "Conditional"
    else:
        verdict = "Ready"

    doc = {"format_version": "0.2",
           "producer": {"team": team, "product": "Cohort Control Tower", "track": "cohort-firewall"},
           "data": {"dataset": "TCGA-BRCA", "releases": [
               {"name": name_A, "files": release_files(dir_A)},
               {"name": name_B, "files": release_files(dir_B)}]},
           "question": f"Did release {name_B} change anything that affects analyses A1 "
                       f"(ER-alpha protein by ER status) and A2 (RNA subtype mix)?",
           "verdict": verdict, "findings": findings, "analyses": analyses, "repairs": repairs}
    return doc, notes
