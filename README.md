# Cohort Control Tower

**Track 1 · Cohort Firewall** · OncoSys Cancer AI Lab Healthcare AI Hackathon, Nov 2026
Built on the organizers' challenge kit: [disleveK/oncosys-hackathon-challenge](https://github.com/disleveK/oncosys-hackathon-challenge)

> **Research use only.** This tool checks data releases. It does not diagnose, recommend treatment, score risk or make any clinical claim.

## User and workflow
**User:** a data engineer or researcher who receives a new release of the TCGA-BRCA data.

**The one workflow:** *"A new release (B) arrived. Before anyone reruns results on it, what changed compared with release A, and can analyses A1 and A2 still be trusted?"*

1. Point the tool at release A and release B.
2. It finds every sample that was removed, copied, renamed, swapped or edited, even when the sample count is unchanged.
3. It says which analysis each change affects **and which it doesn't**, reruns A1 and A2, and gives a verdict.
4. It writes `findings.json` in the shared format (v0.2), ready for a Repair Review Desk.

| Analysis | Question | Reads |
|---|---|---|
| **A1** | Do ER-positive tumors have more ER-alpha protein than ER-negative ones? | clinical ER status + protein `ERALPHA` |
| **A2** | What share of primary tumors is each PAM50 subtype? | clinical `PAM50Call_RNAseq` |

## Team
- Shravya Ushake
- Nishant Naravarajula
- Sai Guna Vardhan Vudatala

## Setup
Needs Python 3.10+.
```bash
git clone https://github.com/nish-naravarajula/cohort-control-tower.git
cd cohort-control-tower
pip install -r requirements.txt          # or just: pip install pandas jsonschema cryptography
python data/tcga-brca/fetch.py           # downloads release A (~70 MB) and checks checksums
```

## Run it
From the repo root:
```bash
python control_tower/run.py <release_A_folder> <release_B_folder> findings.json
python format/validate.py findings.json
```

**Try it on the practice envelope (same steps as judging):**
```bash
python envelopes/unseal.py practice orange meadow pilot lantern
python control_tower/run.py data/tcga-brca/raw envelopes/opened/practice/release findings.json
python format/validate.py findings.json
```
On Windows, if `unseal.py` fails with `UnicodeEncodeError`, run `$env:PYTHONIOENCODING="utf-8"` first.

**Expected output (about 2 seconds):**
```
Verdict: Conditional   (3 findings)
F1 [medium] 5 sample(s) present in release A are missing from protein_rppa.tsv.gz in release B.
F2 [high]   5 new sample ID(s) in protein_rppa.tsv.gz carry profiles identical to existing samples.
F3 [high]   1 pair(s) of samples traded their profiles: TCGA-A8-A09K-01 <-> TCGA-B6-A0WX-01.
A1: changed
A2: unchanged
```

## How it works
```
release A ─┐
           ├─► detect.py ─► impact.py ─► report.py ─► findings.json
release B ─┘  what changed  who is hurt   write it up
```
| File | What it does |
|---|---|
| `control_tower/detect.py` | Fingerprints (SHA-256) each file, then each **sample's** values; pairs samples by ID and classifies each difference |
| `control_tower/impact.py` | Reruns the organizers' A1/A2 and decides which changes reach which analysis |
| `control_tower/report.py` | Groups changes into findings with evidence, repairs, before/after results and a verdict |
| `control_tower/run.py` | One command that runs it all |

**Detecting changes.** Each sample's row of values gets one fingerprint. Looking up B's fingerprint in A tells us where the values came from:

| What we see | Classified as |
|---|---|
| ID in A, not in B | **removed** |
| New ID, fingerprint matches a sample still in B | **duplicate** |
| New ID, fingerprint matches a sample now gone | **renamed** |
| New ID, fingerprint new | **added** |
| Same ID, now holds a sample's values and vice versa | **swap** |
| Same ID, holds another sample's values (no trade) | **moved** |
| Same ID, fingerprint new | **edited** (and which column changed) |
| Same ID on two rows | **repeated ID** |
| Columns appear or disappear | **columns added/removed** |

Reordered rows or columns are **not** reported: the same data in a different order gives the same fingerprint.

**Deciding impact.** A change affects an analysis only if all three hold:
1. the analysis **reads that file**,
2. the sample is **in the analysis's own sample list** (in A or B),
3. for an edit, the changed **column is one the analysis reads** (A1: ER status, ERALPHA; A2: PAM50).

**Verdict:** `Ready` (nothing reaches A1/A2) · `Conditional` (problems found, each with a proposed repair) · `Insufficient Evidence` (an analysis couldn't be rerun, e.g. a file is missing).

## One result you can verify
On the practice release, protein profiles of `TCGA-A8-A09K-01` (ER-positive) and `TCGA-B6-A0WX-01` (ER-negative) are swapped: ERALPHA `3.43` ↔ `-1.07`. A1 moves from a difference of **3.051 → 3.019**, and A2 is unchanged. Check it yourself:
```bash
python data/tcga-brca/reference_analyses.py data/tcga-brca/raw
python data/tcga-brca/reference_analyses.py envelopes/opened/practice/release
```
Our output matches the organizers' reference (`format/examples/control-tower-practice-B.json`): same samples, kinds, severities, affected analyses and verdict.

**Tested on 17 releases:** 8 made with `make_release.py` (ER and PAM50 relabels, swaps, duplicates, irrelevant edits, changes in mutation and copy-number files) found **44/44 changed samples with 0 false alarms**. 9 hand-made edge cases (reordered columns, missing file, renamed ID, repeated rows, dropped column, 3-way rotation) were all handled. Every output passes `validate.py`.

## Limitations
- **Exact matches only.** A near-copy with one value changed is reported as an *edit*, not a duplicate.
- A swap whose values were **also altered** shows as two edits.
- It shows **that** releases differ, not **which one is correct**, so every finding is marked *needs human review*.
- Formatting changes (e.g. rounding `3.4327` → `3.43`) count as edits.
- Missing RNA in judge releases is expected and reported as a note. Neither A1 nor A2 uses RNA.

## Out of scope
No diagnosis, treatment advice, risk scores or clinical claims. No re-identification or linking to outside records. No API keys, and no raw TCGA files are committed (they're fetched).

## Demo video
<link>

## Third-party code and data
- [pandas](https://pandas.pydata.org) (BSD-3), [jsonschema](https://github.com/python-jsonschema/jsonschema) (MIT), [cryptography](https://github.com/pyca/cryptography) (Apache-2.0 / BSD)
- Organizers' kit: `reference_analyses.py`, `fetch.py`, `make_release.py`, `envelopes/unseal.py`, `format/` ([challenge repo](https://github.com/disleveK/oncosys-hackathon-challenge))
- Data: TCGA-BRCA via UCSC Xena. *"The results shown here are in whole or part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga."* Full credits in [DATA_LICENSES.md](DATA_LICENSES.md).