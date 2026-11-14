import gzip, hashlib, os
import pandas as pd

FILES = ["clinical.tsv", "survival.tsv", "mutation_gene_level.tsv.gz",
         "copy_number_gistic.tsv.gz", "rna_hiseq.tsv.gz", "protein_rppa.tsv.gz"]
TALL = {"clinical.tsv": "sampleID", "survival.tsv": "sample"}


def file_hash(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as f:
        data = f.read().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def load_samples(path, name):
    """One row per sample, one column per field / gene / protein."""
    if name in TALL:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).set_index(TALL[name])
    return pd.read_csv(path, sep="\t", index_col=0, dtype=str, keep_default_na=False).T


def hash_samples(df):
    rows = df.to_numpy()                       # plain array: much faster than df.apply
    return pd.Series([hashlib.sha256("|".join(r).encode()).hexdigest() for r in rows], index=df.index)


def owners(hashes):
    """hash -> list of samples with that profile (several samples can legitimately share one)."""
    out = {}
    for s, h in hashes.items():
        out.setdefault(h, []).append(s)
    return out


def patient(s):
    return s[:12]


def compare(hash_A, hash_B):
    own_A = owners(hash_A)
    names_A, names_B = set(hash_A.index), set(hash_B.index)
    result = {"removed": [], "duplicates": [], "renamed": [], "added_new": [],
              "swaps": [], "moved": [], "edited": []}

    renamed_from = set()
    for s in sorted(names_B - names_A):
        same = own_A.get(hash_B[s], [])
        gone = [o for o in same if o not in names_B]
        src = [o for o in same if patient(o) != patient(s)]   # a same-patient sibling may share data
        if gone:
            result["renamed"].append((s, gone[0]))       # same profile, old ID disappeared
            renamed_from.add(gone[0])
        elif src:
            result["duplicates"].append((s, src[0]))     # copy of a sample that is still there
        else:
            result["added_new"].append(s)                # new profile (or a same-patient sibling)

    result["removed"] = sorted((names_A - names_B) - renamed_from)

    seen = set()
    for s in sorted(n for n in names_A & names_B if hash_A[n] != hash_B[n]):
        src = [o for o in own_A.get(hash_B[s], []) if o != s]
        if not src:
            result["edited"].append(s)
            continue
        partner = next((o for o in src if o in names_B and hash_B[o] == hash_A[s]), None)
        if partner:
            if partner not in seen:
                result["swaps"].append((s, partner))
                seen.update({s, partner})
        else:
            result["moved"].append((s, src[0]))
    return result


def changed_columns(df_A, df_B, sample):
    a, b = df_A.loc[sample], df_B.loc[sample]
    return {c: (a[c], b[c]) for c in df_A.columns if a[c] != b[c]}


def compare_releases(dir_A, dir_B):
    report = {}
    for name in FILES:
        in_A = os.path.exists(f"{dir_A}/{name}")
        in_B = os.path.exists(f"{dir_B}/{name}")
        if not in_A and not in_B:
            continue
        if not (in_A and in_B):
            report[name] = "only in A" if in_A else "only in B"
            continue
        if file_hash(f"{dir_A}/{name}") == file_hash(f"{dir_B}/{name}"):
            report[name] = "unchanged"
            continue

        df_A = load_samples(f"{dir_A}/{name}", name)
        df_B = load_samples(f"{dir_B}/{name}", name)

        # the same sample ID appearing more than once in a file
        repeated = sorted(set(df_B.index[df_B.index.duplicated()]))
        if df_A.index.has_duplicates:
            df_A = df_A[~df_A.index.duplicated()]
        if repeated:
            df_B = df_B[~df_B.index.duplicated()]

        # columns (clinical fields, genes, proteins) added, removed or reordered:
        # compare only the shared ones, in release A's order (skipped when nothing differs, to stay fast)
        cols_A, cols_B = set(df_A.columns), set(df_B.columns)
        cols_removed = [c for c in df_A.columns if c not in cols_B]
        cols_added = [c for c in df_B.columns if c not in cols_A]
        if list(df_A.columns) != list(df_B.columns):
            shared = [c for c in df_A.columns if c in cols_B]
            df_A, df_B = df_A[shared], df_B[shared]

        res = compare(hash_samples(df_A), hash_samples(df_B))
        res["edited"] = {s: changed_columns(df_A, df_B, s) for s in res["edited"]}
        res["repeated_ids"] = repeated
        res["columns_removed"] = cols_removed
        res["columns_added"] = cols_added
        if not any(res.values()):
            res = "unchanged (layout only)"            # reordered rows/columns, same data
        report[name] = res
    return report
