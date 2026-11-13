import gzip, hashlib
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
    if name in TALL:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        return df.set_index(TALL[name])
    return pd.read_csv(path, sep="\t", index_col=0, dtype=str, keep_default_na=False).T


def hash_samples(df):
    return df.apply(lambda row: hashlib.sha256("|".join(row.values).encode()).hexdigest(), axis=1)


def compare(hash_A, hash_B):
    owner_in_A = {h: s for s, h in hash_A.items()}
    names_A, names_B = set(hash_A.index), set(hash_B.index)
    result = {"removed": sorted(names_A - names_B), "duplicates": [], "added_new": [],
              "swaps": [], "moved": [], "edited": []}
    for s in sorted(names_B - names_A):
        src = owner_in_A.get(hash_B[s])
        if src:
            result["duplicates"].append((s, src))
        else:
            result["added_new"].append(s)
    changed = sorted(s for s in names_A & names_B if hash_A[s] != hash_B[s])
    for s in changed:
        src = owner_in_A.get(hash_B[s])
        if src is None:
            result["edited"].append(s)
        elif owner_in_A.get(hash_B.get(src)) == s:
            if s < src:
                result["swaps"].append((s, src))
        else:
            result["moved"].append((s, src))
    return result


def compare_releases(dir_A, dir_B):
    report = {}
    for name in FILES:
        if file_hash(f"{dir_A}/{name}") == file_hash(f"{dir_B}/{name}"):
            report[name] = "unchanged"
            continue
        a = hash_samples(load_samples(f"{dir_A}/{name}", name))
        b = hash_samples(load_samples(f"{dir_B}/{name}", name))
        report[name] = compare(a, b)
    return report