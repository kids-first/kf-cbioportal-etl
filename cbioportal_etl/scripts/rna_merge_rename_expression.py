#!/usr/bin/env python3
"""Merge rsem gene outputs for cBio usage.

Converts to wide-format
Replace BS ID with cBio ID
Repeat gene names merge to take highest mean expression
Final result also z scored across cohort as log2(FPKM + 1) or as log2(TPM + 1)
"""

import argparse
import concurrent.futures
import io
import json
import multiprocessing as mp
import os
import sys
import tarfile
import time
import warnings

import numpy as np
import pandas as pd

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths

# Per-worker globals (set once per process by _init_worker)
_REF_INDEX = None
_EXPR = None
_RSEM_DIR = None
_W = {}  # set in the parent before the pool forks


def _init_worker(ref_index, expr_type, rsem_dir):
    global _REF_INDEX, _EXPR, _RSEM_DIR
    _REF_INDEX, _EXPR, _RSEM_DIR = ref_index, expr_type, rsem_dir


def load_rsem_file(rsem_file: str):
    """Read one RSEM file; return the expression column as float32 aligned to the reference gene index."""
    try:
        s = pd.read_csv(
            os.path.join(_RSEM_DIR, rsem_file),
            sep="\t",
            index_col="gene_id",
            usecols=["gene_id", _EXPR],
            dtype={_EXPR: np.float32},
        )[_EXPR]
        if not s.index.equals(_REF_INDEX):
            s = s[~s.index.duplicated()].reindex(_REF_INDEX)
        return s.to_numpy(dtype=np.float32)
    except Exception as e:
        print(f"{e}\nFailed to process {rsem_file}", file=sys.stderr)
        return None


def gene_symbols(gene_ids) -> np.ndarray:
    """'ENSG000_TP53' -> 'TP53'"""
    return pd.Series(np.asarray(gene_ids)).str.split("_", n=1).str[1].to_numpy()


def dedupe_by_max(symbols: np.ndarray, mean_expr: np.ndarray) -> np.ndarray:
    """Row positions to keep: one per symbol (highest mean), ordered by descending mean."""
    tmp = pd.DataFrame({"sym": symbols, "m": mean_expr})
    tmp = tmp.sort_values("m", ascending=False, kind="stable").drop_duplicates("sym", keep="first")
    return tmp.index.to_numpy()


def zscore_rows(a: np.ndarray, block: int = 4096) -> np.ndarray:
    """Row-wise z-score (ddof=0, NaN-aware). NaN / zero-variance results become 0."""
    out = np.empty(a.shape, dtype=np.float32)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for s in range(0, a.shape[0], block):
            x = a[s:s + block]
            if np.isnan(x).any():
                mu = np.nanmean(x, axis=1, keepdims=True, dtype=np.float64)
                sd = np.nanstd(x, axis=1, keepdims=True, dtype=np.float64)
            else:
                mu = x.mean(axis=1, keepdims=True, dtype=np.float64)
                sd = x.std(axis=1, keepdims=True, dtype=np.float64)
            sd[sd == 0] = np.inf  # zero variance -> 0
            z = out[s:s + block]
            np.subtract(x, mu, out=z, casting="same_kind")
            z /= sd
            z[np.isnan(z)] = 0
    return out


def load_healthy_stats(archive_path: str, expr_type: str, needed_keys: set) -> dict:
    """Return {key: (mean Series, std Series)} indexed by Hugo symbol, for needed references only."""
    stats_by_key = {}
    if not needed_keys:
        return stats_by_key
    with tarfile.open(archive_path, "r:gz") as tar:
        for member in tar:  # iterate lazily instead of getmembers()
            if not member.isfile():
                continue
            name = os.path.basename(member.name)
            if not name.endswith(f"_{expr_type}.tsv"):
                continue
            key = name.replace(f"_log_{expr_type}.tsv", "") + f"_{expr_type}"
            if key not in needed_keys:
                continue
            file_obj = tar.extractfile(member)
            if not file_obj:
                continue
            print(f"{time.time()} Loading healthy reference {key}", file=sys.stderr)
            df = pd.read_csv(io.TextIOWrapper(file_obj), sep="\t", index_col=0)
            symbols = gene_symbols(df.index)
            keep = dedupe_by_max(symbols, df.mean(axis=1).to_numpy())
            sub = df.iloc[keep]
            idx = pd.Index(symbols[keep])
            mu, sd = sub.mean(axis=1), sub.std(axis=1)
            mu.index, sd.index = idx, idx
            stats_by_key[key] = (mu, sd)
            del df, sub
            if len(stats_by_key) == len(needed_keys):
                break
    return stats_by_key


def _fmt_chunk(span):
    s, e = span
    rows, cols, mat = _W["rows"], _W["cols"], _W["mat"]
    if rows is None:
        block = mat[s:e][:, cols]
        g = _W["genes"][s:e]
    else:
        r = rows[s:e]
        block = mat[np.ix_(r, cols)]
        g = _W["genes"][r]
    return pd.DataFrame(block, index=g).to_csv(
        sep="\t", header=False, float_format=_W["ff"]
    )


def write_projects(mat, genes, project_samples, col_pos, path_fmt,
                   rows=None, float_format=None, chunk=2000, workers=None):
    n = mat.shape[0] if rows is None else len(rows)
    spans = [(s, min(s + chunk, n)) for s in range(0, n, chunk)]
    ctx = mp.get_context("fork")
    for project, names in project_samples.items():
        _W.update(mat=mat, genes=genes, rows=rows,
                  cols=np.array([col_pos[x] for x in names]), ff=float_format)
        with ctx.Pool(workers) as pool, open(path_fmt.format(project=project), "w") as fh:
            fh.write("\t".join(["Hugo_Symbol", *names]) + "\n")
            for text in pool.imap(_fmt_chunk, spans):
                fh.write(text)


def match_type_for(library_type, default_match_type):
    if library_type is None:
        key = None
    else:
        key = str(library_type).strip().lower()
    if key == "rrna depletion":
        return "totalRNA"
    if key in {"poly-t enrichment", "hybrid selection"}:
        return "polyA"
    return default_match_type if default_match_type != "none" else None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge rsem gene files using manifest file info.")
    parser.add_argument("-t", "--table", dest="table",
                        help="Table with cbio project, kf bs ids, cbio IDs, and file names")
    parser.add_argument("-r", "--rsem-dir", dest="rsem_dir", help="rsem file directory")
    parser.add_argument("-et", "--expression-type", dest="expression_type",
                        choices=["TPM", "FPKM"], default="TPM",
                        help="Which expression value to use: TPM or FPKM. Default is TPM.")
    parser.add_argument("-sc", "--study-config", dest="study_config", help="cbio study config file.")
    parser.add_argument(
        "-dmt", "--default-match-type", dest="default_match_type",
        choices=["polyA", "totalRNA", "none"], default="none",
        help="Default match type for samples with unknown RNA library type for z-score calculations. "
        "Use 'polyA' or 'totalRNA' to override fallback to intra-cohort z-score.",
    )
    args = parser.parse_args()

    expr = args.expression_type
    rsem_dir = args.rsem_dir.rstrip("/")
    out_dir = "merged_rsem/"
    os.makedirs(out_dir, exist_ok=True)

    all_file_meta = pd.read_csv(args.table, sep="\t", dtype={"cbio_sample_name": str})
    rna_subset = all_file_meta[all_file_meta["etl_file_type"] == "rsem"].copy()
    pairs = rna_subset[["file_name", "cbio_sample_name"]].drop_duplicates(subset="cbio_sample_name")
    files = pairs["file_name"].tolist()
    samples = pairs["cbio_sample_name"].tolist()

    # ---- Load RSEM files into one preallocated float32 matrix (genes x samples) ----
    print(f"{time.time()} Reading RSEM files...", file=sys.stderr)
    ref_index = pd.Index(
        pd.read_csv(os.path.join(rsem_dir, files[0]), sep="\t", usecols=["gene_id"])["gene_id"]
    )
    master = np.empty((len(ref_index), len(files)), dtype=np.float32)
    loaded = []
    with concurrent.futures.ProcessPoolExecutor(
        initializer=_init_worker, initargs=(ref_index, expr, rsem_dir)
    ) as executor:
        for i, vals in enumerate(executor.map(load_rsem_file, files, chunksize=8), 1):
            if vals is not None:
                master[:, len(loaded)] = vals
                loaded.append(samples[i - 1])
            if i % 50 == 0:
                print(f"{time.time()} Loaded {i} files", file=sys.stderr)
    master = master[:, : len(loaded)]
    col_pos = {s: i for i, s in enumerate(loaded)}

    # ---- Collapse duplicate gene symbols (keep highest mean expression) ----
    symbols = gene_symbols(ref_index)
    row_mean = (np.nanmean if np.isnan(master).any() else np.mean)(master, axis=1, dtype=np.float64)
    keep = dedupe_by_max(symbols, row_mean)
    master = master[keep]  # contiguous copy of the kept rows only
    genes = pd.Index(symbols[keep], name="Hugo_Symbol")
    del row_mean, keep, symbols

    project_samples = {
        proj: [s for s in dict.fromkeys(grp["cbio_sample_name"]) if s in col_pos]
        for proj, grp in rna_subset.groupby("cbio_project", sort=False)
    }

    # ---- Raw expression output ----
    print(f"{time.time()} Outputting {expr} expression results", file=sys.stderr)
    write_projects(master, genes, project_samples, col_pos,
                   out_dir + "{project}.rsem_merged." + expr + ".txt")

    # ---- log2(x + 1) in place ----
    master += 1
    np.log2(master, out=master)
    log_master = master

    print(f"{time.time()} Calculating z-scores...", file=sys.stderr)
    if "etl_experiment_strategy" in rna_subset.columns:
        # sample lists per library type (None == missing)
        groups = {}
        for s, lib in zip(rna_subset["cbio_sample_name"], rna_subset["etl_experiment_strategy"]):
            groups.setdefault(None if pd.isna(lib) else lib, {})[s] = None

        # column positions per library type (only samples that loaded)
        group_cols = {}
        for lib, sample_dict in groups.items():
            names = [s for s in sample_dict if s in col_pos]
            if names:
                group_cols[lib] = np.array([col_pos[s] for s in names])

        # ---- Intra-cohort z-score per library type: always written ----
        intra = np.zeros(log_master.shape, dtype=np.float32)
        for lib, cols in group_cols.items():
            print(f"{time.time()} Calculating intra-cohort z-score for {lib}", file=sys.stderr)
            intra[:, cols] = zscore_rows(log_master[:, cols])
        write_projects(
            intra, genes, project_samples, col_pos,
            out_dir + "{project}.rsem_merged_tumor_only_zscore_" + expr + ".txt",
            float_format="%.4f",
        )
        del intra  # free before the vs-healthy matrix is allocated

        # ---- Tumor vs healthy z-score ----
        NONE_TYPES = {"rnaseq", "rna-seq", "rna seq", "none", "nan", ""}
        only_lib = next(iter(groups)) if len(groups) == 1 else "<multiple>"
        all_none = len(groups) == 1 and (
            only_lib is None or str(only_lib).strip().lower() in NONE_TYPES
        )
        skip_vs_healthy = all_none and str(args.default_match_type).strip().lower() in {"none", "", "nan"}

        if skip_vs_healthy:
            print(f"{time.time()} Skipping output of vs-healthy z-score because values are equivalent "
                  "to intra-cohort z-score", file=sys.stderr)
        else:
            # Load only the healthy references we actually need, reduced to mean/std
            TOOL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            with open(args.study_config) as f:
                config_data = resolve_config_paths(json.load(f), TOOL_DIR)
            needed = {
                f"healthy_{mt}_{expr}"
                for mt in (match_type_for(lib, args.default_match_type) for lib in group_cols)
                if mt
            }
            healthy = load_healthy_stats(config_data["rsem_ref"]["archive"], expr, needed)

            out = np.zeros(log_master.shape, dtype=np.float32)
            covered = np.zeros(len(genes), dtype=bool)

            for lib, cols in group_cols.items():
                print(f"{time.time()} Calculating tumor vs healthy z-score for {lib}", file=sys.stderr)
                mt = match_type_for(lib, args.default_match_type)
                print(f"{time.time()}  ➤ Match type being used: {mt}", file=sys.stderr)
                ref = healthy.get(f"healthy_{mt}_{expr}") if mt else None

                if ref is not None:
                    mu_s, sd_s = ref
                    pos = np.flatnonzero(genes.isin(mu_s.index))
                    x = log_master[np.ix_(pos, cols)]  # copy, safe to modify in place
                    mu = mu_s.reindex(genes[pos]).to_numpy(np.float32)[:, None]
                    sd = sd_s.reindex(genes[pos]).to_numpy(np.float32)[:, None]
                    with np.errstate(all="ignore"):
                        x -= mu
                        x /= sd
                    x[np.isnan(x)] = 0
                    out[np.ix_(pos, cols)] = x
                    covered[pos] = True
                    del x
                else:
                    print(f"{time.time()} No healthy reference available for {lib} library type, "
                          "using intra-cohort z-score instead", file=sys.stderr)
                    # recomputed rather than kept in memory from the step above
                    out[:, cols] = zscore_rows(log_master[:, cols])
                    covered[:] = True

            rows = np.flatnonzero(covered)
            write_projects(
                out, genes, project_samples, col_pos,
                out_dir + "{project}.rsem_merged_vs_healthy_zscore_" + expr + ".txt",
                rows=rows, float_format="%.4f",
            )

    else:
        print(f"{time.time()} No etl_experiment_strategy column found, using intra-cohort z-score", file=sys.stderr)
        z = zscore_rows(log_master)
        write_projects(z, genes, project_samples, col_pos,
                       out_dir + "{project}.rsem_merged_tumor_only_zscore_" + expr + ".txt",
                       float_format="%.4f")
