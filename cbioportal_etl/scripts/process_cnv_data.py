#!/usr/bin/env python
"""Read in GATK and ControlFreeC CNV data and convert to cBio format.

Gene overlaps are computed in-process with numpy (no bedtools subprocess per sample), and
per-sample results are kept as integer arrays in a genes x samples matrix instead of dicts.
"""

import argparse
import csv
import json
import logging
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths

logger = logging.getLogger(__name__)

PVALUE: float = 0.05
NOT_SEEN: int = int(np.iinfo(np.int32).min)  # marks genes with no CNV overlap in a sample
SEG_HEADER = "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean"
WRITE_CHUNK_ROWS = 2000  # genes formatted per write when producing wide tables


class SampleResult(NamedTuple):
    """Per-sample output of process_sample."""

    raw_cn: np.ndarray  # per-gene CN in GeneIndex order, NOT_SEEN where no CNV overlaps the gene
    ploidy: int
    seg_text: str  # tab separated seg rows, newline terminated


class GeneIndex:
    """Gene intervals from a 4 column BED (chrom, start, end, gene), indexed per chromosome."""

    def __init__(self, bed_path: Path) -> None:
        ref = pd.read_csv(
            bed_path, sep="\t", header=None, usecols=[0, 1, 2, 3], comment="#",
            names=["chrom", "start", "end", "gene"],
            dtype={"chrom": str, "start": np.int64, "end": np.int64, "gene": str},
        )
        codes, uniques = pd.factorize(ref["gene"])
        ref["gene_id"] = codes
        self.genes: np.ndarray = np.asarray(uniques, dtype=object)  # output order: first appearance
        self.n_genes: int = len(self.genes)
        # per chromosome: starts (sorted), ends, running max of ends, gene ids
        self._by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
        for chrom, grp in ref.sort_values("start", kind="stable").groupby("chrom", sort=False):
            ends = grp["end"].to_numpy()
            self._by_chrom[chrom] = (
                grp["start"].to_numpy(), ends, np.maximum.accumulate(ends), grp["gene_id"].to_numpy(),
            )

    def overlaps(self, segments: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Find segment/gene overlaps (BED half-open, at least 1 bp).

        Returns:
            Row positions in segments, gene ids, and overlap sizes in bp, ordered by segment row.

        """
        a_parts: list[np.ndarray] = []
        g_parts: list[np.ndarray] = []
        o_parts: list[np.ndarray] = []
        starts_all = segments["start"].to_numpy()
        ends_all = segments["end"].to_numpy()
        for chrom, sel in segments.groupby("chrom", sort=False).indices.items():
            ref = self._by_chrom.get(chrom)
            if ref is None:
                continue
            r_start, r_end, r_cummax, r_gene = ref
            a_start, a_end = starts_all[sel], ends_all[sel]
            # candidates: ref start < a_end, and not ruled out by a running max end <= a_start
            hi = np.searchsorted(r_start, a_end, side="left")
            lo = np.searchsorted(r_cummax, a_start, side="right")
            counts = np.maximum(hi - lo, 0)
            total = int(counts.sum())
            if not total:
                continue
            rep = np.repeat(np.arange(len(sel)), counts)
            b = lo[rep] + np.arange(total) - np.repeat(np.cumsum(counts) - counts, counts)
            size = np.minimum(a_end[rep], r_end[b]) - np.maximum(a_start[rep], r_start[b])
            keep = size > 0
            a_parts.append(sel[rep[keep]])
            g_parts.append(r_gene[b[keep]])
            o_parts.append(size[keep])
        if not a_parts:
            empty = np.empty(0, dtype=np.int64)
            return empty, empty, empty
        a_idx = np.concatenate(a_parts)
        order = np.argsort(a_idx, kind="stable")  # back to segment order, ref order within a segment
        return a_idx[order], np.concatenate(g_parts)[order], np.concatenate(o_parts)[order]


_GENE_INDEX: GeneIndex | None = None


def _init_worker(gene_index: GeneIndex) -> None:
    """Give each worker process the gene index once, instead of with every task."""
    global _GENE_INDEX  # noqa: PLW0603
    _GENE_INDEX = gene_index


def get_gene_cnv_array(segments: pd.DataFrame, gene_index: GeneIndex, ploidy: int = 2) -> np.ndarray:
    """Assign a CN to each gene overlapped by the segments.

    For repeat entries, keep largest overlap with greatest non-neutral CN.

    Returns:
        int32 array in gene_index order, NOT_SEEN for genes without overlap

    """
    a_idx, gene_ids, sizes = gene_index.overlaps(segments)
    cns = segments["cn"].to_numpy()[a_idx]
    best_size = [-1] * gene_index.n_genes
    best_cn = [NOT_SEEN] * gene_index.n_genes
    for gene, size, cn in zip(gene_ids.tolist(), sizes.tolist(), cns.tolist(), strict=True):
        cur = best_size[gene]
        if (
            cur < 0
            or (size > cur and cn != ploidy)
            or (size == cur and abs(ploidy - cn) > abs(ploidy - best_cn[gene]))
        ):
            best_size[gene] = size
            best_cn[gene] = cn
    return np.array(best_cn, dtype=np.int32)


def gistic_from_cn(cn: np.ndarray, ploidy: np.ndarray, high_gain: int) -> np.ndarray:
    """Convert CN to GISTIC-style values relative to each sample's ploidy.

    high_gain is a relative threshold above ploidy; CN ABOVE ploidy + high_gain is an
    amplification (2). For example high_gain 4 with ploidy 3 makes the threshold 7.
    """
    return np.select(
        [cn > ploidy + high_gain, cn > ploidy, cn == ploidy, cn > 0], [2, 1, 0, -1], default=-2
    )


def _segments(chrom: pd.Series, start: pd.Series, end: pd.Series, cn: pd.Series) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": chrom.astype(str).to_numpy(),
            "start": pd.to_numeric(start).astype(np.int64).to_numpy(),
            "end": pd.to_numeric(end).astype(np.int64).to_numpy(),
            "cn": pd.to_numeric(cn).astype(np.int64).to_numpy(),
        }
    )


def _read_str_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    """Read a TSV with every field as an unmodified string."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def _seg_text(seg: pd.DataFrame) -> str:
    """Seg rows as tab separated text; the first column is expected to be the sample ID."""
    return seg.to_csv(sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE, lineterminator="\n")


def read_and_process_ctrlfreec_seg(seg_fname: Path, sample_id: str) -> str:
    """Process ControlFreeC seg file into seg rows with the cBio sample ID and no leading chr."""
    seg = _read_str_tsv(seg_fname)
    seg.iloc[:, 0] = sample_id
    seg.iloc[:, 1] = seg.iloc[:, 1].str.removeprefix("chr")
    return _seg_text(seg)


def get_ctrlfreec_ploidy(cfree_info_fname: Path) -> int:
    """Process ControFreeC Info file for Ploidy information.

    Used to set a default value for CN if sample does not have a one for a specific gene
    """
    ploidy = 2
    with cfree_info_fname.open() as info:
        for data in info:
            datum = data.rstrip("\n")
            if datum.startswith("Output_Ploidy"):
                try:
                    (_key, value) = datum.split("\t")
                    ploidy = int(value)
                except ValueError:
                    logger.warning("Entry %r in %s could not be split", datum, cfree_info_fname)
                break
    return ploidy


def read_and_process_ctrlfreec_pval(pval_fname: Path) -> pd.DataFrame:
    """Filter ControlFreeC CNVs by p-values (both <= PVALUE) and return segments."""
    df = _read_str_tsv(pval_fname)
    wilcox = pd.to_numeric(df["WilcoxonRankSumTestPvalue"], errors="coerce")
    ks = pd.to_numeric(df["KolmogorovSmirnovPvalue"], errors="coerce")
    bed = df.iloc[:, :4][(wilcox <= PVALUE) & (ks <= PVALUE)]  # "NA" becomes NaN, so is dropped
    return _segments(bed.iloc[:, 0], bed.iloc[:, 1], bed.iloc[:, 2], bed.iloc[:, 3])


def read_and_process_cnvkit_theta2_cns(cns_fname: Path) -> pd.DataFrame:
    """Process CNS file from CNVkit/theta2 and return segments, with leading chr stripped."""
    df = _read_str_tsv(cns_fname, usecols=[0, 1, 2, 6])
    return _segments(df.iloc[:, 0].str.removeprefix("chr"), df.iloc[:, 1], df.iloc[:, 2], df.iloc[:, 3])


def read_and_process_gatk_cnv(gatk_seg_fname: Path, sample_id: str) -> tuple[pd.DataFrame, str]:
    """Process GATK CNV seg file.

    Convert MEAN_LOG2_COPY_RATIO to CN (assuming ploidy 2), skip chrM, and strip leading chr.
    The input is already a seg file, so seg rows are returned as well.

    Returns:
        Segments and seg file text

    """
    # skip interval list-style headers by finding the CONTIG header line
    with gatk_seg_fname.open() as f:
        header_line = next((i for i, line in enumerate(f) if line.startswith("CONTIG")), None)
    if header_line is None:
        err = f"No CONTIG header line found in {gatk_seg_fname}"
        raise ValueError(err)
    try:
        df = _read_str_tsv(gatk_seg_fname, header=None, skiprows=header_line + 1)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=range(6), dtype=str)
    df = df[~df[0].str.startswith("chrM")]
    df[0] = df[0].str.removeprefix("chr")
    # ratio is typically log2(cn/ploidy)
    cn = np.rint(2.0 ** df.iloc[:, -2].astype(float) * 2)
    segments = _segments(df[0], df[1], df[2], pd.Series(cn, index=df.index))
    seg = df.iloc[:, :5].copy()
    seg.insert(0, "ID", sample_id)
    return segments, _seg_text(seg)


def process_sample(
    cbio_sample: str,
    cnv_samp_attr: dict[str, str],
    info_samp_attr: dict[str, str] | None,
    seg_samp_attr: dict[str, str] | None,
) -> SampleResult:
    """Process one sample's CNV data (runs in a worker process)."""
    ploidy = 2
    try:
        ftype = cnv_samp_attr["manifest_ftype"]
        if ftype == "ctrlfreec_pval" or ftype.endswith("_calls"):
            logger.debug("Processing %s %s file", cbio_sample, ftype)
            if not seg_samp_attr:
                err = f"No seg file for {cbio_sample}, cannot process CNV data"
                raise FileNotFoundError(err)
            seg_text = read_and_process_ctrlfreec_seg(
                Path(seg_samp_attr["manifest_ftype"]) / seg_samp_attr["fname"], cbio_sample
            )
            # ControlFreeC info file is optional, but if it exists, use it to set ploidy
            if info_samp_attr:
                ploidy = get_ctrlfreec_ploidy(Path(info_samp_attr["manifest_ftype"]) / info_samp_attr["fname"])
            else:
                logger.warning("No info file for %s, using default ploidy of 2", cbio_sample)
            cnv_path = Path(ftype) / cnv_samp_attr["fname"]
            # only a slight difference in cnvkit and pval inputs, the rest should be the same
            if ftype == "ctrlfreec_pval":
                segments = read_and_process_ctrlfreec_pval(cnv_path)
            else:
                segments = read_and_process_cnvkit_theta2_cns(cnv_path)
        else:
            logger.debug("Processing %s GATK CNV file", cbio_sample)
            segments, seg_text = read_and_process_gatk_cnv(Path(ftype) / cnv_samp_attr["fname"], cbio_sample)
        raw_cn = get_gene_cnv_array(segments, _GENE_INDEX, ploidy)
    except Exception as e:
        err = f"{e} processing {cbio_sample}"
        raise RuntimeError(err) from e
    return SampleResult(raw_cn, ploidy, seg_text)


def output_wide_format_tables(
    raw: np.ndarray,
    ploidy: np.ndarray,
    samples: list[str],
    genes: np.ndarray,
    raw_filename: Path,
    gistic_filename: Path,
    high_gain: int,
) -> None:
    """Write the raw and GISTIC wide tables (Hugo_Symbol by sample).

    Only genes overlapped by a CNV in at least one sample are written. Samples without a CNV for
    a written gene get their own ploidy (raw) or 0 (GISTIC).

    Args:
        raw: genes x samples CN matrix, NOT_SEEN where there was no overlap
        ploidy: ploidy per sample
        samples: sample names (columns)
        genes: gene names (rows)
        raw_filename: Output file for raw CN
        gistic_filename: Output file for GISTIC-style values
        high_gain: Value over ploidy that results in a GISTIC-style high gain designation

    """
    header = "\t".join(["Hugo_Symbol", *samples]) + "\n"
    rows = np.flatnonzero((raw != NOT_SEEN).any(axis=1))
    with raw_filename.open("w") as raw_out, gistic_filename.open("w") as gistic_out:
        raw_out.write(header)
        gistic_out.write(header)
        for start in range(0, len(rows), WRITE_CHUNK_ROWS):
            idx = rows[start : start + WRITE_CHUNK_ROWS]
            block = raw[idx]
            cn = np.where(block == NOT_SEEN, ploidy, block)
            gistic = gistic_from_cn(cn, ploidy, high_gain)
            for out, values in ((raw_out, cn), (gistic_out, gistic)):
                out.write(
                    "".join(
                        gene + "\t" + "\t".join(map(str, row)) + "\n"
                        for gene, row in zip(genes[idx].tolist(), values.tolist(), strict=True)
                    )
                )


def prioritize_cnvs(
    cbio_tbl: str, config_data: dict[str, list[str]]
) -> dict[str, dict[str, dict[str, dict[str, str]]]]:
    """Prioritize CNVs based on config data.

    Args:
        cbio_tbl: cBio etl filename
        config_data: CNVs defined by cnv_priority

    Returns:
        Dict of cnv types with prioritized file names and manifest file types.

    """
    # lists and dicts of CNV types
    cnv_ftypes: dict[str, dict[str, dict[str, dict[str, str]]]] = {"cnv": {}, "seg": {}, "info": {}}
    with Path(cbio_tbl).open() as f:
        cbio_reader = csv.reader(f, delimiter="\t")

        header = next(cbio_reader)
        # get key fields from header
        cbio_project_idx = header.index("cbio_project")
        cbio_sample_idx = header.index("cbio_sample_name")
        etl_ftype_idx = header.index("etl_file_type")
        etl_exp_idx = header.index("etl_experiment_strategy")
        manifest_ftype_idx = header.index("file_type")
        fname_idx = header.index("file_name")
        for entry in cbio_reader:
            project = entry[cbio_project_idx]
            cbio_sample = entry[cbio_sample_idx]
            etl_ftype = entry[etl_ftype_idx]
            etl_exp = entry[etl_exp_idx]
            manifest_ftype = entry[manifest_ftype_idx]
            fname = entry[fname_idx]

            if etl_ftype not in cnv_ftypes:
                continue
            samples = cnv_ftypes[etl_ftype].setdefault(project, {})
            if cbio_sample not in samples:
                samples[cbio_sample] = {"fname": fname, "manifest_ftype": manifest_ftype}
            # If the etl_ftype is already set for that cbio_sample
            # check if the new manifest_ftype is higher priority
            elif config_data[etl_exp].index(manifest_ftype) < config_data[etl_exp].index(
                samples[cbio_sample]["manifest_ftype"]
            ):
                samples[cbio_sample] = {"fname": fname, "manifest_ftype": manifest_ftype}
    return cnv_ftypes


def main() -> None:
    """Convert GATK and ControlFreeC CNV data to cBio format."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    parser = argparse.ArgumentParser(
        description="Merges files in gene <tab> entrez id <tab> copy number format into a genes-by-sample copy number table"
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        dest="table",
        help="Table with cbio project, kf bs ids, cbio IDs, and file names",
    )
    parser.add_argument(
        "-j",
        "--config",
        action="store",
        dest="config_file",
        help="json config file with data types and data locations",
    )
    args = parser.parse_args()
    tool_dir = Path(__file__).resolve().parent.parent

    with Path(args.config_file).open() as f:
        config_data = json.load(f)
    config_data: dict = resolve_config_paths(config_data, str(tool_dir))
    logger.info("Prioritizing CNV calls base on %s", config_data["cnv_priority"])
    # subset cnv and seg data using priority from config file
    prioritized_cnv_meta = prioritize_cnvs(args.table, config_data["cnv_priority"])

    gene_index = GeneIndex(Path(config_data["bed_genes"]))
    high_gain: int = config_data["cnv_high_gain"]
    out_dir = Path("merged_cnvs")
    out_dir.mkdir(exist_ok=True, parents=True)

    with ProcessPoolExecutor(initializer=_init_worker, initargs=(gene_index,)) as executor:
        # cnv meta should have ALL IDs to be processed
        for project, samp_data in prioritized_cnv_meta["cnv"].items():
            samples = list(samp_data)
            # genes x samples, filled as results arrive; ploidy per sample is the default CN
            raw = np.full((gene_index.n_genes, len(samples)), NOT_SEEN, dtype=np.int32)
            ploidy = np.empty(len(samples), dtype=np.int32)
            # info and seg are non-GATK specific
            info_meta = prioritized_cnv_meta["info"].get(project, {})
            seg_meta = prioritized_cnv_meta["seg"].get(project, {})
            logger.info("Processing %d %s samples", len(samples), project)
            try:
                with (out_dir / f"{project}.merged_seg.txt").open("w") as out_seg:
                    out_seg.write(SEG_HEADER + "\n")
                    # map yields in sample order, so columns are deterministic
                    results = executor.map(
                        process_sample,
                        samples,
                        [samp_data[s] for s in samples],
                        [info_meta.get(s) for s in samples],
                        [seg_meta.get(s) for s in samples],
                        chunksize=max(1, len(samples) // ((os.cpu_count() or 1) * 4)),
                    )
                    for j, result in enumerate(results):
                        raw[:, j] = result.raw_cn
                        ploidy[j] = result.ploidy
                        out_seg.write(result.seg_text)
                        if (j + 1) % 100 == 0:
                            logger.info("Processed %d of %d %s samples", j + 1, len(samples), project)
            except KeyboardInterrupt:
                logger.error("Keyboard interrupt, exiting")
                executor.shutdown(wait=False, cancel_futures=True)
                sys.exit(1)
            except Exception:
                logger.exception("Failed processing %s CNV data", project)
                executor.shutdown(wait=False, cancel_futures=True)
                sys.exit(1)
            # Write raw and GISTIC results to wide format
            logger.info("Writing %s raw and GISTIC CNV data to wide format", project)
            output_wide_format_tables(
                raw=raw,
                ploidy=ploidy,
                samples=samples,
                genes=gene_index.genes,
                raw_filename=out_dir / f"{project}.predicted_cnv.txt",
                gistic_filename=out_dir / f"{project}.discrete_cnvs.txt",
                high_gain=high_gain,
            )
            logger.info("Finished processing %s CNV data", project)


if __name__ == "__main__":
    main()
