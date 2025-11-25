# Publication Study Load chdm_gse230168: Molecular Characterization of Chordomas (GSE230168, publication)
As part of the Open Chordoma effort, this dataset was loaded from GEO.
Unfortunately, as of 11/25/2025, the author has only made available counts to GENCODE 39, but geenrated using `Rsubread`.
Due to patient confidentiality we cannot reprocess BAMs and running RSEM to make the dataset load compatible with current Open Chordoma studies has not happened.
Therefore, this will exist as a standalone.

## Data processing
### Inputs
- [Counts](./ChorCountsv39.csv.gz): `Rsubread` counts from author, generated from STAR alignment using a GENCODE 39 reference
- [Metadata](https://github.com/SBaluszek/chordoma_RNA_met/blob/main/spl_clust.csv): Cluter metadata file
### Formatting
- For `data_clinical_*`, data was copied into an Excel sheet, and formatted to give the most basic info from what was provided
- Input `Counts` had two samples not in `Metadata` file, so they were dropped:
   ```sh
   cut -f 1-2,4-28,30- ChorCountsv39.csv -d "," > ChorCountsv39_subset.csv
   ```
- Custom [parse_counts.py](./parse_counts.py) script was used to convert ENSEMBL IDs to Hugo gene symbols via biomart. Produces TSV `counts_w_gene_sym.tsv`
   ```sh
   python3 ~/tools/kf-cbioportal-etl/PUBLICATIONS/chdm_gse230168/parse_counts.py ChorCountsv39_subset.csv 2> parse.log
   ```
- Custom [collapse_and_zscore_cts](./collapse_and_zscore_cts.py) script run to create collapsed (repeat gene symbol, highest mean expr) and cohort zscore
   ```sh
   python3 ~/tools/kf-cbioportal-etl/PUBLICATIONS/chdm_gse230168/collapse_and_zscore_cts.py -t counts_w_gene_sym.tsv
   ```
- Subdir made with harlinks to RNA data, case lists and supporting metadata "hand made"
```
chdm_gse230168
├── case_lists
│   ├── cases_RNA_Seq_v2_mRNA.txt
│   └── cases_all.txt
├── data_clinical_patient.txt
├── data_clinical_sample.txt
├── data_rna_seq_v2_mrna.txt
├── data_rna_seq_v2_mrna_median_Zscores.txt
├── meta_clinical_patient.txt
├── meta_clinical_sample.txt
├── meta_rna_seq_v2_mrna.txt
├── meta_rna_seq_v2_mrna_median_Zscores.txt
└── meta_study.txt
```
## Data load
Uploaded to `s3://kf-strides-232196027141-cbioportal-studies/studies/chdm_gse230168/`.
Study loaded following typical git commit process outlined [here](https://github.com/d3b-center/aws-infra-pedcbioportal-import?tab=readme-ov-file#importdelete-triggers)