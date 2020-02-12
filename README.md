# Outline on ETL for converting data from cavatica and data service to pedcbioportal format
## Software Prerequisites

+ `python3` v3.5.3+
  + `sevenbridges` package (https://sevenbridges-python.readthedocs.io/en/latest/installation/)
  + `numpy`, `pandas`, `scipy`
+ `bedtools` (https://bedtools.readthedocs.io/en/latest/content/installation.html)
+ Seven bridges credentials file created

## Starting file inputs
+ `dx index file`

This tsv file is a dictionary for converting existing diagnoses into cbioportal studies. The portal has pre-defined disease types
that will likely overlap with your sample's existing diagnosis.  The last column is a long form description that you can customize. The ETL is designed to either categorized your data
inputs into seprate studies depending on thier diagnosis, like this:

```
CBTTC Diagnosis	MB cbioportal id	MB cbioportal name
Adenoma	aptad	Atypical Pituitary Adenoma
Atypical Teratoid Rhabdoid Tumor (ATRT)	atrt	Atypical Teratoid/Rhabdoid Tumor
Brainstem glioma- Diffuse intrinsic pontine glioma	dipg	Diffuse Intrinsic Pontine Glioma
Chordoma	chdm	Chordoma
Choroid plexus carcinoma	cpc	Choroid Plexus Carcinoma
Choroid plexus papilloma	cpp	Choroid Plexus Papilloma
Craniopharyngioma	cranio	Craniopharyngioma
Dysembryoplastic neuroepithelial tumor (DNET)	dnt	Dysembryoplastic Neuroepithelial Tumor
Dysplasia/Gliosis	other	Other
Ependymoma	epmt	Ependymomal Tumor
Ganglioglioma	gng	Ganglioglioma
Germinoma	gmn	Germinoma
Glial-neuronal tumor NOS	gnos	Glioma, NOS
Gliomatosis Cerebri	gnos	Glioma, NOS
High-grade glioma/astrocytoma (WHO grade III/IV)	phgg	Pediatric High Grade Gliomas
Langerhans Cell histiocytosis	lch	Langerhans Cell Histiocytosis
Low-grade glioma/astrocytoma (WHO grade I/II)	plgg	Pediatric Low Grade Gliomas
Malignant peripheral nerve sheath tumor (MPNST)	mpnst	Malignant Peripheral Nerve Sheath Tumor
Medulloblastoma	mbl	Medulloblastoma
Meningioma	mng	Meningioma
Metastatic secondary tumors	other	Other
Neurocytoma	cnc	Central Neurocytoma
Neurofibroma/Plexiform	nfib	Neurofibroma
Non-germinomatous germ cell tumor	bmgct	Mixed Germ Cell Tumor
Oligodendroglioma	odg	Oligodendroglioma
Other	other	Other
Not Reported	other	Other
Pineoblastoma	pbl	Pineoblastoma
Primary CNS lymphoma	pcnsl	Primary DLBCL of the central nervous system
Schwannoma	schw	Schwannoma
Supratentorial or Spinal Cord PNET	pnet	Primitive Neuroectodermal Tumor
Teratoma	tt	Teratoma
Cavernoma	other	Other
Ewing's Sarcoma	es	Ewing Sarcoma
Ganglioneuroblastoma	gnbl	Ganglioneuroblastoma
Hemangioblastoma	hmbl	Hemangioblastoma
Neuroblastoma	nbl	Neuroblastoma
Rhabdomyosarcoma	rms	Rhabdomyosarcoma
Sarcoma	sarcnos	Sarcoma, NOS
Subependymal Giant Cell Astrocytoma (SEGA)	epm	Ependymoma
Pilocytic Astrocytoma	lggnos	Low-Grade Glioma, NOS
```

OR combine into a single large study, with the cbioportal disease being of a more general oncotree value:
```
CBTTC Diagnosis	MB cbioportal id	MB cbioportal name
Adenoma	brain	Pediatric Brain Tumor Atlas
Atypical Teratoid Rhabdoid Tumor (ATRT)	brain	Pediatric Brain Tumor Atlas
Brainstem glioma- Diffuse intrinsic pontine glioma	brain	Pediatric Brain Tumor Atlas
Chordoma	brain	Pediatric Brain Tumor Atlas
Choroid plexus carcinoma	brain	Pediatric Brain Tumor Atlas
Choroid plexus papilloma	brain	Pediatric Brain Tumor Atlas
Craniopharyngioma	brain	Pediatric Brain Tumor Atlas
Dysembryoplastic neuroepithelial tumor (DNET)	brain	Pediatric Brain Tumor Atlas
Dysplasia/Gliosis	brain	Pediatric Brain Tumor Atlas
Ependymoma	brain	Pediatric Brain Tumor Atlas
Ganglioglioma	brain	Pediatric Brain Tumor Atlas
Germinoma	brain	Pediatric Brain Tumor Atlas
Glial-neuronal tumor NOS	brain	Pediatric Brain Tumor Atlas
Gliomatosis Cerebri	brain	Pediatric Brain Tumor Atlas
High-grade glioma/astrocytoma (WHO grade III/IV)	brain	Pediatric Brain Tumor Atlas
Langerhans Cell histiocytosis	brain	Pediatric Brain Tumor Atlas
Low-grade glioma/astrocytoma (WHO grade I/II)	brain	Pediatric Brain Tumor Atlas
Malignant peripheral nerve sheath tumor (MPNST)	brain	Pediatric Brain Tumor Atlas
Medulloblastoma	brain	Pediatric Brain Tumor Atlas
Meningioma	brain	Pediatric Brain Tumor Atlas
Metastatic secondary tumors	brain	Pediatric Brain Tumor Atlas
Neurocytoma	brain	Pediatric Brain Tumor Atlas
Neurofibroma/Plexiform	brain	Pediatric Brain Tumor Atlas
Non-germinomatous germ cell tumor	brain	Pediatric Brain Tumor Atlas
Oligodendroglioma	brain	Pediatric Brain Tumor Atlas
Other	brain	Pediatric Brain Tumor Atlas
Not Reported	brain	Pediatric Brain Tumor Atlas
Pineoblastoma	brain	Pediatric Brain Tumor Atlas
Primary CNS lymphoma	brain	Pediatric Brain Tumor Atlas
Schwannoma	brain	Pediatric Brain Tumor Atlas
Supratentorial or Spinal Cord PNET	brain	Pediatric Brain Tumor Atlas
Teratoma	brain	Pediatric Brain Tumor Atlas
Cavernoma	brain	Pediatric Brain Tumor Atlas
Ewing's Sarcoma	brain	Pediatric Brain Tumor Atlas
Ganglioneuroblastoma	brain	Pediatric Brain Tumor Atlas
Hemangioblastoma	brain	Pediatric Brain Tumor Atlas
Neuroblastoma	brain	Pediatric Brain Tumor Atlas
Rhabdomyosarcoma	brain	Pediatric Brain Tumor Atlas
Sarcoma	brain	Pediatric Brain Tumor Atlas
Subependymal Giant Cell Astrocytoma (SEGA)	brain	Pediatric Brain Tumor Atlas
Pilocytic Astrocytoma	brain	Pediatric Brain Tumor Atlas
```
+ Cavativa file manifest(s)

This is a manifest of all files to loaded onto the portal.  Currently, likely data sources are consensus calls, ControlFreeC p value cnv file outputs, rsem gene counts, and annotated fusion output from annoFuse. You may have to pull from multiple cavatica projects; a helper script exists to merge all manifests if needed which is explained later.

+ Data processing config file

This is a json formatted file that has tool paths, reference paths, and run time params.  An example is given in `REFS/data_processing_config.json`

+ Metadata processing config file

This is a json config file with file descriptions and case lists required by the cbioportal.

## Initial setup
Here, we parse the cavatica manifests, get clinical info from the data service, and create cbioportal datasheets with cbioportal-compatible names.

### scripts/1a_merge_cavatica_manifests.py
Run if more than one cavatica file manifest was used/needed
```
usage: 1a_merge_cavatica_manifests.py [-h] [-m MANIFEST]

Get all relevant analyzed file outputs from cavatica manifest.

optional arguments:
  -h, --help            show this help message and exit
  -m MANIFEST, --manifest MANIFEST
                        cavatica manifest file csv list
```

### scripts/1b_query_file_manifest.py
This will create a master sheet bs ids, analyte type, project location, and file names.
```
usage: 1b_query_file_manifest.py [-h] [-o OUTPUT] [-m MANIFEST] [-b BLACKLIST]
                                 [-c CONFIG_FILE]

Get all relevant analyzed file outputs from cavatica manifest.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output basename name
  -m MANIFEST, --manifest MANIFEST
                        cavatica csv manifest file, may be from step 1a
  -b BLACKLIST, --blacklist BLACKLIST
                        Optional bs id blacklist
  -c CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
```
Output will look something like this:
```
T/CL BS ID  Norm BS ID  Task Name Task ID Analyte Type  Relevant Outputs  Source Project
BS_VXDGXQKZ BS_D48QXYW6 DNA_TASK  BS_VXDGXQKZBS_D48QXYW6  DNA 1b4fe9f0-bcd8-48d5-85f8-1a618ed29e50.consensus_somatic.vep.maf,7cda5662-1287-43ac-8f4b-5fdc9181cf2e.controlfreec.CNVs.p.value.txt kfdrc-harmonization/pbta-consensus-calls,kfdrc-harmonization/sd-bhjxbdqk-08
```
Data processing config file is used, with relevant fields `dna_ext_list` and `rna_ext_list`

### scripts/2_query_ds_by_bs_id.py
This must be run on the chop network. It grabs all the clinical data from the data service and creates a table for tumor data (`tum_bs_ds_info.txt`), and one for normal data (`tum_bs_ds_info.txt`).  
```
usage: 2_query_ds_by_bs_id.py [-h] [-u URL] [-c CAV] [-j CONFIG_FILE]

Script to walk through data service and grab all relevant meetadata by bs id.

optional arguments:
  -h, --help            show this help message and exit
  -u URL, --kf-url URL  Kids First data service url, i.e. https://kf-api-
                        dataservice.kidsfirstdrc.org/
  -c CAV, --cavatica CAV
                        file with task info from cavatica (see step 1)
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
```
Example tumor file output:
```
BS_ID	PT_ID	external_aliquot_id	analyte_type	source_text_tissue_type	source_text_tumor_descriptor	composition	external_sample_id	external_id	gender	ethnicity	race	source_text_diagnosis	source_text_tumor_location	age_at_event_days
BS_J8EK6RNF	PT_KBFM551M	556480	DNA	Tumor	Progressive Disease Post-Mortem	Solid Tissue	7316-3237-T.WGS	P-04	Male	Not Hispanic or Latino	White	Brainstem glioma- Diffuse intrinsic pontine glioma	Brain Stem	3285
BS_5968GBGT	PT_KTRJ8TFY	SF5606-19	DNA	Tumor	Progressive Disease Post-Mortem	Solid Tissue	7316-3230-T-SF5606-19.WGS	P-13	Male	Hispanic or Latino	Reported Unknown	Brainstem glioma- Diffuse intrinsic pontine glioma	Brain Stem	1825
```
Example normal file output:
```
BS_ID	PT_ID	external_aliquot_id	analyte_type	source_text_tissue_type	source_text_tumor_descriptor	composition	external_sample_id	external_id	gender	ethnicity	race	source_text_diagnosis	source_text_tumor_location	age_at_event_days
BS_HJ7HYZ7N	PT_KBFM551M	556481	DNA	Normal	NULL	Solid Tissue	7316-3237-N-Brain.WGS	P-04	Male	Not Hispanic or Latino	White
BS_SNRF1RKC	PT_KTRJ8TFY	A09700	DNA	Normal	NULL	Peripheral Whole Blood	7316-3223-N-A09687.WGS	P-13	Male	Hispanic or Latino	Reported Unknown
```
Data processing config file is used, with relevant field `threads` as python requests queries are multithreaded

### scripts/3a_create_data_sheets.py
This file creates the clinical information files required for cbioportal. A directory called `datasheets` will be created with one sub-directory for each unique value from the `MB cbioportal id` in the `dx index sheet`. It uses a helper script called `scripts/sample_id_builder_helper.py` to for cbioportal-friendly names. Edit that script if a different set of information from the ds_info sheets is needed.
```
usage: 3a_create_data_sheets.py [-h] [-c CAV] [-t TUMOR] [-n NORMAL]
                                [-j CONFIG_FILE] [-s CL_SUPP]

Convert metadata info in data patient and sample sheets for cbio portal

optional arguments:
  -h, --help            show this help message and exit
  -c CAV, --cavatica CAV
                        file with task info from cavatica (see step 1)
  -t TUMOR, --tumor TUMOR
                        file with tumor metadata (see step 2)
  -n NORMAL, --normal NORMAL
                        file with normal metadata (see step 2)
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
  -s CL_SUPP, --cell-line-supplement CL_SUPP
                        supplemental file with cell line meta data - bs
                        id<tab>type. optional
```
Data processing config file is used, with relevant fields:
+ `dx_tbl_fn`: name of dx index file.  Should be in cwd.
+ `ds_pt_desc`: patient file headers, see example config
+ `ds_samp_desc`: sample file headers, see example config
+ `dna_samp_id_style`: Naming style defintion to use for DNA tumor samples in `scripts/sample_id_builder_helper.py`
+ `dna_norm_id_style`: Naming style defintion to use for DNA normal samples in `scripts/sample_id_builder_helper.py`
+ `rna_samp_id_style`: Naming style defintion to use for DNA tumor samples in `scripts/sample_id_builder_helper.py`
+ `rna_flag`: 1 or 0 value if there is or is not RNA data

### scripts/3b_create_cbio_id_fname_tbl.py
This is another organizing script that collates the cavatica data from step 1b and the datasheets created from step 3a. Outputs dfile with name `cbio_id_fname_table.txt`
```
usage: 3b_create_cbio_id_fname_tbl.py [-h] [-c CAV]

Create temp table with cbio ids, bs ids, and file locations and types for
downstream file merging. Should be run one level above datasheets dir.

optional arguments:
  -h, --help            show this help message and exit
  -c CAV, --cavatica CAV
                        file with task info from cavatica (see step 1b)
```
Example output:
```
Cbio_project	T_CL_BS_ID	Norm_BS_ID	File_Type	Cbio_Tumor_Name	Cbio_Matched_Normal_Name	File_Name
brain	BS_AQMKA8NC	BS_CTEM6SYF	maf	7316-2577-T-463571	BS_CTEM6SYF	f6d71edb-b7e0-4611-8070-b987dff88feb.consensus_somatic.vep.maf
brain	BS_AQMKA8NC	BS_CTEM6SYF	cnv	7316-2577-T-463571	BS_CTEM6SYF	f9100849-4049-4d59-af1b-005666ac896e.controlfreec.CNVs.p.value.txt
brain	BS_FEPRNEXX	NA	rsem	7316-2577-T-463571	NA	2c98cd42-5c3f-4a7e-8e53-7fce6b317222.rsem.genes.results.gz
```

### scripts/4_merge_maf.py
*Prerequisite: Download all maf files to be merged into a directory, using the manifest and wget, sbg api, or xargs and sb cli.* Merges maf files by cbio disease type.  Will remove calls of certain categories as outlined [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#mutation-data)
```
usage: 4_merge_maf.py [-h] [-t TABLE] [-i HEADER] [-m MAF_DIR]
                      [-j CONFIG_FILE]

Merge and filter mafs using cavatica task info and datasheets.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file
                        names
  -i HEADER, --header HEADER
                        File with maf header only. Example in REFS/ dir
  -m MAF_DIR, --maf-dir MAF_DIR
                        maf file directory
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
```
Data processing config file is used, with relevant fields:
+ `cpus`: Number of files to process concurrently using `ProcessPoolExecutor`.  8 recommended
Output created: `merged_mafs` directory with maf files named based on cbio disease type.

### scripts/5a_cnv_genome2gene.py
*Prerequisite: Download all ControlFreeC p value output files to be converted into a directory, using the manifest and wget, sbg api, or xargs and sb cli.* Converts copy number genome coordinate calls to gene copy number calls
```
usage: 5a_cnv_genome2gene.py [-h] [-d CNV_DIR] [-j CONFIG_FILE]

Convert controlFreeC cnv genome coords to gene level

optional arguments:
  -h, --help            show this help message and exit
  -d CNV_DIR, --cnv-dir CNV_DIR
                        cnv as genome coords file directory
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
```
Data processing config file is used, with relevant fields:
+ `cnv_min_len`: Minimum copy number length (in base pairs) to be considered a candidate.  50000 recommended.
+ [`dna_ext_list`][`copy_number`]: File extension of inputs with no leading `.`, i.e. `controlfreec.CNVs.p.value.txt`
+ `cpus`: Number of files to process concurrently using `ProcessPoolExecutor`.  8 recommended
+ `bedtools`: Location of bedtools executable.  Can just be `bedtools` if in PATH
+ `cp_only_script`: Location of helper perl script that merges bedtools results. Script name `get_cbio_copy_only_num.pl`
+ `bed_genes`: Bed file with gene coordinates.  Example:
```
17	7581964	7584072	AC113189.5
17	7583529	7592789	MPDU1
17	7588178	7590170	SOX15
17	7591230	7614871	FXR2
17	7613946	7633383	SHBG
17	7626234	7627876	SAT2
17	7646627	7657768	ATP1B2
17	7661779	7687550	TP53
```
+ `hugo_tsv`: tsv file with HUGO gene symbols and entrez IDs for perl copy script
Output created: `converted_cnvs` directory with converted files, extension `CNVs.Genes.copy_number`, with HUGO gene symbol, entrez ID,and raw copy number.  Example:
```
DCC     1630    4
DCD     117159  3
DCDC1   341019  3
DCHS1   8642    3
DCHS2   54798   3
DCK     1633    3
DCLK2   166614  3
```

### scripts/5b_merge_cnv.py
Merge files generated from step 5a into a genes-by-sample raw copy number matrix
```
usage: 5b_merge_cnv.py [-h] [-t TABLE] [-n CNV_DIR] [-j CONFIG_FILE]

Merge cnv files using cavatica task info and datasheets.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file
                        names
  -n CNV_DIR, --cnv-dir-gene CNV_DIR
                        cnv as gene file directory
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
```
Data processing config file is used, with relevant fields:
+ [`dna_ext_list`][`copy_number`]: File extension of inputs (used in 5a) with no leading `.`, i.e. `controlfreec.CNVs.p.value.txt`
+ `cpus`: Number of files to process concurrently using `ProcessPoolExecutor`.  8 recommended
Output created: `merged_cnvs` directory with merged raw copy number files named based on cbio disease type. Will have extension `predicted_cnv.txt`

### scripts/5c_cnv_discrete.py
*Prerequisite: A cavatica file manifest of ControlFreeC `info.txt` files* Converts raw copy number to gistic-style discrete values. Uses info.txt files to adjust for calculated ploidy
```
usage: 5c_cnv_discrete.py [-h] [-d MERGED_CNV_DIR] [-j CONFIG_FILE]
                          [-m MANIFEST] [-t TABLE] [-p PROFILE]

Convert merged cnv values to discrete coded values.

optional arguments:
  -h, --help            show this help message and exit
  -d MERGED_CNV_DIR, --merged_cnv_dir MERGED_CNV_DIR
                        merged cnv dir
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
  -m MANIFEST, --info_manifest MANIFEST
                        cavatica cfree info file manifest
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file
                        names
  -p PROFILE, --profile PROFILE
                        cavatica profile name. requires
                        .sevenbridges/credentials file be present
```
Data processing config file is used, with relevant fields:
+ `cnv_high_gain`: Value cutoff to assign gistic high amplification gain value of 2, after ploidy adjustment.  For example, a value of `4` (recommended) would be equivalent of setting anything with value greater than (ploidy + 4) to a value of 2, repesenting high gain
+ `threads`: Number of files to process concurrently using `ThreadPoolExecutor`. 40 recommended
Output created: In `merged_cnvs` directory with merged discrete copy number files named based on cbio disease type. Will have extension `discrete_cnvs.txt`

### scripts/6_merge_rename_rsem.py
*Prerequisite: Download all rsem files to be merged into a directory, using the manifest and wget, sbg api, or xargs and sb cli.* Merges rsem files into a gene-by-sample matrix.  Repeat HUGO symbols are resolved by taking the row with the highest mean expression and dropping the rest. Also creates a merged z score table.  Z scores are calculated before breaking up into individual cbio disease files; FPKM counts have a pseudocount added, then are log2 transformed before z score is calculated
```usage: 6_merge_rename_rsem.py [-h] [-t TABLE] [-r RSEM_DIR] [-j CONFIG_FILE]

Merge rsem files using cavatica file info.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file
                        names
  -r RSEM_DIR, --rsem-dir RSEM_DIR
                        rsem file directory
```
Output created: `merged_rsem` directory with by cbio disease `rsem_merged.txt` and `rsem_merged_zscore.txt` files.

### scripts/7_convert_PBTA_fusion.py
*Prerequisite: openPBTA fusion file from releases [here](https://s3.console.aws.amazon.com/s3/buckets/kf-openaccess-us-east-1-prd-pbta/data/?region=us-east-1&tab=overview)* Sort of a hack until we make annoFuse part of the workflow, to incorporate merged and filtered fusions from annoFuse created for the openPBTA project
```usage: 7_convert_PBTA_fusion.py [-h] [-t TABLE] [-f FUSION_FILE] [-s SQ_FILE]
                                [-j CONFIG_FILE]

Convert openPBTA fusion table to cbio format.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file
                        names
  -f FUSION_FILE, --fusion-file FUSION_FILE
                        openPBTA fusion file
  -s SQ_FILE, --center-file SQ_FILE
                        File with BS IDs and sequencing centers
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
```
Data processing config file is used, with relevant fields:
+ `entrez_tsv`: tsv file that is basically the revers of th ehugo tsv file.  Should unfiy that...
Example format:
```
ENTREZ_GENE_ID	HUGO_GENE_SYMBOL
100506173	1060P11.3
6775087	12S RRNA
8923213	12S RRNA
8923219	16S RRNA
109951028	A-GAMMA3'E
1	A1BG
503538	A1BG-AS1
```
Outputs: `merged_fusion` directory with merged fusion files by cbio disease as outlined [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#fusion-data). Files have extension `fusions.txt`