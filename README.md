# Outline on ETL for converting data from cavatica and data service to pedcbioportal format
In general, we are creating upload packages converting our data and metadata to satisfy the requirements outlined [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats).
Further general loading notes can be found in this [Notion page](https://www.notion.so/d3b/Cbioportal-Study-Load-SOP-58812479fabe4d2fa9f72242e331b5ee).
See [below](#collaborative-and-publication-workflows) for special cases like publications or collaborative efforts
## Software Prerequisites

+ `python3` v3.5.3+
  + `sevenbridges-python` package (https://sevenbridges-python.readthedocs.io/en/latest/installation/)
  + `numpy`, `pandas`, `scipy`
+ `bedtools` (https://bedtools.readthedocs.io/en/latest/content/installation.html)
+ `chopaws` https://github.research.chop.edu/devops/aws-auth-cli needed for saml key generation for s3 upload
+ access to https://aws-infra-jenkins-service.kf-strides.org to start cbio load into QA and/or prod using the `d3b-center-aws-infra-pedcbioportal-import` task
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
```python
usage: 1a_merge_cavatica_manifests.py [-h] [-m MANIFEST]

Get all relevant analyzed file outputs from cavatica manifest.

optional arguments:
  -h, --help            show this help message and exit
  -m MANIFEST, --manifest MANIFEST
                        cavatica manifest file csv list
```

### scripts/1b_query_file_manifest.py
This will create a master sheet bs ids, analyte type, project location, and file names.
```python
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
```python
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
```python
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
```python
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

## Helper/Patch scripts
Additional metadata may need to be patched as they may not have existed before as a standard. Skip this section if you have all you need
### utilities/patch_sample_id.py
Used in the `pbta_all` load in April 2021, it takes data exported from the D3b data warehouse to creat valid sample IDs
```python
usage: patch_sample_id.py [-h] [-i INFO] [-d DATA]

Quickly patch sample id from step2 with data warehouse entries

optional arguments:
  -h, --help            show this help message and exit
  -i INFO, --info INFO  A sample info sheet from step 2
  -d DATA, --data DATA  Data Warehouse tsv with biospecimen ID, master aliquot
                        id, and sample id
```
Writes to stdout. Output has the sample IDs updated from the sheet inputted from step2. Run for each sheet.
### utilities/get_survival_from_ds.py
This is the most likely helper script needed. This gets survival data from the data service to be patched to the data_clinical_patient.txt sheets.
```python
usage: get_survival_from_ds.py [-h] [-u URL] [-p PT_LIST]

Get vital status by PT ID

optional arguments:
  -h, --help            show this help message and exit
  -u URL, --kf-url URL  Kids First data service url, i.e. https://kf-api-
                        dataservice.kidsfirstdrc.org/
  -p PT_LIST, --pt-list PT_LIST
                        pt_id list
```
Output: Headerless tsv file called `outcomes.txt`. First column PT ID, second column vitality status, third column age in days at last status update

### utilities/patch_survival_to_data_patient.py
Adds survival data obtained from `outcomes.txt` and calculates overall survival in months using age at diagnosis and age at status with formula 

## Prepare mutation data
### scripts/4_merge_maf.py
*Prerequisite: Download all maf files to be merged into a directory, using the manifest and wget, sbg api, or xargs and sb cli.* Merges maf files by cbio disease type.  Will remove calls of certain categories as outlined [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#mutation-data)
```python
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

## Prepare Copy number data (if avaliable)
### scripts/5a_cnv_genome2gene.py
*Prerequisite: Download all ControlFreeC p value output files to be converted into a directory, using the manifest and wget, sbg api, or xargs and sb cli.* Converts copy number genome coordinate calls to gene copy number calls
```python
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
*Prerequisite: A cavatica file manifest of ControlFreeC `info.txt` files* Merge files generated from step 5a into a genes-by-sample raw copy number matrix
```python
usage: 5b_merge_cnv.py [-h] [-t TABLE] [-n CNV_DIR] [-m MANIFEST]
                       [-j CONFIG_FILE] [-p PROFILE]

Merge cnv files using cavatica task info and datasheets.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file
                        names
  -n CNV_DIR, --cnv-dir-gene CNV_DIR
                        cnv as gene file directory
  -m MANIFEST, --info_manifest MANIFEST
                        cavatica cfree info file manifest
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
  -p PROFILE, --profile PROFILE
                        cavatica profile name. requires
                        .sevenbridges/credentials file be present
```
Data processing config file is used, with relevant fields:
+ [`dna_ext_list`][`copy_number`]: File extension of inputs (used in 5a) with no leading `.`, i.e. `controlfreec.CNVs.p.value.txt`
+ `cpus`: Number of files to process concurrently using `ProcessPoolExecutor`.  8 recommended
Output created: `merged_cnvs` directory with merged raw copy number files named based on cbio disease type. Will have extension `predicted_cnv.txt`

### scripts/5c_cnv_discrete.py
*Prerequisite: A cavatica file manifest of ControlFreeC `info.txt` files* Converts raw copy number to gistic-style discrete values. Uses info.txt files to adjust for calculated ploidy
```python
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

### scripts/5d_merge_seg.py
Simply merges seg files as-is
```python
usage: 5d_merge_seg.py [-h] [-t TABLE] [-i HEADER] [-m SEG_DIR]
                       [-j CONFIG_FILE]

Merge seg files

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file
                        names
  -i HEADER, --header HEADER
                        File with seg header only
  -m SEG_DIR, --seg-dir SEG_DIR
                        seg file directory
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
```
Output created: In `merged_cnvs` directory with merged discrete copy number files named based on cbio disease type. Will have extension `merged_seg.txt`

## Prepare RNA expression and fusion data (if avaliable)
### scripts/6_merge_rename_rsem.py
*Prerequisite: Download all rsem files to be merged into a directory, using the manifest and wget, sbg api, or xargs and sb cli.* Merges rsem files into a gene-by-sample matrix.  Repeat HUGO symbols are resolved by taking the row with the highest mean expression and dropping the rest. Also creates a merged z score table.  Z scores are calculated before breaking up into individual cbio disease files; FPKM counts have a pseudocount added, then are log2 transformed before z score is calculated
```python
usage: 6_merge_rename_rsem.py [-h] [-t TABLE] [-r RSEM_DIR] [-j CONFIG_FILE]

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

### scripts/7_convert_fusion.py
*Prerequisite: openPBTA fusion file from releases [here](https://s3.console.aws.amazon.com/s3/buckets/kf-openaccess-us-east-1-prd-pbta/data/?region=us-east-1&tab=overview) OR annoFuse results downloaded into a directory.* Can use openPBTA results or annoFuse results for cbio fusion file table generation.

```python
usage: 7_convert_fusion.py [-h] [-t TABLE] [-f FUSION_RESULTS] [-m MODE]
                           [-s SQ_FILE] [-j CONFIG_FILE]

Convert openPBTA fusion table OR list of annofuse files to cbio format.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file
                        names
  -f FUSION_RESULTS, --fusion-results FUSION_RESULTS
                        openPBTA fusion file OR annoFuse results dir
  -m MODE, --mode MODE  describe source, pbta or annofuse
  -s SQ_FILE, --center-file SQ_FILE
                        File with BS IDs and sequencing centers. Should have
                        headered columns: BS_ID\SQ_Value
```

Outputs: `merged_fusion` directory with merged fusion files by cbio disease as outlined [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#fusion-data). Files have extension `fusions.txt`. Entrez ID is left blank to allows for reversible fusion gene searching

## Create upload package
### scripts/8_organize_upload_packages.py
Create cases lists, meta files, and organize data for cbio upload. It is
assumed you are at the dir level of all input data files
```python
usage: 8_organize_upload_packages.py [-h] [-o OUT_DIR] [-d DX_FILE]
                                     [-c CONFIG_FILE]

Create cases lists, meta files, and organize data for cbio upload. It is
assumed you are at the dir level of all input data files

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_DIR, --output_dir OUT_DIR
                        output directory name
  -d DX_FILE, --dx-file DX_FILE
                        Dx index file with cbio short and long names
  -c CONFIG_FILE, --config CONFIG_FILE
                        json config file with meta information; see
                        REFS/case_meta_config.json example
```
Case configuration file is used.  It is best to start with one of the pre-configured case config files in REFS/case_{descriptive text}. 
### Single study load
A good example of a *single study load* is in `REFS/case_pbta_combined_meta_config.json`. Likely personalized edits would occur in the following fields:
+ `merged_{data type}`: The `profile_description` key in each is a good place to describe any algorthim or nuances used to generate the data of that type
+ `study`: Here is where you set the overall study description, it's the banner text that people will see in the study overview page that gives them a sense of what the data is.
  + `description`: This field is set up as an array so that a generic form of "text describing" "disease" "more text describing" can be used. Put another way, element one is whatever you want to say about the disease/study until you are ready to mention the disease/study, element two anythnig you may optionally wish to add
  + `groups`: These are access groups defined is cBioportal.  Default is `PUBLIC`, but another can be named is restrictions are needed.  Need to work with Devops for custom groups
  + `cancer_study_identifier`: This is the short name that you create for the study.  It will be the name of the study load folder and will be used by cBioportal to find all relevant information for that study.  *Be sure not to name it the same as an existing study that you do not wish to edit/replace, or it will be overwritten!*
  + `dir_suffix`: Leave blank for single study.  Only relevant for multi-study mode
  + `name_append`: Additional descriptice text to add to study name display.  Study name display is obtained from the `dx_index.txt` file, from the `MB cbioportal name` column

### Multi-study load
A good example of a *multi-study load* is in `REFS/case_pbta_by_dx_meta_config.json`. This is when a large study has multiple studies/diseases within it, that ought to be broken up into seprate studies/disease as determined in the `dx_index.txt` file. Likely personalized edits would occur in the following fields:
+ `merged_{data type}`: The `profile_description` key in each is a good place to describe any algorthim or nuances used to generate the data of that type
+ `study`: Here is where you set the overall study description, it's the banner text that people will see in the study overview page that gives them a sense of what the data is.
  + `description`: This field is set up as an array so that a generic form of "text describing" "disease" "more text describing" can be used. Put another way, element one is whatever you want to say about the disease/study until you are ready to mention the disease/study, element two anythnig you may optionally wish to add
  + `groups`: These are access groups defined is cBioportal.  Default is `PUBLIC`, but another can be named is restrictions are needed.  Need to work with Devops for custom groups
  + `cancer_study_identifier`: Leave this field blank! The `MB cbioportal id` field from the `dx_index.txt` sheet will be used to set the short study name, woth the next field described below used to help differntaite it from other studies of the same disease type.
  + `dir_suffix`: Add a suffix to the short study name to make it less generic and moreeasily identifiable.  For instance, `plgg` from cbttc has `_cbttc` set for this field, so that the final short study name becomes `plgg_cbttc`
  + `name_append`: Additional descriptice text to add to study name display.  Study name display is obtained from the `dx_index.txt` file, from the `MB cbioportal name` column

Outputs: Study directories created a subdirectories in the scripyr output directory.  These subdirectories named as the study short names are to be uploaded.
## Upload the final packages
Upload all of the directories named as study short names to `s3://kf-cbioportal-studies/public/`. You may need to set and/or copy aws your saml key before uploading. Next, edit the file in that bucket called `importStudies.txt` located at `s3://kf-cbioportal-studies/public/importStudies.txt`, with the names of all of the studies you wish to updated/upload. Lastly, go to https://jenkins.kids-first.io/job/d3b-center-aws-infra-pedcbioportal-import/job/master/, click on build. At the `Promotion kf-aws-infra-pedcbioportal-import-asg to QA` and `Promotion kf-aws-infra-pedcbioportal-import-asg to PRD`, the process will pause, click on the box below it to affirm that you want these changes deployed to QA and/or PROD respectively.  If both, you will have to wait for the QA job to finish first before you get the prompt for PROD.
## Congratulations, you did it!

# Collaborative and Publication Workflows
These are highly specialized cases in which all mor most of data comes froma thrid party, and therefore requires specific case-by-case protocols

## OpenTargets
This project is organized much like OpenPBTA in which all genomics data for each assay-type are collated into one giant table.
In general, this fits cBioPortal well.
See `s3://kf-cbioportal-studies/public/ped_opentargets_2021/` for final load packages, `https://pedcbioportal.kidsfirstdrc.org/study/summary?id=ped_opentargets_2021` for final product.

### Inputs
Inputs are located in the old Kids First AWS account (`538745987955`) in this general bucket location: `s3://kf-openaccess-us-east-1-prd-pbta/open-targets/pedcbio/`.
Clinical data are obtained from the `histologies.tsv` file.
Genomic data generally obtained as such:
 - Somatic variant calls: merged maf
 - Copy number: tsv file with copy number, ploidy, and GISTIC-style information in maf-like format (each call is a row)
 - RNA expression: tpm values from rsem stored an `.rds` object
 - RNA fusion: annoFuse output

 ### File Transformation

#### 1. COLLABORATIONS/openTARGETS/clinical_to_datasheets.py
 ```
usage: clinical_to_datasheets.py [-h] [-f HEAD] [-c CLIN] [-s CL_SUPP]

Script to convert clinical data to cbio clinical data sheets

optional arguments:
  -h, --help            show this help message and exit
  -f HEAD, --header-file HEAD
                        tsv file with input file original sample names, output
                        sheet flag, and conversion
  -c CLIN, --clinical-data CLIN
                        Input clinical data sheet
  -s CL_SUPP, --cell-line-supplement CL_SUPP
                        supplemental file with cell line meta data - bs
                        id<tab>type. optional
 ```
 - `-f` Header file example can be found here: `COLLABORATIONS/openTARGETS/header_desc.tsv`
 - `-s` cell line supplemental file here: `REFS/cell_line_supplemental.txt`
 - `-c` `histologies.tsv`

 Outputs a `data_clinical_sample.txt` and `data_clinical_patient.txt` for the cBio package, and a `bs_id_sample_map.txt` mapping file to link BS IDs to gnerated cBioPortal IDs based on the rules for creating a proper somatic event using column `parent_aliquot_id`

_NOTE:_ histologies.tsv was subset on `PBTA` for the initial run

Example run:
`python3 COLLABORATIONS/openTARGETS/clinical_to_datasheets.py -f COLLABORATIONS/openTARGETS/header_desc.tsv -s REFS/cell_line_supplemental.txt -c histologies.tsv`

#### 2. COLLABORATIONS/openTARGETS/rename_filter_maf.py
Rename IDs in `Tumor_Sample_Barcode`
```
usage: rename_filter_maf.py [-h] [-m MAPPING_FILE] [-v MAF_FILE]

Script to pre-filter entries on usually removed criteria except TERT
promoter, convert BS IDs to cBio names

optional arguments:
  -h, --help            show this help message and exit
  -m MAPPING_FILE, --mapping-file MAPPING_FILE
                        tsv file with header and bs_id, sample type, cbio ID mappings
  -v MAF_FILE, --maf-file MAF_FILE
                        openX maf file
```

Example run:
`python3 COLLABORATIONS/openTARGETS/rename_filter_maf.py -m bs_id_sample_map.txt -v snv-consensus-plus-hotspots.maf.tsv.gz`

#### 3. COLLABORATIONS/openTARGETS/cnv_to_tables.py
Convert cnv table to cBio format - genes as rows, samples as cols, one for absolute CN, another for GISTIC-style
```
usage: cnv_to_tables.py [-h] [-m MAPPING_FILE] [-c CNV_TBL]

Script to convert openX cnv table to cBio format

optional arguments:
  -h, --help            show this help message and exit
  -m MAPPING_FILE, --mapping-file MAPPING_FILE
                        tsv file with header and bs_id, sample type, cbio ID mappings
  -c CNV_TBL, --copy-number CNV_TBL
                        openX table
```

Example run:
`python3 COLLABORATIONS/openTARGETS/cnv_to_tables.py -m bs_id_sample_map.txt  -c consensus_seg_annotated_cn_autosomes_xy.tsv.gz`

#### 4. COLLABORATIONS/openTARGETS/rename_export_rsem.R
This script is more of an outline of how exactly it was run rather than quite able to be run as standalone.
A future update will include using the existing mapping file and being able to take args

#### 5. COLLABORATIONS/openTARGETS/case_list_from_datasheet.py
```
usage: case_list_from_datasheet.py [-h] [-d DATASHEET] [-s STUDY_ID]

Generate extra cases lists based on cBio datasheet

optional arguments:
  -h, --help            show this help message and exit
  -d DATASHEET, --datasheet DATASHEET
                        cBio sample datasheet
  -s STUDY_ID, --study-id STUDY_ID
                        cBio cancer_study_identifier
```

Example run:
`python3 COLLABORATIONS/openTARGETS/case_list_from_datasheet.py -d data_clinical_sample.txt -s ped_opentargets_2021`

### TODO
 - Improve step 4
 - Automate or store meta data files
 - Automate standard case list generation