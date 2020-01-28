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


