# Outline on ETL for converting data from cavatica and data service to pedcbioportal format
## Prerequisites
+ `R`
  + `rjson` package
+ `python3`
  + `sevenbridges` package (https://sevenbridges-python.readthedocs.io/en/latest/installation/)
+ `Rscript`
+ `bedtools` (https://bedtools.readthedocs.io/en/latest/content/installation.html)

### Quick Start Example
You're a BAMF with no time for READMEs.  It's recommended to use the wrapper script as it will check through configs and inputs to ensure as much as possible is in place.
 Below is a skeleton example of how a project from start start to finish could go down:

```bash
# all run within same dir, scripts in datahub-cbttc/scripts/ dir of repo
./0_ETL_wrapper.py -s 1,2,3 -c cbttc_complete_tasks.txt -p cavatica -u https://kf-api-dataservice.kidsfirstdrc.org -j cbttc_complete_config.json [-s cell_line_metadata.txt] 2> step0.errs > step0.out &
./0_ETL_wrapper.py -s 4,5a,5b,6,7 -c cbttc_complete_tasks.txt -j cbttc_complete_config.json -p cavatica -i <repo_loc>/datahub-cbttc/scripts/REFS/maf_header_std.txt -d ./CNV -r ./rsem_raw -m ../mafs 2> step0.errs > step0.out &

```
Sample config file:
```json
{
  "cav_dna_projects":["kfdrc-harmonization/sd-bhjxbdqk-07", "kfdrc-harmonization/sd-bhjxbdqk-08"],
  "cav_rna_projects":["kfdrc-harmonization/sd-bhjxbdqk-06", "zhangb1/cbttc-disease-expression"],
  "tn_order": "TN",
  "data_dir": "/home/ubuntu/OpenDIPG-PNOC_DIPG_WGS",
  "bedtools": "bedtools",
  "cp_only_script": "/home/ubuntu/datahub-cbttc/scripts/get_cbio_copy_only_num.pl",
  "bed_genes": "/home/ubuntu/datahub-cbttc/scripts/REFS/GRCh38.84.gtf_genes.bed",
  "hugo_tsv": "/home/ubuntu/datahub-cbttc/scripts/REFS/HUGO_EntrezID.tsb",
  "dna_samp_id_style": "cbttc_dna_std",
  "rna_samp_id_style": "cbttc_rna_std",
  "ens_gene_list":"ensemble_gene_list.txt",
  "script_dir": "/home/ubuntu/datahub-cbttc/scripts/",
  "group": "CBTTC_TEST",
  "cancStudyID": "OpenDIPG-PNOC_DIPG_WGS",
  "dx_tbl_fn": "dipg_dx.txt",
  "mut_desc": "Mutation data derived from WGS data from Nantomics (https://nantomics.com/).",
  "study_short_desc": "(Open DIPG, CBTTC) WGS",
  "study_desc_p1": "Genomic Characterization of",
  "study_desc_p2": "samples provided by the <a href=\"http://CBTTC.org\">Children's Brain Tumor Tissue Consortium</a> and its partners via the <a href=\"http://kidsfirstdrc.org\">Gabriella Miller Kids First Data Resource Center</a>. The Open DIPG Initiative has collected, generated and annotated the largest cohort of DIPG genome data to date. It is part of a larger effort known as the Pediatric Brain Tumor Atlas, which aims to uncover the molecular basis of childhood cancers.",
  "fpkmFileLoc": "/home/ubuntu/cbttc_RNA_res/expr/rsem_merged.txt",
  "mapFileLoc": "/home/ubuntu/cbttc_RNA_res/expr/mappingFile.txt",
  "genesFileLoc": "/home/ubuntu/tools/datahub-cbttc/scripts/REFS/genes.tsv",
  "prepend_cancStudyID_flag": 0,
  "rna_flag": 0,
  "cna_flag": 1,
  "threads": 40,
  "cpus": 8,
  "ds_samp_desc": ["#Patient Identifier\tSample Identifier\tSPECIMEN_ID\tCANCER_TYPE\tCANCER_TYPE_DETAILED\tTUMOR_TISSUE_SITE\tTUMOR_TYPE\tSAMPLE_TYPE\tMATCHED_NORMAL_SAMPLE_ID\tMATCHED_NORMAL_SPECIMEN_ID\tCBTTC_PAIRED_IDS",
    "#Patient identifier\tSample Identifier using external_sample_id\tkfdrc tumor biopsecimen ID\tStudy-defined cancer type\tStudy-defined cancer type detail\ttumor tissue location\tprimary v metastatic tumor designation\tpatient tissue sample or cell line\tmatched normal external_sample_id\tkfdrc matched normal biospecimen ID\tOriginal CBTTC DNA/RNA pair ids, where applicable",
    "#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING",
    "#1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1",
    "PATIENT_ID\tSAMPLE_ID\tSPECIMEN_ID\tCANCER_TYPE\tCANCER_TYPE_DETAILED\tTUMOR_TISSUE_SITE\tTUMOR_TYPE\tSAMPLE_TYPE\tMATCHED_NORMAL_SAMPLE_ID\tMATCHED_NORMAL_SPECIMEN_ID\tCBTTC_PAIRED_IDS\n"],
  "ds_pt_desc": ["#Patient Identifier\tExternal Patient Identifier\tGENDER\tAGE\tTUMOR_SITE\tRACE\tETHNICITY",
    "#Patient identifier\tPatient ID used by generator of data\tGender or sex of the patient\tAge at which the condition or disease was first diagnosed, in years\tTumor location\tracial demographic\tethnic demographic",
    "#STRING\tSTRING\tSTRING\tNUMBER\tSTRING\tSTRING\tSTRING",
    "#1\t0\t1\t1\t1\t1\t1",
    "PATIENT_ID\tEXTERNAL_PATIENT_ID\tGENDER\tAGE\tTUMOR_SITE\tRACE\tETHNICITY\n"]
}
```

### Step-by-step inputs and outputs
Step 0)
If you wish to run the wrapper script, the help gives a decent amount of info as to what's needed for what part of the support and final product needed for file load
```bash
./0_ETL_wrapper.py -h
usage: 0_ETL_wrapper.py [-h] [-s STEPS] [-c CAV] [-p PROFILE] [-j CONFIG_FILE]
                        [-u URL] [-t TUMOR] [-n NORMAL] [-i HEADER]
                        [-m MAF_DIR] [-d CNV_DIR] [-e CNV_DIR_GENE]
                        [-r RSEM_DIR]

Wrapper script to run all or some of ETL component scripts.

optional arguments:
  -h, --help            show this help message and exit
  -s STEPS, --steps STEPS
                        csv list of steps to execute, valid choices are: "all", "1", "2", "3", "4", "5a", "5b", "6", "7".  Refers to scripts:
                        1_query_cavatica.py
                        2_query_ds_by_bs_id.py
                        3_create_data_sheets.py
                        4_merge_maf.py
                        5a_cnv_genome2gene.py
                        5b_merge_cnv.py
                        6_merge_rsem.py
                        7_createLoadingFiles.R
  -c CAV, --cavatica CAV
                        cavatica output file name, relevant steps: 1, 2, 3, 4, 5b, 6
  -p PROFILE, --profile PROFILE
                        cavatica profile name. requires .sevenbridges/credentials file be present. Relevant steps: 1,6
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
  -u URL, --kf-url URL  Kids First data service url, i.e. https://kf-api-dataservice.kidsfirstdrc.org/. Relevant step: 2
  -t TUMOR, --tumor TUMOR
                        file with tumor metadata (see step 2).  Relevant step: 3, not needed if step 2 run
  -n NORMAL, --normal NORMAL
                        file with normal metadata (see step 2) Relevant step: 3, not needed if step 2 run
  -i HEADER, --header HEADER
                        File with maf header only. Relevant step: 4
  -m MAF_DIR, --maf-dir MAF_DIR
                        maf file directory. Relevant step: 4
  -d CNV_DIR, --cnv-dir CNV_DIR
                        cnv as genome coords file directory. Relevant step: 5a
  -e CNV_DIR_GENE, --cnv-dir-gene CNV_DIR_GENE
                        cnv as gene file directory.Relevant step: 5b, not needed if step5a run
  -r RSEM_DIR, --rsem-dir RSEM_DIR
                        rsem file directory. Relevant step: 6
```
Step 1)
```bash
./1_query_cavatica.py -h
usage: 1_query_cavatica.py [-h] [-o OUT] [-p PROFILE] [-c CONFIG_FILE]

Get all relevant analyzed file outputs from projects in cavatica.

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --output OUT  output file name
  -p PROFILE, --profile PROFILE
                        cavatica profile name. requires
                        .sevenbridges/credentials file be present
  -c CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
```
Relevant config file entries:

  + `cav_dna_projects` array of username/project ids from cavactica containing your DNA data
  + `cav_rna_projects` array of username/project ids from cavactica containing your RNA data.  If empty, no RNA data will be search for
  + `tn_order` tumor-normal order for bs ids in cavatica task name.  TN or NT
  + `cna_flag` `0` for no cnv data, `1` if there is cnv data.  For instance, a WGS project probably would have some, a WES probably not
  + `cpus` integer number of cpus to split task walking in cavatica
Outputs:
  + stderr
  + Output tsv file set by input param `-o`.  Contains a header like this, with related info: 
  `T/CL BS ID      Norm BS ID      Task Name       Task ID Analyte Type    Relevant Outputs        Source Project`

Step 2)
```bash
./2_query_ds_by_bs_id.py -h
usage: 2_query_ds_by_bs_id.py [-h] [-u URL] [-c CAV]

Script to walk through data service and grab all relevant meetadata by bs id.

optional arguments:
  -h, --help            show this help message and exit
  -u URL, --kf-url URL  Kids First data service url, i.e. https://kf-api-
                        dataservice.kidsfirstdrc.org/
  -c CAV, --cavatica CAV
                        file with task info from cavatica (see step 1)

```
For now, this is the only step that can't be run on a vm.
Relevant config file entries:
  + `threads` int number of threads to query the data service with.  40 is pretty good.
  
Outputs:
  + stderr
  + two tsv files with the header:
  `BS_ID   PT_ID   external_aliquot_id     analyte_type    source_text_tissue_type source_text_tumor_descriptor    composition     external_sample_id      external_id     gender  ethnicity       race    source_text_diagnosis   source_text_tumor_location      age_at_event_days`
    + `tum_bs_ds_info.txt` Tumor/Cell Line sample related info for both DNA and RNA
    + `norm_bs_ds_info.txt` Normal sample related info for both DNA

Step 3)
```bash
./3_create_data_sheets.py -h
usage: 3_create_data_sheets.py [-h] [-c CAV] [-t TUMOR] [-n NORMAL]
                               [-j CONFIG_FILE]

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
                        id<tab>type. OPTIONAL
``` 
Relevant config file entries:

  + `dna_samp_id_style` Name of style to use describing how the DNA sample names will be constructed from data service
  + `rna_samp_id_style` Name of style to use describing how the RNA sample names will be constructed from data service. Can be blank if no RNA to be loaded.
  + `ds_pt_desc` starting 5 lines of `data_clinical_patient.txt` sheets to be used for project.  Shared by all, fourth line determines what will be shown to maintain relevance
  + `ds_samp_desc` as above for `data_clinical_sample.txt`
  + `dx_tbl_fn` Name of file containing dictionary for converting data ervice diagnoses to cbio disease types Example content:
    ```text
    CBTTC Diagnosis MB cbioportal id        MB cbioportal name
    Adenoma aptad   Atypical Pituitary Adenoma
    Atypical Teratoid Rhabdoid Tumor (ATRT) atrt    Atypical Teratoid/Rhabdoid Tumor
    Brainstem glioma- Diffuse intrinsic pontine glioma      dipg    Diffuse Intrinsic Pontine Glioma
    Chordoma        chdm    Chordoma
    Choroid plexus carcinoma        cpc     Choroid Plexus Carcinoma
    Choroid plexus papilloma        cpp     Choroid Plexus Papilloma
    Craniopharyngioma       cranio  Craniopharyngioma
    Dysembryoplastic neuroepithelial tumor (DNET)   dnt     Dysembryoplastic Neuroepithelial Tumor
    Dysplasia/Gliosis       other   other
    ```
   + `rna_flag` `0` for no RNA data, `1` if RNA data is present.  
Outputs:
  + stderr - IMPORTANT: This file will contain warnings of which DNA and RNA samples did not meet naming conventions listed in the Description of scripts section.
  + `datasheets` directory created in cwd.  Has subdirs for each portal project.
  + `mappingFile.txt` If rna data are present, this file will be created for later use to tie rsem output names to sample names to be used in project

Step 4)
```bash
./4_merge_maf.py -h
usage: 4_merge_maf.py [-h] [-c CAV] [-i HEADER] [-m MAF_DIR]

Merge and filter mafs using cavatica task info and datasheets.

optional arguments:
  -h, --help            show this help message and exit
  -c CAV, --cavatica CAV
                        file with task info from cavatica (see step 1)
  -i HEADER, --header HEADER
                        File with maf header only
  -m MAF_DIR, --maf-dir MAF_DIR
                        maf file directory

```
Relevant config file entries:
  + `cpus`
  
Outputs:
  + stderr
  + output dir in cwd called `merged_mafs/`
  
Step 5a)
```bash
./5a_cnv_genome2gene.py -h
usage: 5a_cnv_genome2gene.py [-h] [-d CNV_DIR] [-j CONFIG_FILE]

Convert controlFreeC cnv genome coords to gene level

optional arguments:
  -h, --help            show this help message and exit
  -d CNV_DIR, --cnv-dir CNV_DIR
                        cnv as genome coords file directory
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
```  
Relevant config file entries:
  + `cpus`
  
Outputs:
  + stderr
  + output dir in cwd called `converted_cnvs/`

Step 5b)
```bash
./5b_merge_cnv.py -h
usage: 5b_merge_cnv.py [-h] [-c CAV] [-n CNV_DIR]

Merge cnv files using cavatica task info and datasheets.

optional arguments:
  -h, --help            show this help message and exit
  -c CAV, --cavatica CAV
                        file with task info from cavatica (see step 1)
  -n CNV_DIR, --cnv-gene-dir CNV_DIR
                        cnv as gene file directory
```
Relevant config file entries:
  + `cpus`

Outputs:
  + stderr
  + output dir in cwd called `merged_cnvs/`
 
Step 6)
```bash
./6_merge_rsem.py -h
usage: 6_merge_rsem.py [-h] [-c CAV] [-r RSEM_DIR] [-j CONFIG_FILE]

Merge rsem files using cavatica task info.

optional arguments:
  -h, --help            show this help message and exit
  -c CAV, --cavatica CAV
                        file with task info from cavatica (see step 1)
  -r RSEM_DIR, --rsem-dir RSEM_DIR
                        rsem file directory
  -j CONFIG_FILE, --config CONFIG_FILE
                        json config file with data types and data locations
  -p PROFILE, --profile PROFILE
                        cavatica profile name. requires .sevenbridges/credentials file be present
```
Relevant config file entries:
  + `ens_gene_list` ENSEMBL gene ID gene symbol tsv list.  i.e.:
    ```text
    gene_id gene_symbol
    ENSG00000000003.14      TSPAN6
    ENSG00000000005.5       TNMD
    ENSG00000000419.12      DPM1
    ENSG00000000457.13      SCYL3
    ENSG00000000460.16      C1orf112
    ENSG00000000938.12      FGR
    ENSG00000000971.15      CFH
    ENSG00000001036.13      FUCA2
    ENSG00000001084.10      GCLC
    ```
  + `cpus`
  
Outputs:
  + stderr
  + output file in subdir created in cwd called `merged_rsem/rsem_merged.txt`

Step 7)
```bash
Rscript 7_createLoadingFiles.R <config file>
```
Relevant config file entries:
  + `data_dir` Location where you have been running steps 1-7
  + `rna_flag` 0 if no RNA to be loaded, 1 if RNA to be loaded
  + `cna_flag` 0 if no CNV data to be loaded, 1 if CNV data to be loaded
  + `dx_tbl_fn` Name of file containing dictionary for converting data ervice diagnoses to cbio disease types. See Step 3 description above
  + `mut_desc` Description of where DNA data came from relating to the sequencing core or company used to generate data
  + `cancStudyID` Short no-spaces name describing main group(s), hospital, foundation, etc that generated the data.  Generally in format `{source}-{disease_type}`
  + `prepend_cancStudyID_flag` `0` if single disease and `cancStudyID` will be used as the single descriptor/directory for output in `processed` folder, `1` if cbio short disease name to preprend `cancStudyID` and each will be used as descriptor/directory for output in `processed` folder
  + `group` authorization group.  Use `PUBLIC` if intended to be displayed to ALL, `CBTTC_TEST` for example to show only to D3b members.
  + `fpkmFileLoc` location of output file from Step 6
  + `mapFileLoc` location of `mappingFile.txt` from Step 3
  + `genesFileLoc` Location of `genes.tsv` file - can be found in `REFS` of this repo.  Example format:   
    ```text
    ENTREZ_GENE_ID  HUGO_GENE_SYMBOL        GENETIC_ENTITY_ID       TYPE    CYTOBAND        LENGTH
    -1977   PRKCA_PS664     45407   phosphoprotein  17q22-q23.2     NULL
    -1976   MAPK8_PT183_Y185        35453   phosphoprotein  10q11.22        NULL
    -1975   WWTR1_PS89      60392   phosphoprotein  3q23-q24        NULL
    -1974   ERBB3_PY1298    8317    phosphoprotein  12q13   NULL
    -1973   EGFR_PY992      7876    phosphoprotein  7p12    NULL
    -1972   EIF4EBP1_PT37   8025    phosphoprotein  8p12    NULL
    -1971   RET_PY905       46695   phosphoprotein  10q11.2 NULL
    -1970   GYS1_PS641      11365   phosphoprotein  19q13.3 NULL
    -1969   PRRT2_PS660     45597   phosphoprotein  16p11.2 NULL
    ```
   + `study_short_desc` A suffix to be added to the long-form disease description.
   + `study_desc_p1` Long form study description to prepend to the disease name
   + `study_desc_p2` Long form study description to append to the disea name.
Outputs:
  + stdout
  + stderr
  + `processed` directory.  Copy it's contents recursively to the `public` directory in the repo and push as instructed in the beginnign to run validation and pull request to load a study.  If updating an existing study, may need to remove directories to be updated first

## Helpful tips
  + *IMPORTANT*: This entire api has been designed to run all in the working directory.  All output files and subdirs are intended to be created in the same directory in which you run each and every script from.
  + If you want to create a project for only s subset of samples, filter out the cavatica task file from Step 1, tumor and normal tables from Step 2, and run the rest of the steps on those inputs
  + Common sense - if the relevants inputs exists already, no need to run that step again!
  + cBio Display Behaviors
    - For display name, the second column of the `dx_tbl_fn` file and the `study_short_desc` variable combined are what set that, i.e. `Pediatric High Grade Glioma (CBTTC)`
    - If you want a subdir for each cbio short name for loading and display, set `prepend_cancStudyID_flag` to `1` i.e. `Pediatric High Grade Glioma (CBTTC)`, `processed` subdir `phgg_cbttc`
    - If a study should be "all one disease", set all possible `dx_tbl_fn` defs for each data service diagnosis to the desired cbio disease and set `prepend_cancStudyID_flag` to `0` i.e. `Diffuse Intrinsic Pontine Glioma (Open DIPG, PNOC) clinical sequencing`, `processed` subdir `OpenDIPG-PNOC_DIPG_WES`

### Description of scripts
This summarizes what each script is doing in this repo

0) A wrapper script that will run either `all` (not quite possible as data service is not accessible by vm, and a local machine might not have the memory/storage for processing larger datasets) or a csv list of steps, represented by the alphanumeric start of each script below to run.
1) All DNA and RNA somatic pipeline runs are collated extracting task IDs, relevant file names, and associated BS IDs outputting to a table
2) Query data service for each biospcimen ID getting all related diagnoses, patient info including KF PT IDs, sample IDs (external sample id) assigned during the project and output a table with tumor sample nifo and a separate table for normals
3) Convert outputs from step 2 into data_clinical_patient.txt and data_clinical_sample.txt sheets fit for cbioportal and placing each in a disease-specific directory
    - A diagnosis dictionary is used to convert CBTTC-specific diagnoses into something that exists in cbioportal
    - Tumor(or CL)/normal pairing is gleaned from the task names from the cavatica jobs by BS IDs
    - Original kf biospecimen IDs are `;` separated, and located in the `SPECIMEN_ID` column
    - `external_sample_id` is used to tie DNA samples to their RNA sample counterparts.  If the `composition` is `Derived Cell Line`, the sample id is augmented with a `-CL` suffix. Caveats:
        - `sample_id_builder_helper.py` script has defs so that other project sample names can be defined and derived.  If a rule doesn't exist in that script that doesn't fit your needs, create a new one and update the config file to use that rule.
        - If two samples of the same analyte (DNA, RNA) type have the same sample ID, they are excluded. Tumor samples mapped against multiple normals are also excluded
        - Matched DNA and RNA *should* have the same sample IDs.  However, some legacy samples only line up by their 7316-NNN prefixes.  In those instances, those first two parts are used to link the samples and become the new sample ID.  The original sample IDs are recorded in the last column of the data sheet in the `CBTTC_PAIRED_IDS` column, so that it can be traced back if needed
        - For patient earliest age of diagnosis, the minimum age associated with that diagnosis for the sample, is used from the data service
        - For cell lines, if you have a file that can help describe them, supply that file to the optional `-s` arg.  Format of the file should be tsv, headerless, `bs_id<tab>desc` `desc` should be abbreviated and have no spaces.  
    - stderr file will tell you what samples were left out and why, as well as which samples (by first two parts prefix) were renamed.
4) If needed, download from cavatica any `.strelka.vep.maf` files required using the sheet from step one.  Ths script will then concatenate individual mafs by disease, while also excluding standard cbioportal excluded `Variant_Classification`.
    - It will also ignore `Entrez_Gene_ID` as the reference seems to have a bug in ID assignment to HUGO gene symbol.
    - Also, the genome reference build is ignored as cbioportal current won't accept anything other than build 37.
    - BS IDs are replaced with the samples IDs used in the data sheets
5) If needed, download from cavatica any `.CNVs` files required using the sheet from step one. 5a and 5b will do the following respectively:
    - Convert genome-coordinated cnv info to gene-based.
    - Same as number 4, except with copy number files.  They are merged into a table, samples renamed.
6) If needed, download from cavatica any `.genes.results.gz` files required using the sheet from step one. Merge all rsem files into a table, similar to CNV table.  Will be used as input for Load scripts

7) Created by Pichai, these scripts will reformat cnv tables to continuous and discrete values, create all case lists and project (disease) descriptions, and create RNA expression and Z score files based on the RNA table from step 3 above.
- It will also copy and organized all the created and reformatted DNA files.
- Will output completed project into a `processed` subdirectory
- Uses `CNALoadingHelper.R` and `RNALoadingHelper.R` when project is flagged as having these data
 
Last step is to basically copy all new file output into the repo, and push as described above.  Kartik has implemented automatic checks to see if validation would pass for cbioportal upon pull request.
