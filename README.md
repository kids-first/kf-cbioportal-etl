# Outline on ETL for converting data from cavatica and data service to pedcbioportal format
In general, we are creating upload packages converting our data and metadata to satisfy the requirements outlined [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats).
Further general loading notes can be found in this [Notion page](https://www.notion.so/d3b/Cbioportal-Study-Load-SOP-58812479fabe4d2fa9f72242e331b5ee).
See [below](#collaborative-and-publication-workflows) for special cases like publications or collaborative efforts
## I have everything and I know I am doing
Below assumes you have already created the necessary tables from dbt
1. Run commands as outlined in `scripts/get_study_metadata.py`. Copy/move those files to the cBio loader ec2 instance
1. Copy over the approriate aws account key and download files. Example using `pbta_all` study:

   ```sh
   python3 ~/tools/kf-cbioportal-etl/scripts/get_files_from_manifest.py -m genomics_file_manifest.txt -f RSEM_gene,annofuse_filtered_fusions_tsv,annotated_public_outputs,ctrlfreec_pval,ctrlfreec_info,ctrlfreec_bam_seg -p saml 2> pbta_dl.log &
   python3 ~/tools/kf-cbioportal-etl/scripts/get_files_from_manifest.py -m dgd_genomics_file_manifest.txt -f DGD_MAF,DGD_FUSION -p d3b 2> dgd_dl.log &
   ```

1. Copy and edit `REFS/data_processing_config.json` and `REFS/pbta_all_case_meta_config.json` as needed
1. Run pipeline script - ignore manifest section, it is a placeholder for a better function download method

   ```sh
   scripts/genomics_file_cbio_package_build.py -t cbio_file_name_id.txt -c pbta_all_case_meta_config.json -d data_processing_config.json -f both
   ```
1. Check logs and outputs for errors, especially `validator.errs` and `validator.out`, assuming everything else went fine, to see if any `ERROR` popped up that would prevent the pakcage from loading properly once pushed to the bucket and Jenkins import job is run

### Final output example

In the end, if you named your output dir `processed`, you'll end up with this example output from `pbta_all` study:
```sh
processed
└── pbta_all
    ├── case_lists
    │   ├── cases_3way_complete.txt
    │   ├── cases_RNA_Seq_v2_mRNA.txt
    │   ├── cases_all.txt
    │   ├── cases_cna.txt
    │   ├── cases_cnaseq.txt
    │   ├── cases_sequenced.txt
    │   └── cases_sv.txt
    ├── data_CNA.txt -> /home/ubuntu/mount/pbta_all/merged_cnvs/pbta_all.discrete_cnvs.txt
    ├── data_clinical_patient.txt -> /home/ubuntu/mount/pbta_all/datasheets/data_clinical_patient.txt
    ├── data_clinical_sample.txt -> /home/ubuntu/mount/pbta_all/datasheets/data_clinical_sample.txt
    ├── data_cna.seg.txt -> /home/ubuntu/mount/pbta_all/merged_cnvs/pbta_all.merged_seg.txt
    ├── data_fusions.txt -> /home/ubuntu/mount/pbta_all/merged_fusion/pbta_all.fusions.txt
    ├── data_gene_matrix_CHOP.txt -> /home/ubuntu/mount/pbta_all/datasheets/data_gene_matrix_CHOP.txt
    ├── data_linear_CNA.txt -> /home/ubuntu/mount/pbta_all/merged_cnvs/pbta_all.predicted_cnv.txt
    ├── data_mutations_extended.txt -> /home/ubuntu/mount/pbta_all/merged_mafs/pbta_all.maf
    ├── data_rna_seq_v2_mrna.txt -> /home/ubuntu/mount/pbta_all/merged_rsem/pbta_all.rsem_merged.txt
    ├── data_rna_seq_v2_mrna_median_Zscores.txt -> /home/ubuntu/mount/pbta_all/merged_rsem/pbta_all.rsem_merged_zscore.txt
    ├── meta_CNA.txt
    ├── meta_FUSION.txt
    ├── meta_clinical_patient.txt
    ├── meta_clinical_sample.txt
    ├── meta_cna.seg.txt
    ├── meta_gene_matrix_CHOP.txt
    ├── meta_linear_CNA.txt
    ├── meta_mutations_extended.txt
    ├── meta_rna_seq_v2_mrna.txt
    ├── meta_rna_seq_v2_mrna_median_Zscores.txt
    └── meta_study.txt
```
## Details
Use this section as a reference in case your overconfidence got the best of you.
## Software Prerequisites

+ `python3` v3.5.3+
  + `numpy`, `pandas`, `scipy`
+ `bedtools` (https://bedtools.readthedocs.io/en/latest/content/installation.html)
+ `chopaws` https://github.research.chop.edu/devops/aws-auth-cli needed for saml key generation for s3 upload
+ access to https://aws-infra-jenkins-service.kf-strides.org to start cbio load into QA and/or prod using the `d3b-center-aws-infra-pedcbioportal-import` task
+ Access to the `postgres` D3b Warehouse database at `d3b-warehouse-aurora-prd.d3b.io`. Need at least read access to tables with the `bix_workflows` schema
+ [cbioportal git repo](https://github.com/cBioPortal/cbioportal) needed to validate the final study output

## Starting file inputs
Most starting files are exported from the D3b Warehouse. An example of file exports can be found here `scripts/export_clinical.sh`, we now use `scripts/get_study_metadata.py` to get the files.
However, a python wrapper script that leverages the `x_case_meta_config.json` is recommended to use for each study.

### scripts/get_study_metadata.py
```
usage: get_study_metadata.py [-h] [-d DB_INI] [-p PROFILE] [-c CONFIG_FILE]

Pull clinical data and genomics file etl support from D3b data warehouse.

optional arguments:
  -h, --help            show this help message and exit
  -d DB_INI, --db-ini DB_INI
                        Database config file - formatting like aws or sbg creds
  -p PROFILE, --profile PROFILE
                        ini profile name
  -c CONFIG_FILE, --config CONFIG_FILE
                        json config file with meta information; see REFS/pbta_all_case_meta_config.json example
```

### From D3b Warehouse
#### - Genomic files manifest
This is a s3 manifest of all files to loaded onto the portal.
It is generally created by Bix-Ops and loaded into the D3b Warehouse.
If the study is combining a KF/PBTA study with DGD, you may need to download a second manifest.
#### - Data clinical sample sheet
This is the cBioportal-formatted sample sheet that follows guidelines from [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#clinical-sample-columns)
#### - Data clinical patient sheet
This is the cBioportal-formatted patient sheet that follows guidelines from [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#clinical-patient-columns)
#### - Genomics metadata file
Seemingly redundant, this file contains the file locations, BS IDs, file type, and cBio-formatted sample IDs of all inputs.
It helps simplify the process to integrate better into the downstream tools.
This is the file that goes in as the `-t` arg in all the data collating tools
#### - Sequencing center info resource file
This is a simple file this BS IDs and sequencing center IDs and locations.
It is necessary to patch in a required field for the fusion data
#### - Data gene matrix - *OPTIONAL*
This is only required if you have a custom panel - like the DGD does
### User-edited
#### - Data processing config file

This is a json formatted file that has tool paths, reference paths, and run time params.
An example is given in `REFS/data_processing_config.json`.
This section here:
```json
"file_loc_defs": {
    "_comment": "edit the values based on existing/anticipated source file locations, relative to working directory of the script being run",
    "mafs": {
      "kf": "annotated_public_outputs",
      "dgd": "DGD_MAF",
      "header": "/home/ubuntu/tools/kf-cbioportal-etl/REFS/maf_KF_CONSENSUS.txt"
    },
    "cnvs": {
      "pval": "ctrlfreec_pval",
      "info": "ctrlfreec_info",
      "seg": "ctrlfreec_bam_seg"
    },
    "rsem": "RSEM_gene",
    "fusion": "annofuse_filtered_fusions_tsv",
    "dgd_fusion": "DGD_FUSION",
    "fusion_sq_file": ""
  },
  "dl_file_type_list": ["RSEM_gene","annofuse_filtered_fusions_tsv","annotated_public_outputs",
    "ctrlfreec_pval","ctrlfreec_info","ctrlfreec_bam_seg", "DGD_MAF"],
```
Will likely need the most editing existing based on your input, and should only need to updated if something changes after initial load.

#### - Metadata processing config file

This is a json config file with file descriptions and case lists required by the cbioportal.
An example is given in `REFS/pbta_all_case_meta_config.json`.
Within this file is a `_doc` section with a decent explanation of the file format and layout.
Be sure to review all data types to be loaded by review all `meta_*` to see if they match incoming data.
Likely personalized edits would occur in the following fields:
+ `merged_{data type}`: The `profile_description` key in each is a good place to describe any algorithm or nuances used to generate the data of that type. Also be sure to remove any data types not being loaded, as that determines what genomic file collation steps are run.
+ `study`: Here is where you set the overall study description, it's the banner text that people will see in the study overview page that gives them a sense of what the data is.
  + `description`: This field is set up as an array so that a generic form of "text describing" "disease" "more text describing" can be used. Put another way, element one is whatever you want to say about the disease/study until you are ready to mention the disease/study, element two anything you may optionally wish to add
  + `groups`: These are access groups defined is cBioportal.  Default is `PUBLIC`, but another can be named is restrictions are needed.  Need to work with Devops for custom groups
  + `cancer_study_identifier`: This is the short name that you create for the study.  It will be the name of the study load folder and will be used by cBioportal to find all relevant information for that study.
  + `type_of_cancer`: This is the oncotree code used to categorize the study to a disease type that best summarizes all samples in the study. These are the default codes: http://oncotree.mskcc.org/#/home. Internally, we have added `phgg` and `plgg`. If your study doesn't fit, propose a new one to be added
  + `display_name`: This is what will show as a long form title on the site home page
  + `short_name`: This is the short version. By default, should be the same as `cancer_study_identifier`


## Pipeline script
After downloading the genomic files and files above as needed, and properly editing config files as needed, this script should generate and validate the cBioportal load package

### scripts/genomics_file_cbio_package_build.py
```
usage: genomics_file_cbio_package_build.py [-h] [-t TABLE] [-m MANIFEST] [-c CBIO_CONFIG] [-d DATA_CONFIG] [-f [{both,kf,dgd}]]

Download files (if needed), collate genomic files, organize load package.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file names
  -m MANIFEST, --manifest MANIFEST
                        Download file manifest, if needed
  -c CBIO_CONFIG, --cbio-config CBIO_CONFIG
                        cbio case and meta config file
  -d DATA_CONFIG, --data-config DATA_CONFIG
                        json config file with data types and data locations
  -f [{both,kf,dgd}], --dgd-status [{both,kf,dgd}]
                        Flag to determine load will have pbta/kf + dgd(both), kf/pbta only(kf), dgd-only(dgd)
```
+ `-t` would be the [Genomics metadata file](#genomics-metadata-file)
+ `-c` would be the [Metadata processing config file](#metadata-processing-config-file)
+ `-d` would be the [Data processing config file](#data-processing-config-file)

Check the pipeline log output for any errors that might have occurred.

## Upload the final packages
Upload all of the directories named as study short names to `s3://kf-cbioportal-studies/public/`. You may need to set and/or copy aws your saml key before uploading. Next, edit the file in that bucket called `importStudies.txt` located at `s3://kf-cbioportal-studies/public/importStudies.txt`, with the names of all of the studies you wish to updated/upload. Lastly, go to https://jenkins.kids-first.io/job/d3b-center-aws-infra-pedcbioportal-import/job/master/, click on build. At the `Promotion kf-aws-infra-pedcbioportal-import-asg to QA` and `Promotion kf-aws-infra-pedcbioportal-import-asg to PRD`, the process will pause, click on the box below it to affirm that you want these changes deployed to QA and/or PROD respectively.  If both, you will have to wait for the QA job to finish first before you get the prompt for PROD.
## Congratulations, you did it!

# Collaborative and Publication Workflows
These are highly specialized cases in which all or most of the data come from a third party, and therefore requires specific case-by-case protocols.

## OpenTargets
This project is organized much like OpenPBTA in which all genomics data for each assay-type are collated into one giant table.
In general, this fits cBioPortal well.
Input files mostly come from a "subdirectory" from within `s3://kf-openaccess-us-east-1-prd-pbta/`, consisting of:
 - `histologies.tsv`
 - `snv-consensus-plus-hotspots.maf.tsv.gz`
 - `consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz`
 - `consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz`
 - `gene-expression-rsem-tpm-collapsed.rds`
 - `fusion-putative-oncogenic.tsv`

See `https://pedcbioportal.kidsfirstdrc.org/study/summary?id=ped_opentargets_2021` for final product.

## Prep work
The histologies file needs `formatted_sample_id` added and likely a blacklist from the D3b Warehouse or some other source to supress duplicate RNA libraries from different sequencing methods.
Since we are not handling `Methylation` yet, it is recommneded those entries be removed ahead of time.
To create the histologies file, recommended method is to:
1. `docker pull pgc-images.sbgenomics.com/d3b-bixu/open-pedcan:latest` if you haven't already
1. Pull the OpenPedCan repo (should probably make the script more flexible pulling a 12GB repo for such a small task is a bit overkill): https://github.com/PediatricOpenTargets/OpenPedCan-analysis
1. Export from D3b Warehouse the latest existing cBio IDs to use for population. Ensure that the output is csv double-quoted. Copy that into the OpenPedCan-analysis repo in `analyses/pedcbio-sample-name/input/cbio_name.csv`. Currently that can be obtained using the sql command:
    ```sql

    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.aml_sd_pet7q6f2_2018_cbio_sample
    union
    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.aml_sd_z6mwd3h0_2018_cbio_sample
    union
    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.x01_fy16_nbl_maris_cbio_sample
    union
    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.cbtn_cbio_sample
    union
    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.dgd_cbio_sample
    union
    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.pnoc_cbio_sample
    union
    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.oligo_nation_cbio_sample
    union
    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.os_sd_zxjffmef_2015_cbio_sample

    ```
1. Run an interactive docker, and ensure to mount a volume that will have the repo and whatever input histologies file you end up using, i.e. `docker run -it --mount type=bind,source=/home/ubuntu,target=/WORK pgc-images.sbgenomics.com/d3b-bixu/open-pedcan:latest /bin/bash`
1. In that container, go to the location of `analyses/pedcbio-sample-name/pedcbio_sample_name_col.R`
1. Get a blacklist from D3b Warehouse, exporting table `bix_workflows.cbio_hide_reasons`
1. Run `Rscript --vanilla pedcbio_sample_name_col.R --hist_dir path-to-hist-dir`. Histologies file must be `histologies.tsv`, modify file name or create sym link if needed. Results will be in `results` as `histologies-formatted-id-added.tsv`

### Inputs
Inputs are located in the old Kids First AWS account (`538745987955`) in this general bucket location: `s3://kf-openaccess-us-east-1-prd-pbta/open-targets/`.
Clinical data with cBio names are obtained from the `histologies-formatted-id-added.tsv` file, as noted in [Prep Work section](#prep-work).
Genomic data generally obtained as such:
 - Somatic variant calls: merged maf
 - Copy number: tsv file with copy number, ploidy, and GISTIC-style information in maf-like format (each call is a row)
 - RNA expression: tpm values from rsem stored an `.rds` object
 - RNA fusion: annoFuse output

### File Transformation
It's recommended to datasheets in a dir called `datasheets`, and the rest of the outputs into it's own dir to keep things sane and also be able to leverage existing study build script in `scripts/organize_upload_packages.py`
#### 1. COLLABORATIONS/openTARGETS/clinical_to_datasheets.py
 ```
usage: clinical_to_datasheets.py [-h] [-f HEAD] [-c CLIN] [-s CL_SUPP]

Script to convert clinical data to cbio clinical data sheets

optional arguments:
  -h, --help            show this help message and exit
  -f HEAD, --header-file HEAD
                        tsv file with input file original sample names, output
                        sheet flag, and conversion
  -b BLACKLIST, --blacklist BLACKLIST
                        because every club needs a bouncer. Headered tsv file with BS ID and reason
  -c CLIN, --clinical-data CLIN
                        Input clinical data sheet
  -s CL_SUPP, --cell-line-supplement CL_SUPP
                        supplemental file with cell line meta data - bs
                        id<tab>type. optional
 ```
 - `-f` Header file example can be found here: `COLLABORATIONS/openTARGETS/header_desc.tsv`
 - `-s` cell line supplemental file here: `REFS/cell_line_supplemental.txt`
 - `-b` blacklist exported from D3b Warehouse
 - `-c` `histologies-formatted-id-added.tsv`

 Outputs a `data_clinical_sample.txt` and `data_clinical_patient.txt` for the cBio package, and a `bs_id_sample_map.txt` mapping file to link BS IDs to gnerated cBioPortal IDs based on the rules for creating a proper somatic event using column `parent_aliquot_id`

Example run:
`python3 COLLABORATIONS/openTARGETS/clinical_to_datasheets.py -f COLLABORATIONS/openTARGETS/header_desc.tsv -c histologies-formatted-id-added.tsv -b cbio_hide_reasons.tsv 2> clin2.errs`

#### 2. COLLABORATIONS/openTARGETS/rename_filter_maf.py
_NOTE_ for v11 input, I ran the following command `zcat snv-dgd.maf.tsv.gz | perl -e '$skip = <>; $skip= <>; while(<>){print $_;}' | gzip -c >> snv-consensus-plus-hotspots.maf.tsv.gz`
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
`python3 COLLABORATIONS/openTARGETS/cnv_to_tables.py -m bs_id_sample_map.txt  -c consensus_wgs_plus_cnvkit_wxs.tsv.gz`

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