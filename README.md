# Outline on ETL for converting data from cavatica and data service to pedcbioportal format
In general, we are creating upload packages converting our data and metadata to satisfy the requirements outlined [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats).
Further general loading notes can be found in this [Notion page](https://www.notion.so/d3b/Cbioportal-Study-Load-SOP-58812479fabe4d2fa9f72242e331b5ee).
See [below](#collaborative-and-publication-workflows) for special cases like publications or collaborative efforts
## I have everything and I know I am doing
Below assumes you have already created the necessary tables from dbt
1. Run commands as outlined in `scripts/get_study_metadata.py`. Copy/move those files to the cBio loader ec2 instance
1. Copy over the appropriate aws account key and download files. Example using `pbta_all` study:

   ```sh
    python3 scripts/get_files_from_manifest.py -m cbtn_genomics_file_manifest.txt,pnoc_genomics_file_manifest.txt,x01_genomics_file_manifest.txt,dgd_genomics_file_manifest.txt -f RSEM_gene,annofuse_filtered_fusions_tsv,annotated_public_outputs,ctrlfreec_pval,ctrlfreec_info,ctrlfreec_bam_seg,annotated_public -t aws_buckets_key_pairs.txt -s turbo -c cbio_file_name_id.txt
   ```
  `aws_bucket_key_pairs.txt` is a headerless tsv file with bucket name and aws profile name pairs

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
    ├── data_linear_CNA.txt -> /home/ubuntu/mount/pbta_all/merged_cnvs/pbta_all.predicted_cnv.txt
    ├── data_mutations_extended.txt -> /home/ubuntu/mount/pbta_all/merged_mafs/pbta_all.maf
    ├── data_rna_seq_v2_mrna.txt -> /home/ubuntu/mount/pbta_all/merged_rsem/pbta_all.rsem_merged.txt
    ├── data_rna_seq_v2_mrna_median_Zscores.txt -> /home/ubuntu/mount/pbta_all/merged_rsem/pbta_all.rsem_merged_zscore.txt
    ├── data_sv.txt -> /home/ubuntu/mount/pbta_all/pbta_all.fusions.txt
    ├── meta_CNA.txt
    ├── meta_SV.txt
    ├── meta_clinical_patient.txt
    ├── meta_clinical_sample.txt
    ├── meta_cna.seg.txt
    ├── meta_linear_CNA.txt
    ├── meta_mutations_extended.txt
    ├── meta_rna_seq_v2_mrna.txt
    ├── meta_rna_seq_v2_mrna_median_Zscores.txt
    └── meta_study.txt
```
# Details
Use this section as a reference in case your overconfidence got the best of you

## REFS
In case you want to use different reference inputs...
 - From data_processing_config.json `bed_genes`:
   - This is used to collate ControlFreeC results into gene hits
   - For VEP 105, gtf was downloaded from `https://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz`
   - Then, using bedops and a perl-one-liner: 
   ```sh
   cat Homo_sapiens.GRCh38.105.chr.gtf | perl -e 'while(<>){@a=split /\t/; if($a[2] eq "gene" && $a[8] =~ /gene_name/){print $_;}}'  | convert2bed -i gtf --attribute-key=gene_name  > Homo_sapiens.GRCh38.105.chr.gtf_genes.bed
   ```
To get aws bucket prefixes to add key (meaning aws profile names) to:
```sh
cat *genomic* | cut -f 15 | cut -f 1-3 -d "/" | sort | uniq > aws_bucket_key_pairs.txt
```
Just remove the `s3_path` and `None` entries


## Software Prerequisites

+ `python3` v3.5.3+
  + `numpy`, `pandas`, `scipy`
+ `bedtools` (https://bedtools.readthedocs.io/en/latest/content/installation.html)
+ `chopaws` https://github.research.chop.edu/devops/aws-auth-cli needed for saml key generation for s3 upload
+ access to https://github.com/d3b-center/aws-infra-pedcbioportal-import repo. To start a load job:
  + Create a branch and edit the `import_studies.txt` file with the study name you which to load. Can be an MSKCC datahub link or a local study name
  + Push the branch to remote - this will kick off a github action to load into QA
  + To load into prod, make a PR. On merge, load to prod will kick off
  + aws `stateMachinePedcbioImportservice` Step function service is used to view and mangage running jobs
  + To repeat a load, click on the ▶️ icon in the git repo to select the job you want to re-run
  + *Note*, if your branch importStudies.txt is the same as main, you may have tot rigger it yourself. To do so, go to [actions](https://github.com/d3b-center/aws-infra-pedcbioportal-import/actions), on the left panel choose which action you want, then from the drop down in the right panel, pick which branch you want that action to run on 
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
DEPRECATED and will be removed from future releases
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
      "header": "/home/ubuntu/tools/kf-cbioportal-etl/REFS/maf_KF_CONSENSUS.txt"
    },
    "cnvs": {
      "pval": "ctrlfreec_pval",
      "info": "ctrlfreec_info",
      "seg": "ctrlfreec_bam_seg"
    },
    "rsem": "RSEM_gene",
    "fusion": "annofuse_filtered_fusions_tsv",
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

### scripts/get_files_from_manifest.py
Currently, file locations are still too volatile to trust to make downloading part of the pipeline. Using various combinations of buckets and sbg file ID pulls will eventually get you everything
```
usage: get_files_from_manifest.py [-h] [-m MANIFEST] [-f FTS] [-p PROFILE] [-s SBG_PROFILE] [-c CBIO] [-a] [-d]

Get all files for a project.

optional arguments:
  -h, --help            show this help message and exit
  -m MANIFEST, --manifest-list MANIFEST
                        csv list of of genomic file location manifests
  -f FTS, --file-types FTS
                        csv list of workflow types to download
  -p PROFILE, --profile PROFILE
                        aws profile name. Leave blank if using sbg instead
  -s SBG_PROFILE, --sbg-profile SBG_PROFILE
                        sbg profile name. Leave blank if using AWS instead
  -c CBIO, --cbio CBIO  Add cbio manifest to limit downloads
  -a, --active-only     Set to grab only active files. Recommended.
  -d, --debug           Just output manifest subset to see what would be grabbed
```
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
  -l, --legacy          If set, will run legacy fusion output
```
+ `-t` would be the [Genomics metadata file](#genomics-metadata-file)
+ `-c` would be the [Metadata processing config file](#metadata-processing-config-file)
+ `-d` would be the [Data processing config file](#data-processing-config-file)

Check the pipeline log output for any errors that might have occurred.

## Upload the final packages
Upload all of the directories named as study short names to `s3://kf-cbioportal-studies/public/`. You may need to set and/or copy aws your saml key before uploading. Next, edit the file in that bucket called `importStudies.txt` located at `s3://kf-cbioportal-studies/public/importStudies.txt`, with the names of all of the studies you wish to updated/upload. Lastly, follow the directions reference in [Software Prerequisites](#software-prerequisites) to load the study.
## Congratulations, you did it!

# Collaborative and Publication Workflows
These are highly specialized cases in which all or most of the data come from a third party, and therefore requires specific case-by-case protocols.

## OpenTargets
This project is organized much like OpenPBTA in which all genomics data for each assay-type are collated into one giant table.
In general, this fits cBioPortal well.
Input files mostly come from a "subdirectory" from within `s3://d3b-openaccess-us-east-1-prd-pbta/open-targets/`, consisting of:
 - `histologies.tsv`
 - `snv-consensus-plus-hotspots.maf.tsv.gz`
 - `consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz`
 - `consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz`
 - `gene-expression-rsem-tpm-collapsed.rds`
 - `fusion-putative-oncogenic.tsv`

See `https://pedcbioportal.kidsfirstdrc.org/study/summary?id=ped_opentargets_2021` for final product.

## Prep work
The histologies file needs `formatted_sample_id` added and likely a blacklist from the D3b Warehouse or some other source to supress duplicate RNA libraries from different sequencing methods.
Since we are not handling `Methylation` yet, it is recommended those entries be removed ahead of time.
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
	  union
    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.chdm_sd_7spqtt8m_cbio_sample
	  union
    select participant_id, formatted_sample_id, specimen_id, analyte_types, normal_bs_id, normal_sample_id
    from prod_cbio.chdm_phs001643_2018_cbio_sample

    ```
1. Run an interactive docker, and ensure to mount a volume that will have the repo and whatever input histologies file you end up using, i.e. `docker run -it --mount type=bind,source=/home/ubuntu,target=/WORK pgc-images.sbgenomics.com/d3b-bixu/open-pedcan:latest /bin/bash`
1. In that container, go to the location of `analyses/pedcbio-sample-name/pedcbio_sample_name_col.R`
1. Get a blacklist from D3b Warehouse, exporting table `bix_workflows.cbio_hide_reasons`
1. Run `Rscript --vanilla pedcbio_sample_name_col.R --hist_dir path-to-hist-dir`. Histologies file must be `histologies.tsv`, modify file name or create sym link if needed. Results will be in `results` as `histologies-formatted-id-added.tsv`

### Inputs
Inputs are located in the old D3b AWS account (`684194535433`) in this general bucket location: `s3://d3b-openaccess-us-east-1-prd-pbta/open-targets/`.
Clinical data with cBio names are obtained from the `histologies-formatted-id-added.tsv` file, as noted in [Prep Work section](#prep-work).
Genomic data generally obtained as such:
 - Somatic variant calls: merged maf
 - Copy number: tsv file with copy number, ploidy, and GISTIC-style information in maf-like format (each call is a row)
 - RNA expression: tpm values from rsem stored an `.rds` object
 - RNA fusion: annoFuse output
For example, for v12, bucket s3://d3b-openaccess-us-east-1-prd-pbta/open-targets/v12/:
```
consensus_wgs_plus_cnvkit_wxs.tsv.gz
fusion-dgd.tsv.gz
fusion-putative-oncogenic.tsv
gene-expression-rsem-tpm-collapsed.rds
tcga-gene-expression-rsem-tpm-collapsed.rds
snv-consensus-plus-hotspots.maf.tsv.gz
snv-dgd.maf.tsv.gz
```

### File Transformation
It's recommended to put datasheets in a dir called `datasheets`, and the rest of the outputs into it's own dir to keep things sane and also be able to leverage existing study build script in `scripts/organize_upload_packages.py`
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
  -s SKIP, --skip SKIP  Skip typical #version header line
  -n TYPE, --study TYPE
                        study name, like "openpbta"
```

Example run:
`python3 COLLABORATIONS/openTARGETS/rename_filter_maf.py -m bs_id_sample_map.txt -v snv-consensus-plus-hotspots.maf.tsv.gz -s 1 -n openpedcan_v11`

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
  -s TYPE, --study TYPE
                        study name, like "openpbta"
```

Example run:
`python3 COLLABORATIONS/openTARGETS/cnv_to_tables.py -m bs_id_sample_map.txt  -c consensus_wgs_plus_cnvkit_wxs.tsv.gz -s openpedcan_v11`

#### 4. COLLABORATIONS/openTARGETS/rename_export_rsem.R
Note, I merged the tcga into the main rds. I also needed an instance with _64GB ram_ in order to calc z scores. Update: Can also achieve by setting up 32GB swap space 
Additionally,can also use C++ implementation to speed up zscore computation and save memory. Note: writing in C++ is slower than write_tsv function in R but have less memory footprint.
```
Usage: /home/ubuntu/tools/kf-cbioportal-etl/COLLABORATIONS/openTARGETS/rename_export_rsem.R [options]


Options:
	--rna_rds=RNA_RDS
		openX rsem rds expression object

	--map_id=MAP_ID
		mapping ID file with headers: BS_ID	Sample Type	Cbio ID

	--type=TYPE
		study name, like 'openpbta'

	--computeZscore=R/C++
		Implementation to use C++ or R to compute zscore and write a tsv file

	-h, --help
		Show this help message and exit
```
Example run:
`Rscript COLLABORATIONS/openTARGETS/rename_export_rsem.R --rna_rds gene_tcga_expression_common_merge.rds --map_id bs_id_sample_map.txt --type openpedcan_v11 --computeZscore R 2> rna_convert.errs`

### DEPRECATED
Leaving here until next major release. Refer to [convert sv as fusion](#5-scriptsconvert_fusion_as_svpy)
#### 5a. scripts/rna_convert_fusion.py
Before running, to leverage an existing fusion conversion, I first ran:
`COLLABORATIONS/openTARGETS/reformat_cbio_sample_index.py -t bs_id_sample_map.txt -n openpedcan_v11 > fusion_sample_name_input.txt`
to reformat the sample name index. Then actually ran the rna fusion script:

```
Convert openPBTA fusion table OR list of annofuse files to cbio format.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file names
  -f FUSION_RESULTS, --fusion-results FUSION_RESULTS
                        openPBTA fusion file OR annoFuse results dir
  -m MODE, --mode MODE  describe source, pbta or annofuse
  -s SQ_FILE, --center-file SQ_FILE
                        File with BS IDs and sequencing centers. Should have headered columns: BS_ID SQ_Value
  -o OUT_DIR, --out-dir OUT_DIR
                        Result output dir. Default is merged_fusion
```
Example run:
`scripts/rna_convert_fusion.py -t fusion_sample_name_input.txt -f fusion-putative-oncogenic.tsv -m pbta -s COLLABORATIONS/openTARGETS/seq_center_info_updated.txt`

#### 5b. /scripts/add_dgd_fusion.py
Since dgd was just added, this appends to the existing file as it has different columns and input format
```
Output fields DGD fusion - meant to be appended to an existing file!

optional arguments:
  -h, --help            show this help message and exit
  -f FUSION_DIR, --fusion-dir FUSION_DIR
                        Fusion file directory
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file names
  -a, --append          Optional - if given will output to stdout to append, else will create new merged file
  -m, --merged          If input is already merged, treat fusion dir as file instead
  -o OUT_DIR, --out-dir OUT_DIR
                        Result output dir. Default is merged_fusion
```
Example run:
`scripts/add_dgd_fusion.py -f fusion-dgd.tsv.gz -t fusion_sample_name_input.txt -a -m >> merged_fusion/openpedcan_v11.fusions.txt`
#### 5. scripts/convert_fusion_as_sv.py
```
usage: convert_fusion_as_sv.py [-h] [-t TABLE] [-f FUSION_RESULTS] [-o OUT_DIR] -m MODE

Convert openPBTA fusion table OR list of annofuse files to cbio format.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file names
  -f FUSION_RESULTS, --fusion-results FUSION_RESULTS
                        annoFuse results dir OR openX merged fusion file
  -o OUT_DIR, --out-dir OUT_DIR
                        Result output dir. Default is merged_fusion
  -m MODE, --mode MODE  describe source, openX or kfprod
```
Example run:
`python3 ~/tools/kf-cbioportal-etl/scripts/convert_fusion_as_sv.py -t fusion_sample_table.txt -m openX -f ../pbta-fusion-putative-oncogenic.tsv`
#### 6. scripts/organize_upload_packages.py
Leverage the existing meta config and package organizer from kids first to create all relevant meta and case files...except for the added case lists achieved by the next step
```
Create cases lists, meta files, and organize data for cbio upload. It is assumed you are at the dir level of all input data files

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_DIR, --output_dir OUT_DIR
                        output directory name
  -c CONFIG_FILE, --config CONFIG_FILE
                        json config file with meta information; see REFS/case_meta_config.json example
```
Example run:
`python3 scripts/organize_upload_packages.py -o processed -c COLLABORATIONS/openTARGETS/openpedcan_v11_case_meta_config.json`

#### 7. COLLABORATIONS/openTARGETS/case_list_from_datasheet.py
Last step before validation and upload
```
usage: case_list_from_datasheet.py [-h] [-d DATASHEET] [-s STUDY_ID]

Generate extra cases lists based on cBio datasheet

optional arguments:
  -h, --help            show this help message and exit
  -d DATASHEET, --datasheet DATASHEET
                        cBio sample datasheet
  -s STUDY_ID, --study-id STUDY_ID
                        cBio cancer_study_identifier
  -c COHORTS, --cohort-csv COHORTS
                        add a special case for a cohort-specific case list using a csv list and omit from hist split
  -m SAMP_MIN, --min-sample SAMP_MIN
                        min number of samples a cohort must have to generate a case list for. recommend 3
```

Example run:
`python3 COLLABORATIONS/openTARGETS/case_list_from_datasheet.py -d data_clinical_sample.txt -s openpedcan_v11 -c GTEx -m 3`

## OpenPBTA
Quite similar to openTargets with some caveats

### Prep Work
1. Obtain files from `s3://d3b-openaccess-us-east-1-prd-pbta/data/release-xxx`, with the `xxx` being a release version of interest.
  Files to get includes:
  ```
  pbta-histologies.tsv
  consensus_seg_annotated_cn_x_and_y.tsv.gz
  consensus_seg_annotated_cn_autosomes.tsv.gz
  pbta-snv-consensus-mutation.maf.tsv.gz
  pbta-snv-scavenged-hotspots.maf.tsv.gz
  pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
  pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
  pbta-fusion-putative-oncogenic.tsv.gz
  pbta-cnv-consensus.seg.gz
  ```
1. Create complete maf:
  `python3 COLLABORATIONS/openPBTA/merge_pbta_file_columns.py -1 pbta-snv-consensus-mutation.maf.tsv.gz -2 pbta-snv-scavenged-hotspots.maf.tsv.gz -o pbta-snv-consensus-plus-hotspots.maf`
1. Create complete rsem rds:
  `Rscript COLLABORATIONS/openPBTA/merge_rsem_rds.R --stranded pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --polya pbta-gene-expression-rsem-fpkm-collapsed.polya.rds`
1. Create complete cnv tsv:
  `python3 COLLABORATIONS/openPBTA/merge_pbta_file_columns.py -1 consensus_seg_annotated_cn_autosomes.tsv.gz -2 consensus_seg_annotated_cn_x_and_y.tsv.gz -o consensus_seg_annotated_cn.tsv`
### Create study
Uses many of the same scripts as open targets
#### 1. Create clinical data sheets
Histologies file is parsed and scplit into `data_clinical_sample.txt` and `data_clinical_patient.txt`
```sh
python3 COLLABORATIONS/openPBTA/clinical_to_datasheets.py -h
usage: clinical_to_datasheets.py [-h] [-f HEAD] [-c CLIN] [-b BLACKLIST] [-p PBTA_DS]

Script to convert clinical data to cbio clinical data sheets

optional arguments:
  -h, --help            show this help message and exit
  -f HEAD, --header-file HEAD
                        tsv file with input file original sample names, output sheet flag, and conversion
  -c CLIN, --clinical-data CLIN
                        Input clinical data sheet
  -b BLACKLIST, --blacklist BLACKLIST
                        because every club needs a bouncer. Headered tsv file with BS ID and reason
  -p PBTA_DS, --pbta-all-data-sample-sheet PBTA_DS
                        pbta_all sample sheet to use as sample name source
```
Example: `python3 ~/tools/kf-cbioportal-etl/COLLABORATIONS/openPBTA/clinical_to_datasheets.py -f ~/tools/kf-cbioportal-etl/COLLABORATIONS/openPBTA/header_desc.tsv -c pbta-histologies.tsv -p /home/ubuntu/mount/pbta_all/processed/pbta_all/data_clinical_sample.txt -b cbio_hide_reasons_subset.tsv 2> clin2.err`
cBio hide reasons table is exported from D3b warehouse
#### 2. Rename and filter maf file
Filter out Silent, Intron, IGR, 3'UTR, 5'UTR, 3'Flank, 5'Flank, RNA and rename samples to cBio IDs
```sh
python3 COLLABORATIONS/openTARGETS/rename_filter_maf.py -h
usage: rename_filter_maf.py [-h] [-m MAPPING_FILE] [-v MAF_FILE] [-s SKIP] [-n TYPE]

Script to pre-filter entries on usually removed criteria except TERT promoter, convert BS IDs to cBio names

optional arguments:
  -h, --help            show this help message and exit
  -m MAPPING_FILE, --mapping-file MAPPING_FILE
                        tsv file with header and bs_id, sample type, cbio ID mappings
  -v MAF_FILE, --maf-file MAF_FILE
                        openX maf file
  -s SKIP, --skip SKIP  Skip typical #version header line
  -n TYPE, --study TYPE
                        study name, like "openpbta"
```
Example: `python3 /home/ubuntu/tools/kf-cbioportal-etl/COLLABORATIONS/openTARGETS/rename_filter_maf.py -m ../bs_id_sample_map.txt -v ../pbta-snv-consensus-plus-hotspots.maf.tsv -n openpbta`
Mapping file is an output from step 1
#### Convert CNV to cBio CNV tables
Predicted: Just plain copy number, discrete: GISTIC-style assignment
```sh
python3 ~/tools/kf-cbioportal-etl/COLLABORATIONS/openTARGETS/cnv_to_tables.py -h
usage: cnv_to_tables.py [-h] [-m MAPPING_FILE] [-c CNV_TBL] [-s TYPE]

Script to convert openX cnv table to cBio format

optional arguments:
  -h, --help            show this help message and exit
  -m MAPPING_FILE, --mapping-file MAPPING_FILE
                        tsv file with header and bs_id, sample type, cbio ID mappings
  -c CNV_TBL, --copy-number CNV_TBL
                        openX tables, listed as csv since they tend to seprate autosomes, xy
  -s TYPE, --study TYPE
                        study name, like "openpbta"
```
Example: `python3 COLLABORATIONS/openTARGETS/cnv_to_tables.py -m ../bs_id_sample_map.txt -c ../consensus_seg_annotated_cn.tsv -s openpbta`
#### Rename rsem
Mostly just rename samples and export as text. Also z score output
```sh
Rscript COLLABORATIONS/openTARGETS/rename_export_rsem.R --help
Loading required package: data.table
Loading required package: optparse
Loading required package: Rcpp
Loading required package: funr
Loading required package: readr
Usage: COLLABORATIONS/openTARGETS/rename_export_rsem.R [options]


Options:
        --rna_rds=RNA_RDS
                openX rsem rds expression object

        --map_id=MAP_ID
                mapping ID file with headers: BS_ID    Sample Type    Cbio ID

        --type=TYPE
                study name, like 'openpbta'

        --computeZscore=COMPUTEZSCORE
                Use C++ Method to compute zscore and write file. Usage: C++ or R

        -h, --help
                Show this help message and exit
```
Example: `Rscript COLLABORATIONS/openTARGETS/rename_export_rsem.R --rna_rds ../pbta-gene-expression-rsem-fpkm-collapsed.rds --map_id ../bs_id_sample_map.txt --computeZscore C++ --type openpbta`
#### Convert fusion as SV
cBio has a "new" format  with a richer SV input. More data can be ported over.
For compatibility, since standard study conversion script is used, first run: `python3 COLLABORATIONS/openTARGETS/reformat_cbio_sample_index.py -t ../bs_id_sample_map.txt -n openpbta > fusion_sample_table.txt`
Then:
```sh
python3 scripts/convert_fusion_as_sv.py -h
usage: convert_fusion_as_sv.py [-h] [-t TABLE] [-f FUSION_RESULTS] [-o OUT_DIR] -m MODE

Convert openPBTA fusion table OR list of annofuse files to cbio format.

optional arguments:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        Table with cbio project, kf bs ids, cbio IDs, and file names
  -f FUSION_RESULTS, --fusion-results FUSION_RESULTS
                        annoFuse results dir OR openX merged fusion file
  -o OUT_DIR, --out-dir OUT_DIR
                        Result output dir. Default is merged_fusion
  -m MODE, --mode MODE  describe source, openX or kfprod
  ```
  Example: `python3 /home/ubuntu/tools/kf-cbioportal-etl/scripts/convert_fusion_as_sv.py -t fusion_sample_table.txt -f ../pbta-fusion-putative-oncogenic.tsv -m openX`
#### Organize and create links
This script creates meta files and case lists needed for cBio loading
```sh
python3 scripts/organize_upload_packages.py -h
usage: organize_upload_packages.py [-h] [-o OUT_DIR] [-c CONFIG_FILE] [-l]

Create cases lists, meta files, and organize data for cbio upload. It is assumed you are at the dir level of all input data files

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_DIR, --output_dir OUT_DIR
                        output directory name
  -c CONFIG_FILE, --config CONFIG_FILE
                        json config file with meta information; see REFS/case_meta_config.json example
  -l, --legacy          If set, will run legacy fusion output
  ```
  Example: `python3 /home/ubuntu/tools/kf-cbioportal-etl/scripts/organize_upload_packages.py -o processed -c /home/ubuntu/tools/kf-cbioportal-etl/COLLABORATIONS/openPBTA/openpbta_case_meta_config.json`