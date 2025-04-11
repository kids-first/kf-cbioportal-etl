# Whole Study Import - Details
The full study import pipeline is designed to load a complete dataset into PedcBioPortal. It involves generating configuration files, gathering all required metadata and genomic input files, validating the package, and preparing it for upload. This workflow downloads and processes every relevant file from the D3b Warehouse, ensuring a comprehensive and consistent study load. While more resource-intensive than [incremental updates](INCREMENTAL_UPDATES.md), it guarantees that all study components are fully rebuilt and aligned with the latest data.

## REFS
In case you want to use different reference inputs:
 - From study_config.json `bed_genes`:
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

## Config file
This is a json formatted file generated from [generate_config.py](../cbioportal_etl/scripts/generate_config.py) that has tool paths, reference paths, run time params, file descriptions and case lists required by the cBioPortal.
An example is given in `cbioportal_etl/STUDY_CONFIGS/pbta_all_case_meta_config.json` and `cbioportal_etl/STUDY_CONFIGS/json_files/pbta_all_data_processing_config.json`.

### cbioportal_etl/scripts/generate_config.py
```
usage: generate_config.py [-h] [-db DB_INI] [-p PROFILE] [-s STUDY] [-ctsv CONFIG_TSV] [-stsv STUDY_TSV]

Pull clinical data and genomics file etl support from D3b data warehouse.

options:
  -h, --help            show this help message and exit
  -db DB_INI, --db-ini DB_INI
                        Database config file - formatting like aws or sbg creds
  -p PROFILE, --profile PROFILE
                        ini profile name
  -s STUDY, --study STUDY
                        Cancer study ID to compare on server
  -ctsv CONFIG_TSV, --config-tsv CONFIG_TSV
                        Path to the input TSV file containing static values (preset in repo)
  -stsv STUDY_TSV, --study-tsv STUDY_TSV
                        Path to the input TSV file containing completed values for json generation
```
### Custom config values
The config file typically does not need to be modified unless you are running the script outside of the provided tool/repo environment or making custom changes to the default file locations and tool paths. 
Within the config file is a `_doc` section with a decent explanation of the file format and layout.
Be sure to review all data types to be loaded and see if they match incoming data.
Likely personalized edits would occur in the following fields:
+ `merged_{data type}`: The `profile_description` key in each is a good place to describe any algorithm or nuances used to generate the data of that type. Also be sure to remove any data types not being loaded, as that determines what genomic file collation steps are run.
+ `study`: Here is where you set the overall study description, it's the banner text that people will see in the study overview page that gives them a sense of what the data is.
  + `description`: This field is set up as an array so that a generic form of "text describing" "disease" "more text describing" can be used. Put another way, element one is whatever you want to say about the disease/study until you are ready to mention the disease/study, element two anything you may optionally wish to add
  + `groups`: These are access groups defined is cBioPortal.  Default is `PUBLIC`, but another can be named is restrictions are needed.  Need to work with DevOps for custom groups
  + `cancer_study_identifier`: This is the short name that you create for the study.  It will be the name of the study load folder and will be used by cBioPortal to find all relevant information for that study.
  + `type_of_cancer`: This is the oncotree code used to categorize the study to a disease type that best summarizes all samples in the study. [These](http://oncotree.mskcc.org/#/home) are the default codes. Internally, we have added `phgg` and `plgg`. If your study doesn't fit, propose a new one to be added
  + `display_name`: This is what will show as a long form title on the site home page
  + `short_name`: This is the short version. By default, should be the same as `cancer_study_identifier`
  + `file_loc_defs`: This is where the tool paths/reference paths are located.


## Starting file inputs
Most starting files are exported from the D3b Warehouse. An example of file exports can be found here `cbioportal_etl/scripts/export_clinical.sh`, we now use `cbioportal_etl/scripts/get_study_metadata.py` to get the file locations.
The generated `study_config.json` is recommended to use for each study.

### cbioportal_etl/scripts/get_study_metadata.py
```
usage: get_study_metadata.py [-h] [-db DB_INI] [-p PROFILE] [-sc STUDY_CONFIG] [-r REF_DIR] [-a]

Pull clinical data and genomics file etl support from D3b data warehouse.

options:
  -h, --help            show this help message and exit
  -db DB_INI, --db-ini DB_INI
                        Database config file - formatting like aws or sbg creds
  -p PROFILE, --profile PROFILE
                        ini profile name
  -sc STUDY_CONFIG, --study-config STUDY_CONFIG
                        json config file with study information; see cbioportal_etl/STUDY_CONFIGS/json_files for examples
  -r REF_DIR, --ref-dir REF_DIR
                        dir name containing template data_clinical* header files
  -a, --all             flag to include all relevant files, not just status=active files, NOT RECOMMENDED
```

### From D3b Warehouse
#### - Data clinical sample sheet
This is the cBioPortal-formatted sample sheet that follows guidelines from [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#clinical-sample-columns)

#### - Data clinical patient sheet
This is the cBioPortal-formatted patient sheet that follows guidelines from [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#clinical-patient-columns)


## Download starting file inputs 
### cbioportal_etl/scripts/get_files_from_manifest.py
```
usage: get_files_from_manifest.py [-h] [-m MANIFEST] [-f FTS] [-at AWS_TBL] [-sp SBG_PROFILE] [-c CBIO] [-ao] [-rm] [-d] [-o]

Get all files for a project.

options:
  -h, --help            show this help message and exit
  -m MANIFEST, --manifest MANIFEST
                        csv list of of genomic file location manifests
  -f FTS, --file-types FTS
                        csv list of workflow types to download
  -at AWS_TBL, --aws-tbl AWS_TBL
                        Table with bucket name and keys to subset on
  -sp SBG_PROFILE, --sbg-profile SBG_PROFILE
                        sbg profile name. Leave blank if using AWS instead
  -c CBIO, --cbio CBIO  Add cbio manifest to limit downloads. Dn NOT use if using cbio_file_name_id file as the manifest
  -ao, --active-only    Set to grab only active files. Recommended.
  -rm, --rm-na          Remove entries where file_id and s3_path are NA.
  -d, --debug           Just output manifest subset to see what would be grabbed
  -o, --overwrite       If set, overwrite if file exists
```
You can run this script to verify that all required starting files have been downloaded. This also serves as an additional check to confirm that you have access to all necessary files.
### cbioportal_etl/scripts/check_downloads.py
```
usage: check_downloads.py [-h] [-ms MANIFEST]

Check that files were downloaded

options:
  -h, --help            show this help message and exit
  -ms MANIFEST, --manifest-subset MANIFEST
                        tsv list of desired genomic files
```

## Generate and validate cBiopotal load package
After downloading the genomic files and files above as needed, this script should generate and validate the cBioPortal load package

### cbioportal_etl/scripts/genomics_file_cbio_package_build.py
```
usage: genomics_file_cbio_package_build.py [-h] [-m MANIFEST] [-sc STUDY_CONFIG] [-dgd [{both,kf,dgd}]] [-l] [-ad]

Download files (if needed), collate genomic files, organize load package.

options:
  -h, --help            show this help message and exit
  -m MANIFEST, --manifest MANIFEST
                        Download file manifest with cbio project, kf bs ids, cbio IDs, and file names
  -sc STUDY_CONFIG, --study-config STUDY_CONFIG
                        cbio study config file
  -dgd [{both,kf,dgd}], --dgd-status [{both,kf,dgd}]
                        Flag to determine load will have pbta/kf + dgd(both), kf/pbta only(kf), dgd-only(dgd)
  -ad, --add-data       Flag to skip validation when running for add_data directory
```

Check the pipeline log output for any errors that might have occurred.

## Final output example
In the end, you'll end up with this example output from `pbta_all` study in `processed` dir:
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
├── data_CNA.txt
├── data_clinical_patient.txt
├── data_clinical_sample.txt
├── data_clinical_timeline_clinical_event.txt
├── data_clinical_timeline_imaging.txt
├── data_clinical_timeline_specimen.txt
├── data_clinical_timeline_surgery.txt
├── data_clinical_timeline_treatment.txt
├── data_cna.seg.txt
├── data_linear_CNA.txt
├── data_mutations_extended.txt
├── data_rna_seq_v2_mrna.txt
├── data_rna_seq_v2_mrna_median_Zscores.txt
├── data_sv.txt
├── meta_CNA.txt
├── meta_clinical_patient.txt
├── meta_clinical_sample.txt
├── meta_clinical_timeline_clinical_event.txt
├── meta_clinical_timeline_imaging.txt
├── meta_clinical_timeline_specimen.txt
├── meta_clinical_timeline_surgery.txt
├── meta_clinical_timeline_treatment.txt
├── meta_cna.seg.txt
├── meta_linear_CNA.txt
├── meta_mutations_extended.txt
├── meta_rna_seq_v2_mrna.txt
├── meta_rna_seq_v2_mrna_median_Zscores.txt
├── meta_study.txt
└── meta_sv.txt
```
Note! Most other studies won't have a timeline set of files.

## Upload the final packages
Upload all of the directories named as study short names to `s3://kf-strides-232196027141-cbioportal-studies/studies/`. Authentication is done via **`aws sts assume-role`**. You can assume the appropriate role with the following command:
```sh
export $(printf "AWS_ACCESS_KEY_ID=%s AWS_SECRET_ACCESS_KEY=%s AWS_SESSION_TOKEN=%s" \
$(aws sts assume-role \
  --role-arn arn:aws:iam::232196027141:role/bix-dev-data-role \
  --role-session-name MySessionName \
  --query "Credentials.[AccessKeyId,SecretAccessKey,SessionToken]" \
  --output text))
```
Once authenticated, upload your processed folder(s) using:
```sh
aws s3 cp processed/ s3://kf-strides-232196027141-cbioportal-studies/studies/ --recursive
```
For more context on required permissions and roles, see the [AWS Infra PedCBioPortal Import](https://github.com/d3b-center/aws-infra-pedcbioportal-import) repo reference in the [Software Prerequisites](../README.md#software-prerequisites) section.

## Load into QA/PROD
An AWS step function exists to load studies on to the QA and PROD servers.
  + Create a branch in the `d3b-center/aws-infra-pedcbioportal-import` git repo (**MUST START WITH `feature/`**) and edit the `import_studies.txt` file with the study name you which to load. Can be an MSKCC datahub link or a local study name
  + Push the branch to remote - this will kick off a github action to load into QA
  + To load into prod, make a PR. On merge, load to prod will kick off
  + aws `stateMachinePedcbioImportservice` Step function service is used to view and manage running jobs
  + To repeat a load, click on the ▶️ icon in the git repo to select the job you want to re-run
  + *Note*, if your branch importStudies.txt is the same as main, you may have tot rigger it yourself. To do so, go to [actions](https://github.com/d3b-center/aws-infra-pedcbioportal-import/actions), on the left panel choose which action you want, then from the drop down in the right panel, pick which branch you want that action to run on 

### Congratulations, you did it!