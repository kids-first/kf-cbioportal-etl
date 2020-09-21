# New Study Load Request

## Description
Please provide the following information:
1. Is this a KF or external study?

1. If KF, please provide the study ID (`SD_122345`) and study title

1. Please fill out out the cbio meta_study info, for instance:
    ```
    type_of_cancer: brain
    cancer_study_identifier: pbta_all
    name: Pediatric Brain Tumor Atlas (PBTA, Provisional)
    description: Genomic characterization of Pediatric Brain Tumor Atlas samples provided by the <a href="http://CBTTC.org">Children's Brain Tumor Tissue Consortium</a>, the <a href="http://www.pnoc.us/">Pacific Neuroncology Consortium (PNOC)</a>, and its partners via the <a href="http://kidsfirstdrc.org">Gabriella Miller Kids First Data Resource Center</a>. Updated Februrary 2020 from last load, July 2019.
    short_name: pbta_all
    ```
    Type of cancer can be found here: https://github.com/kids-first/kf-cbioportal-etl/blob/master/REFS/type_of_cancer
1. If KF, provide the Cavatica link(s) with delivered or harmonized data. Ops can help if needed

1. Provide an example of a sample ID that can be used to tie together DNA and RNA (if applicable)

1. Indicate whether samples are to be loaded into one big study, or broken down into sub-studies


## Data to load
Please check off which data types are to be loaded:

- [ ] SNV
- [ ] CNV
- [ ] RNA expression
- [ ] RNA fusions

## Study load status
Is this study `Provisional` - a KF dataset that has been freshly haromonized with little filtration, or other/`Publication`?

## Sample inclusion/exclusion
Please list any samples that should be excluded and why