# Compare current versus build
This documentation addresses a [QC script](../scripts/diff_studies.py) for clinical metadata. It streamlines the process of identifying and summarizing changes slated to be made.

```sh
python3 scripts/diff_studies.py --help
usage: diff_studies.py [-h] [-u URL] [-s STUDY] [-t TOKEN] [-d DATA_DIR]

Compare local clinical data to server

options:
  -h, --help            show this help message and exit
  -u URL, --url URL     url to search against
  -s STUDY, --study STUDY
                        Cancer study ID to compare on server
  -t TOKEN, --token TOKEN
                        Token file obtained from Web API
  -d DATA_DIR, --datasheet-dir DATA_DIR
                        Directory containing data_clinical_*.txt
```

## INPUTS:
 - `-u, --url`: cBioportal api deployment site. Default: https://pedcbioportal.kidsfirstdrc.org/api/v2/api-docs
 - `-s, --study`: cBioportal cancer study ID, i.e. `pbta_all`
 - `-t, --token`: File obtained from navigating to https://pedcbioportal.kidsfirstdrc.org/webAPI#using-data-access-tokens, then clicking on `Download Token`. File is reusable
 - `-d, --datasheet-dir`: Name of directory containing `data_clinical_patient.txt` and `data_clinical_sample.txt` being vetted for upload

## OUTPUTS:
Essentially two change log files, `patient_portal_v_build.txt` and `sample_portal_v_build.txt`.
For the patient and sample views, each file respectively has:
 - A list, one per line, per ID, per attribute, of what would change if the data were loaded
 - A list of IDs that would be removed from the portal, if any
 - A list of IDs that would be added if any
 - A summary of the number of changes of each attribute type

### patient_portal_v_build.txt example:
```
Per Patient changes:
Patient PT_017WC8PS attribute ETHNICITY would change from NA to Not Available
Patient PT_01HNFSBZ attribute ETHNICITY would change from NA to Not Available
Patient PT_01HNFSBZ attribute GERMLINE_SEX_ESTIMATE would change from Unknown to NA
Patient PT_01HNFSBZ attribute CANCER_PREDISPOSITIONS would change from None documented to NA
Patient PT_01SH4F1X attribute AGE_IN_DAYS would change from 3838 to NA
Patient PT_01SH4F1X attribute OS_MONTHS would change from 45 to 54
Patient PT_0324HWD5 attribute AGE_IN_DAYS would change from 3121 to NA
Patient PT_047YGDRW attribute ETHNICITY would change from NA to Not Available
Patient PT_04V47WFC attribute AGE_IN_DAYS would change from 5717 to NA
Patient PT_08M919BH attribute EFS_MONTHS would change from 44 to 62
Patient PT_08M919BH attribute OS_MONTHS would change from 44 to 62
Patient PT_0BSG3R3N attribute AGE_IN_DAYS would change from 3431 to NA
Patient PT_0BVR16FK attribute ETHNICITY would change from NA to Not Available
Patient PT_0CE0HFYB attribute GERMLINE_SEX_ESTIMATE would change from Male to NA

...

CHANGE SUMMARY:
ETHNICITY has 358 change(s)
GERMLINE_SEX_ESTIMATE has 220 change(s)
CANCER_PREDISPOSITIONS has 29 change(s)
AGE_IN_DAYS has 147 change(s)
OS_MONTHS has 99 change(s)
EFS_MONTHS has 60 change(s)
EFS_STATUS has 18 change(s)
SEX has 9 change(s)
AGE has 6 change(s)
OS_STATUS has 4 change(s)
```

### sample_portal_v_build.txt example:
```
Per Sample changes:
Sample 16510-1 attribute TUMOR_FRACTION would change from 0.349951221921 to 0.34995122192100003
Sample 16510-15 attribute TUMOR_FRACTION would change from 0.892871847605 to 0.8928718476049999
Sample 16510-2 attribute TUMOR_FRACTION would change from 0.242536563786 to 0.24253656378600005
Sample 16510-8 attribute TUMOR_FRACTION would change from 0.557284218924 to 0.5572842189239999
Sample 7316-100 attribute TUMOR_FRACTION would change from 0.270649989118 to 0.27064998911800003
Sample 7316-1017 attribute TUMOR_FRACTION would change from 0.570184695999 to 0.559801637737
Sample 7316-104 attribute TUMOR_FRACTION would change from 0.664343255194 to 0.6643432551940001
Sample 7316-1045 attribute TUMOR_FRACTION would change from 0.477859261757 to 0.496989582389
Sample 7316-105 attribute MOLECULAR_SUBTYPE would change from NA to LGG, BRAF V600E
Sample 7316-105 attribute TUMOR_PLOIDY would change from NA to 2
Sample 7316-105 attribute CANCER_TYPE_DETAILED would change from NA to Low-grade glioma, BRAF V600E
Sample 7316-105 attribute CANCER_GROUP would change from NA to Low-grade glioma
Sample 7316-105 attribute TUMOR_FRACTION would change from NA to 0.823344460708
Sample 7316-105 attribute PATHOLOGY_FREE_TEXT_DIAGNOSIS would change from NA to pilocytic astrocytoma ii

...

CHANGE SUMMARY:
27 Samples in build would be added to the portal: 1235928,1235929,1235930,1235931,1235932,1235933,1235934,1235935,1235936,1235937,1235938,1235939,1235940,1235941,1235981,1240110,1240112,1240114,1240116,1242273,1242274,1242276,1250775,1250776,1250777,1250778,1273223
TUMOR_FRACTION has 1005 change(s)
MOLECULAR_SUBTYPE has 488 change(s)
TUMOR_PLOIDY has 403 change(s)
CANCER_TYPE_DETAILED has 847 change(s)
CANCER_GROUP has 517 change(s)
PATHOLOGY_FREE_TEXT_DIAGNOSIS has 507 change(s)
BROAD_HISTOLOGY has 390 change(s)
CNS_REGION has 410 change(s)
CANCER_TYPE has 8 change(s)
ONCOTREE_CODE has 7 change(s)
TUMOR_TYPE has 8 change(s)
TUMOR_TISSUE_TYPE has 17 change(s)
EXPERIMENT_STRATEGY has 1 change(s)
SPECIMEN_ID has 1 change(s)
CBTN_TUMOR_TYPE has 6 change(s)
SAMPLE_TYPE has 2 change(s)
```