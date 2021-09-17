#!/bin/bash

#############################################################
# This script is more of a guide on how to export all the   #
# necessary clinical files and manifests from D3b Warehouse #
#############################################################

# Genomic files manifest
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from genomics_file_manifest) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > genomics_file_manifest.txt
# Genomic files manifest
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from dgd_genomics_file_manifest) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > dgd_genomics_file_manifest.txt
# Genomics metadata file: often the -t arg input
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from cbio_genomics_file_etl_dict) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > cbio_file_name_id.txt
# Data clinical sample sheet
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from pbta_all_data_clinical_sample) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > datasheets/pbta_all/data_clinical_sample.txt
# Data clinical patient sheet
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from pbta_all_data_clinical_patient) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > datasheets/pbta_all/data_clinical_patient.txt
# Sequencing center info resource file
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from seq_center_resource) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > seq_center_bs_ids.txt
# Data gene matrix
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from dgd_data_gene_matrix) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > data_gene_matrix_CHOP.txt