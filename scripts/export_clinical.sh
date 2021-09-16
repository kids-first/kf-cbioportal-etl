#!/bin/bash

#############################################################
# This script is more of a guide on how to export all the   #
# necessary clinical files and manifests from D3b Warehouse #
#############################################################

# get data clinical sample file
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from pbta_all_data_clinical_sample) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > data_clinical_sample.txt
# get data clinical patient file
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from pbta_all_data_clinical_patient) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > data_clinical_patient.txt
# get genomics metadata file: often the -t arg input
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from cbio_genomics_file_etl_dict) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > cbio_file_name_id.txt
# get sequencing center info resource file
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from seq_center_resource) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > seq_center_bs_ids.txt
# copy this into the load package folder
psql -h d3b-warehouse-aurora-prd.d3b.io -U brownm28 postgres -c "SET search_path TO bix_workflows; COPY(select * from dgd_data_gene_matrix) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > data_gene_matrix_CHOP.txt