#!/usr/bin/env bash
# This script deletes a public virtual study from specified environment.
# Reference: https://docs.cbioportal.org/create-and-publish-virtual-study/#un-publish-virtual-study
set -euo pipefail
default_url="https://pedcbioportal.kidsfirstdrc.org"
read -p "Enter cBioportal URL to delete from. Leave Blank for $default_url: " url
url=${url:-$default_url}
echo "Removing from to $url"
# Prompt for the virtual study ID
read -p "Enter virtual study ID: " study_id
# Prompt for the publication token
read -s -p "Enter publication token: " token
curl \
  -X DELETE \
  -H "X-PUBLISHER-API-KEY: $token" \
  -v "$url/api/public_virtual_studies/$study_id"