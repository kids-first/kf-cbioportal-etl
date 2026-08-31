#!/usr/bin/env python3
"""Merge rsem isoform outputs for cBio usage.

Converts to generic assay format
Replace BS ID with cBio ID
First column is RSEM isoform ID, second is gene ID and symbol, the rest are sample columns
Final result also z scored across cohort as log2(FPKM + 1) or as log2(TPM + 1)
"""

import argparse
import concurrent.futures
import io
import json
import os
import sys
import tarfile

import numpy as np
import pandas as pd
from scipy import stats

from cbioportal_etl.scripts.resolve_config_paths import resolve_config_paths