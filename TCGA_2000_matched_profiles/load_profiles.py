#!/usr/bin/env python3
"""Load the paired TCGA matrices with identical sample and gene ordering.

Install pandas if needed. Place this script beside the extracted data files.
Numeric data are loaded as float32 to reduce memory; NA values stay missing.
The mutation matrix indicates observed records, not confirmed wild-type status.
"""
from pathlib import Path
import pandas as pd

base = Path(__file__).resolve().parent
samples = pd.read_csv(base / 'selected_samples.csv')
sample_ids = samples['sample_id'].tolist()
numeric_types = {sample_id: 'float32' for sample_id in sample_ids}
expression = pd.read_csv(base / 'expression_2000.tsv.gz', sep='\t', index_col=0, dtype=numeric_types)
mutation_presence = pd.read_csv(base / 'mutation_presence_2000.tsv.gz', sep='\t', index_col=0, dtype=numeric_types)

assert expression.columns.tolist() == sample_ids
assert mutation_presence.columns.tolist() == sample_ids
assert expression.index.equals(mutation_presence.index)
assert samples['patient_id'].nunique() == 2000

print('Expression, genes x samples:', expression.shape)
print('Mutation presence, genes x samples:', mutation_presence.shape)

# If your model expects samples as rows and genes as columns:
X_expression = expression.T
X_mutation = mutation_presence.T
print('Model input shape for either modality:', X_expression.shape)

# Full variant records can be read in chunks, retaining FILTER and consequence:
# for chunk in pd.read_csv(base / 'mutations_2000.maf.gz', sep='\t', comment='#',
#                          chunksize=100000, dtype=str, keep_default_na=False):
#     process(chunk)
