
import pandas as pd
import mysql.connector
import numpy as np

# Connect to CancerInsight
connection = mysql.connector.connect(
    host="localhost",
    user="root",
    password=input("MySQL password: "),
    database="cancerinsight"
)

cursor = connection.cursor()

# Dataset location
DATA_DIR = "../TCGA_2000_matched_profiles/"

# Convert pandas missing values into SQL NULL
def clean_value(value):
    if pd.isna(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


# ==========================================
# IMPORT TCGA SAMPLES
# ==========================================

print("Loading sample metadata...")

samples = pd.read_csv(DATA_DIR + "selected_samples.csv")

sample_columns = [
    "sample_id", "patient_id", "cancer_type",
    "sample_order_1based", "same_sample_type_key",
    "exact_vial_match", "cancer_type_match",
    "expression_sample_id", "expression_tsv_column_1based",
    "patient_id_expression", "sample_type_key_15",
    "vial_key_16_expression", "sample_type_code_expression",
    "is_tumour_sample_expression", "cancer_type_expression",
    "cancer_type_source_expression", "cancer_type_QC_expression",
    "cancer_type_conflict_expression",
    "qc_aliquot_found_expression", "qc_Do_not_use_expression",
    "qc_AWG_excluded_because_of_pathology_expression",
    "qc_patient_annotation_expression",
    "qc_sample_annotation_expression",
    "qc_aliquot_annotation_expression",
    "qc_aliquot_annotation_updated_expression",
    "mutation_sample_id", "patient_id_mutation",
    "vial_key_16_mutation", "sample_type_code_mutation",
    "is_tumour_sample_mutation", "cancer_type_mutation",
    "cancer_type_source_mutation", "cancer_type_QC_mutation",
    "cancer_type_conflict_mutation", "qc_aliquot_found_mutation",
    "qc_Do_not_use_mutation",
    "qc_AWG_excluded_because_of_pathology_mutation",
    "qc_patient_annotation_mutation",
    "qc_sample_annotation_mutation",
    "qc_aliquot_annotation_mutation",
    "qc_aliquot_annotation_updated_mutation",
    "mutation_record_count", "pass_mutation_record_count",
    "has_PASS_record"
]

sample_sql = """
INSERT INTO tcga_samples (
    sample_id, patient_id, cancer_type,
    sample_order_1based, same_sample_type_key,
    exact_vial_match, cancer_type_match,
    expression_sample_id, expression_tsv_column_1based,
    patient_id_expression, sample_type_key_15,
    vial_key_16_expression, sample_type_code_expression,
    is_tumour_sample_expression, cancer_type_expression,
    cancer_type_source_expression, cancer_type_QC_expression,
    cancer_type_conflict_expression,
    qc_aliquot_found_expression, qc_Do_not_use_expression,
    qc_AWG_excluded_because_of_pathology_expression,
    qc_patient_annotation_expression,
    qc_sample_annotation_expression,
    qc_aliquot_annotation_expression,
    qc_aliquot_annotation_updated_expression,
    mutation_sample_id, patient_id_mutation,
    vial_key_16_mutation, sample_type_code_mutation,
    is_tumour_sample_mutation, cancer_type_mutation,
    cancer_type_source_mutation, cancer_type_QC_mutation,
    cancer_type_conflict_mutation, qc_aliquot_found_mutation,
    qc_Do_not_use_mutation,
    qc_AWG_excluded_because_of_pathology_mutation,
    qc_patient_annotation_mutation,
    qc_sample_annotation_mutation,
    qc_aliquot_annotation_mutation,
    qc_aliquot_annotation_updated_mutation,
    mutation_record_count, pass_mutation_record_count,
    has_PASS_record
) VALUES (""" + ",".join(["%s"] * len(sample_columns)) + ")"

sample_values = [
    tuple(clean_value(v) for v in row)
    for row in samples[sample_columns].itertuples(index=False, name=None)
]

cursor.executemany(sample_sql, sample_values)
connection.commit()

print("Imported samples:", len(samples))


# ==========================================
# IMPORT TCGA GENES
# ==========================================

print("Loading gene metadata...")

genes = pd.read_csv(DATA_DIR + "gene_mapping.csv")

gene_columns = [
    "gene_id",
    "gene_symbol",
    "entrez_gene_id",
    "mutation_matrix_row_mappable",
    "mutation_mapping_method",
    "mutation_mapping_reason",
    "observed_MAF_Ensembl_ids"
]

gene_sql = """
INSERT INTO tcga_genes (
    gene_id, gene_symbol, entrez_gene_id,
    mutation_matrix_row_mappable,
    mutation_mapping_method, mutation_mapping_reason,
    observed_MAF_Ensembl_ids
) VALUES (%s, %s, %s, %s, %s, %s, %s)
"""

gene_values = [
    tuple(clean_value(v) for v in row)
    for row in genes[gene_columns].itertuples(index=False, name=None)
]

cursor.executemany(gene_sql, gene_values)

connection.commit()
print("Imported genes:", len(genes))


# ==========================================
# VERIFY IMPORTS
# ==========================================

cursor.execute("SELECT COUNT(*) FROM tcga_samples")
print("Samples in database:", cursor.fetchone()[0])

cursor.execute("SELECT COUNT(*) FROM tcga_genes")
print("Genes in database:", cursor.fetchone()[0])

cursor.close()
connection.close()

print("TCGA metadata import complete!")
