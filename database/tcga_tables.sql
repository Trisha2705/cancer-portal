
USE cancerinsight;

-- ==========================================
-- 1. TCGA RESEARCH SAMPLES
-- ==========================================

CREATE TABLE tcga_samples (
    sample_id VARCHAR(30) PRIMARY KEY,
    patient_id VARCHAR(30) NOT NULL,
    cancer_type VARCHAR(10) NOT NULL,

    sample_order_1based INT,
    same_sample_type_key BOOLEAN,
    exact_vial_match BOOLEAN,
    cancer_type_match BOOLEAN,

    expression_sample_id VARCHAR(30),
    expression_tsv_column_1based INT,
    patient_id_expression VARCHAR(30),
    sample_type_key_15 VARCHAR(20),
    vial_key_16_expression VARCHAR(20),
    sample_type_code_expression INT,
    is_tumour_sample_expression BOOLEAN,
    cancer_type_expression VARCHAR(10),
    cancer_type_source_expression VARCHAR(100),
    cancer_type_QC_expression VARCHAR(100),
    cancer_type_conflict_expression BOOLEAN,

    qc_aliquot_found_expression BOOLEAN,
    qc_Do_not_use_expression BOOLEAN,
    qc_AWG_excluded_because_of_pathology_expression VARCHAR(20),
    qc_patient_annotation_expression TEXT,
    qc_sample_annotation_expression VARCHAR(100),
    qc_aliquot_annotation_expression TEXT,
    qc_aliquot_annotation_updated_expression TEXT,

    mutation_sample_id VARCHAR(30),
    patient_id_mutation VARCHAR(30),
    vial_key_16_mutation VARCHAR(20),
    sample_type_code_mutation INT,
    is_tumour_sample_mutation BOOLEAN,
    cancer_type_mutation VARCHAR(10),
    cancer_type_source_mutation VARCHAR(100),
    cancer_type_QC_mutation VARCHAR(100),
    cancer_type_conflict_mutation BOOLEAN,

    qc_aliquot_found_mutation BOOLEAN,
    qc_Do_not_use_mutation VARCHAR(20),
    qc_AWG_excluded_because_of_pathology_mutation VARCHAR(20),
    qc_patient_annotation_mutation TEXT,
    qc_sample_annotation_mutation TEXT,
    qc_aliquot_annotation_mutation TEXT,
    qc_aliquot_annotation_updated_mutation TEXT,

    mutation_record_count INT,
    pass_mutation_record_count INT,
    has_PASS_record BOOLEAN,

    INDEX idx_tcga_patient (patient_id),
    INDEX idx_tcga_cancer_type (cancer_type)
);


-- ==========================================
-- 2. TCGA RESEARCH GENES
-- ==========================================

CREATE TABLE tcga_genes (
    gene_id VARCHAR(50) PRIMARY KEY,
    gene_symbol VARCHAR(100),
    entrez_gene_id BIGINT,
    mutation_matrix_row_mappable BOOLEAN,
    mutation_mapping_method VARCHAR(100),
    mutation_mapping_reason VARCHAR(100),
    observed_MAF_Ensembl_ids TEXT,

    INDEX idx_tcga_gene_symbol (gene_symbol),
    INDEX idx_tcga_entrez (entrez_gene_id)
);


USE cancerinsight;

CREATE TABLE tcga_variants (
    variant_id BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,

    sample_id VARCHAR(30) NOT NULL,
    patient_id VARCHAR(30),
    gene_symbol VARCHAR(100),
    entrez_gene_id VARCHAR(30),

    genome_assembly VARCHAR(20),
    chromosome VARCHAR(10),
    start_position BIGINT,
    end_position BIGINT,
    strand VARCHAR(10),

    variant_classification VARCHAR(100),
    variant_type VARCHAR(30),

    reference_allele TEXT,
    tumor_seq_allele1 TEXT,
    tumor_seq_allele2 TEXT,

    hgvsc TEXT,
    hgvsp TEXT,
    hgvsp_short TEXT,
    transcript_id TEXT,

    t_depth INT,
    t_ref_count INT,
    t_alt_count INT,
    n_depth INT,
    n_ref_count INT,
    n_alt_count INT,

    consequence TEXT,
    clin_sig TEXT,
    impact VARCHAR(50),
    filter_status VARCHAR(100),
    cosmic TEXT,

    INDEX idx_variant_sample (sample_id),
    INDEX idx_variant_patient (patient_id),
    INDEX idx_variant_gene (gene_symbol),
    INDEX idx_variant_classification (variant_classification),

    CONSTRAINT fk_tcga_variant_sample
        FOREIGN KEY (sample_id)
        REFERENCES tcga_samples(sample_id)
);
