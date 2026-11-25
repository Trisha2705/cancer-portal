TCGA: 2,000 RANDOMLY SELECTED MATCHED RNA/DNA PROFILES
Prepared 12 November 2026

QUICK START

Extract this ZIP into one folder. The data files are already complete; gzip
compression can be read directly by Python pandas or R data readers.
Both matrices have 20,531 gene rows and 2,000 sample columns, with identical
gene identifiers, sample names and ordering. sample_id is the shared 16-character
TCGA vial barcode. All 2,000 samples are from different patients.

FILES

expression_2000.tsv.gz
    20,531 x 2,000 expression matrix; first column gene_id.
    Original normalized expression strings and missing values are preserved.
    Column names are the common sample_id values, not the original RNA aliquot
    names. selected_samples.csv retains every original RNA/DNA identifier.

mutation_presence_2000.tsv.gz
    20,531 x 2,000 gene-level observed-mutation matrix, aligned exactly with
    expression_2000.tsv.gz by both gene rows and sample columns.
    1 = at least one MC3 source record maps to this gene/sample by an exact,
        unambiguous gene-symbol match.
    0 = no exact-symbol record is reported for this gene/sample.
    NA = a missing/ambiguous gene-symbol mapping; 60 expression
         gene rows are retained as NA rather than assigned mutation-negative.
    Zero is NOT proof of biological wild-type status or adequate sequencing
    coverage. This is an observed-record indicator, not a pathogenicity score.

mutations_2000.maf.gz
    Full per-variant profiles: 771,664 source MAF records for the
    selected 2,000 DNA aliquots. All 114 original source columns are retained,
    with sample_id and patient_id added at the beginning (116 columns total).
    Includes genomic position, allele, consequence, depths/counts, original
    tumour/normal barcodes, FILTER and other MC3 annotations. Genome build is
    GRCh37. This is the authoritative mutation detail file.

selected_samples.csv
    Exactly 2,000 samples, in the same order as both matrix headers. Includes
    common sample_id, patient_id, original expression_sample_id and
    mutation_sample_id, cancer_type, source quality flags and record counts.

mutation_sample_summary.csv
    Per-sample total/PASS/mapped mutation counts in matching sample order.

gene_mapping.csv
    Expression gene_id, source gene symbol, original positive Entrez identifier,
    and mutation-mapping eligibility/method/reason fields. Full gene_id values
    are preserved, including symbols represented as '?'.

cancer_type_counts.csv
    Eligible and selected patient counts per cancer type. All 33 types occur
    in the resulting random selection; their representation is not balanced.

selection_summary.json, mutation_subset_summary.json, validation_summary.json
    Sampling policy, observed counts, mapping/filter counts and validation.

load_profiles.py
    Example loader using pandas. It asserts alignment and demonstrates
    transposing to 2,000 samples x 20,531 genes for model input.

scripts/ and source_provenance.json
    Extraction code and provenance for the original official source files.
    SHA256SUMS.txt contains hashes for all delivered component files.

EXACT SELECTION METHOD

1. Eligibility: the supplied exact-vial crosswalk, containing 9,605 RNA-DNA
   pairs from 9,538 unique patients. Matching requires identical first 16
   barcode characters: patient, sample type and vial.
2. Sort the unique patient IDs; use Python random.Random(42).sample(..., 2000)
   without replacement. This is a simple random sample of patients, without
   stratification by cancer type. Random seed: 42.
3. For a selected patient with more than one eligible pair, choose the smallest
   (expression_sample_id, mutation_sample_id) pair lexicographically.
4. Sort the resulting 2,000 samples by the common sample_id and use that same
   order in both matrices and sample metadata.

The selection was independently re-derived from the original crosswalk and
matched exactly. Every selected RNA/DNA pair was verified at the vial level.

EXPRESSION VALUES

The source is the TCGA Pan-Cancer Atlas EB++ batch-corrected, upper-quartile
normalized RSEM expression resource. It is not raw read counts or TPM.
No extra log transformation, normalization, gene filtering or imputation was
performed. There are 687,292 missing expression entries in the subset,
preserved as missing. Missing expression is not zero expression.

MUTATION FILTERING AND INTERPRETATION

The MAF retains ALL source records for the selected aliquots: all consequence
classes, including synonymous/noncoding records, and all FILTER values.
Of 771,664 selected records, 666,183 have FILTER exactly PASS.
The observed-mutation matrix uses the same all-records definition.
It does not identify pathogenic or necessarily protein-altering variants.

The source MAF Entrez_Gene_Id field is literal 0 for all selected records,
so it cannot link variants to expression Entrez IDs. Instead, the matrix uses
exact, case-sensitive Hugo_Symbol -> expression gene_symbol matches where the
symbol uniquely identifies one expression row. '?' and duplicated expression
symbols are excluded from mapping. Symbols associated with multiple distinct
populated Ensembl Gene identifiers in the SELECTED MAF records are also treated
as ambiguous. The expression Entrez IDs are preserved for reference.

No alias guessing, silent case conversion, or merging of expression rows was
performed. A unique usable symbol with no observed exact-match records can
have zeros across all samples; its gene-mapping metadata indicates that the
mapping is unverified. Such zeros may also reflect annotation/name differences.
All selected variants remain in the full MAF, including any that do not map to
the expression gene list. Counts and mapping outcomes are recorded in
mutation_subset_summary.json and gene_mapping.csv.

This is an unfiltered random matched subset. Source QC annotations are retained
without applying additional study-specific exclusions. 24 selected RNA
samples carry Do_not_use=True. DNA aliquot QC metadata are unavailable for most
selected aliquots; a blank flag means unknown, not passing QC. There are 149
selected samples with source mutation records but no literal PASS records.
Choose further QC/consequence filtering explicitly for the intended analysis.

VALIDATION

- Exactly 2,000 unique patients and common sample IDs; one pair per patient.
- Both source aliquots in each pair share the same 16-character vial barcode.
- Both matrices have identical 20,531 gene rows and 2,000 sample columns.
- Every row has the expected number of values; mutation codes are 0/1 or NA
  according to each gene's mapping eligibility.
- All selected DNA aliquots are present in the source MAF.
- Every sample's total and PASS MAF counts match the original crosswalk.
- No expression values were numerically transformed; source NAs were retained.
- The complete source gzip files were reassembled with full SHA-256 checks;
  the original GDC downloads had passed their published MD5 checksums.

REPRODUCING EXTRACTION

With the original dataset and matching metadata arranged under SOURCE/raw
and SOURCE/processed, run from a working directory:

    python3 scripts/select_and_extract_expression.py --source-dir SOURCE --output-dir OUTPUT --n 2000 --seed 42

Then run scripts/extract_mutations.py with its --source-maf, --samples,
--gene-map and --output-dir arguments pointing to the source MAF and newly
written selected_samples.csv/gene_mapping.csv. The mutation script requires
NumPy; its --help describes the exact options. The expression script uses only
the Python standard library. The delivered data need no extraction rerun.

OFFICIAL SOURCES

TCGA Pan-Cancer Atlas supplemental datasets:
https://gdc.cancer.gov/about-data/publications/pancanatlas

GDC barcode definitions:
https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

Hoadley et al., Cell 2018, expression and Pan-Cancer Atlas methods:
https://doi.org/10.1016/j.cell.2018.03.022

Ellrott et al., Cell Systems 2018, MC3 somatic mutation methods:
https://doi.org/10.1016/j.cels.2018.03.002
