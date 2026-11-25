
import pandas as pd
import mysql.connector
from pathlib import Path
import time

# Locate the original compressed MAF file
BASE_DIR = Path(__file__).resolve().parent.parent
MAF_FILE = BASE_DIR / "TCGA_2000_matched_profiles" / "mutations_2000.maf.gz"

BATCH_SIZE = 5000

# Connect to MySQL
conn = mysql.connector.connect(
    host="localhost",
    user="root",
    password=input("MySQL password: "),
    database="cancerinsight"
)

cursor = conn.cursor()

# Prevent accidentally importing the same records twice
cursor.execute("SELECT COUNT(*) FROM tcga_variants")
existing_count = cursor.fetchone()[0]

if existing_count > 0:
    print(f"Import stopped: tcga_variants already contains {existing_count} records.")
    print("Do not rerun without checking the existing data first.")
    cursor.close()
    conn.close()
    raise SystemExit

# Original MAF column names mapped to our SQL table
columns = [
    "sample_id",
    "patient_id",
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "HGVSc",
    "HGVSp",
    "HGVSp_Short",
    "Transcript_ID",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
    "n_ref_count",
    "n_alt_count",
    "Consequence",
    "CLIN_SIG",
    "IMPACT",
    "FILTER",
    "COSMIC"
]

# SQL table column names, in matching order
insert_sql = """
INSERT INTO tcga_variants (
    sample_id, patient_id, gene_symbol, entrez_gene_id,
    genome_assembly, chromosome, start_position, end_position,
    strand, variant_classification, variant_type,
    reference_allele, tumor_seq_allele1, tumor_seq_allele2,
    hgvsc, hgvsp, hgvsp_short, transcript_id,
    t_depth, t_ref_count, t_alt_count,
    n_depth, n_ref_count, n_alt_count,
    consequence, clin_sig, impact, filter_status, cosmic
)
VALUES (
    %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
    %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
    %s, %s, %s, %s, %s, %s, %s, %s, %s
)
"""

def clean_text(value):
    if pd.isna(value):
        return None
    return str(value)

def clean_integer(value):
    if pd.isna(value):
        return None
    try:
        return int(value)
    except (ValueError, TypeError):
        return None

# Numeric fields in the source MAF
integer_columns = {
    "Start_Position",
    "End_Position",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
    "n_ref_count",
    "n_alt_count"
}

print("Starting TCGA mutation import...")
print(f"Source file: {MAF_FILE}")
print(f"Batch size: {BATCH_SIZE:,}")

start_time = time.time()
total_imported = 0

try:
    # Read the compressed file in manageable batches
    reader = pd.read_csv(
        MAF_FILE,
        sep="\t",
        compression="gzip",
        dtype=str,
        chunksize=BATCH_SIZE,
        na_values=["."]
    )

    for chunk_number, chunk in enumerate(reader, start=1):

        # Check that all expected source columns exist
        missing_columns = [col for col in columns if col not in chunk.columns]

        if missing_columns:
            raise ValueError(f"Missing MAF columns: {missing_columns}")

        # Select only the columns we need
        selected = chunk[columns]

        rows = []

        for record in selected.itertuples(index=False, name=None):
            cleaned_record = []

            for column, value in zip(columns, record):
                if column in integer_columns:
                    cleaned_record.append(clean_integer(value))
                else:
                    cleaned_record.append(clean_text(value))

            rows.append(tuple(cleaned_record))

        # Insert this batch into MySQL
        cursor.executemany(insert_sql, rows)
        conn.commit()

        total_imported += len(rows)

        print(
            f"Batch {chunk_number}: "
            f"{total_imported:,} records imported"
        )

    elapsed = time.time() - start_time

    cursor.execute("SELECT COUNT(*) FROM tcga_variants")
    database_count = cursor.fetchone()[0]

    print("\nImport complete!")
    print(f"Records imported: {total_imported:,}")
    print(f"Records in database: {database_count:,}")
    print(f"Time taken: {elapsed:.1f} seconds")

except Exception:
    conn.rollback()
    print("\nAn error occurred during import.")
    print("Previously committed batches may still be in the database.")
    raise

finally:
    cursor.close()
    conn.close()
