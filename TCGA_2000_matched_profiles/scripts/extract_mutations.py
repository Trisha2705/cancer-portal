#!/usr/bin/env python3
"""Extract TCGA MC3 mutation profiles for a fixed, paired sample manifest.

The MAF subset is the authoritative mutation output: all selected source rows,
consequence classes, FILTER flags and source fields are retained. The aligned
matrix is an observed-record indicator, not a pathogenicity or wild-type call.
Its columns follow selected_samples.csv and its rows follow gene_mapping.csv.
The delivered matrix uses exact, unambiguous Hugo symbols because this MC3
release populates Entrez_Gene_Id with zero. Unusable or ambiguous expression
symbols are NA; absent exact-symbol records are zero, with mapping status
documented. No gene-symbol aliases or case conversions are inferred.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
import gzip
import io
import json
from pathlib import Path
import re
import time

import numpy as np


def positive_entrez(value: str) -> int | None:
    value = str(value).strip()
    if not re.fullmatch(r"[0-9]+", value):
        return None
    number = int(value)
    return number if number > 0 else None


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


@contextmanager
def gzip_text_output(path: Path):
    """Deterministic gzip wrapper without changing input record text."""
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw,
                           compresslevel=3, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                yield text


def run(args) -> dict:
    samples = read_table(args.samples)
    genes = read_table(args.gene_map)
    if len(samples) != args.expected_samples:
        raise ValueError(f"Expected {args.expected_samples} samples, found {len(samples)}")
    if len(genes) != args.expected_genes:
        raise ValueError(f"Expected {args.expected_genes} gene rows, found {len(genes)}")
    if not samples or not genes:
        raise ValueError("Sample and gene manifests must both be nonempty")
    sample_required = {"sample_id", "patient_id", "expression_sample_id", "mutation_sample_id"}
    if not sample_required.issubset(samples[0]):
        raise ValueError("Sample manifest lacks required identifier columns")
    if not {"gene_id", "gene_symbol", "entrez_gene_id"}.issubset(genes[0]):
        raise ValueError("Gene map lacks required identifier columns")

    sample_ids = [row["sample_id"] for row in samples]
    mutation_ids = [row["mutation_sample_id"] for row in samples]
    patient_ids = [row["patient_id"] for row in samples]
    for label, values in [("sample", sample_ids), ("DNA aliquot", mutation_ids),
                          ("patient", patient_ids)]:
        if len(set(values)) != len(values):
            raise ValueError(f"Manifest has duplicate {label} identifiers")
    for row in samples:
        if row["expression_sample_id"][:16] != row["mutation_sample_id"][:16]:
            raise ValueError("A selected pair is not an exact-vial match")
        if row["sample_id"] != row["mutation_sample_id"][:16]:
            raise ValueError("Common sample_id does not equal the vial-level barcode")
        if row["patient_id"] != row["mutation_sample_id"][:12]:
            raise ValueError("Patient identifier disagrees with mutation barcode")
    dna_to_column = {barcode: column for column, barcode in enumerate(mutation_ids)}

    gene_ids = [row["gene_id"] for row in genes]
    gene_to_rows = defaultdict(list)
    valid_gene_row = np.zeros(len(genes), dtype=bool)
    for index, row in enumerate(genes):
        entrez = positive_entrez(row["entrez_gene_id"])
        if entrez is not None:
            gene_to_rows[entrez].append(index)
            valid_gene_row[index] = True
    presence = np.zeros((len(genes), len(samples)), dtype=np.uint8)

    # Each counter is indexed by the selected sample's final matrix column.
    total_by_sample = np.zeros(len(samples), dtype=np.int64)
    pass_by_sample = np.zeros(len(samples), dtype=np.int64)
    mapped_by_sample = np.zeros(len(samples), dtype=np.int64)
    invalid_entrez_by_sample = np.zeros(len(samples), dtype=np.int64)
    outside_expression_by_sample = np.zeros(len(samples), dtype=np.int64)
    filters, consequences, pass_consequences, builds = Counter(), Counter(), Counter(), Counter()
    observed_entrez, mapped_entrez, outside_entrez = set(), set(), set()
    source_rows = selected_rows = source_comments = 0
    start = last_update = time.monotonic()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    maf_out = args.output_dir / "mutations_2000.maf.gz"
    matrix_out = args.output_dir / "mutation_presence_2000.tsv.gz"
    header = None
    required_maf = {"Tumor_Sample_Barcode", "Entrez_Gene_Id", "Hugo_Symbol",
                    "FILTER", "Variant_Classification"}
    with gzip.open(args.source_maf, "rt", encoding="utf-8", newline="") as source:
        with gzip_text_output(maf_out) as output:
            for line in source:
                if line.startswith("#"):
                    source_comments += 1
                    output.write(line)
                    continue
                if not line.strip():
                    continue
                if header is None:
                    header = line.rstrip("\r\n").split("\t")
                    if not required_maf.issubset(header):
                        raise ValueError("Unexpected source MAF header")
                    if len(header) != len(set(header)):
                        raise ValueError("Duplicate MAF column names")
                    indices = {name: header.index(name) for name in required_maf}
                    build_index = header.index("NCBI_Build") if "NCBI_Build" in header else None
                    added_columns = [name for name in ("sample_id", "patient_id") if name not in header]
                    output.write("\t".join(added_columns + header) + "\n")
                    continue
                source_rows += 1
                # Only split the leading columns for unselected MAF records.
                tumor_index = indices["Tumor_Sample_Barcode"]
                leading = line.split("\t", tumor_index + 1)
                if len(leading) <= tumor_index:
                    raise ValueError(f"Malformed MAF record at data row {source_rows}")
                column = dna_to_column.get(leading[tumor_index])
                if column is None:
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != len(header):
                    raise ValueError(f"Selected MAF record has {len(fields)} fields, expected {len(header)}")
                row = samples[column]
                prefix = "\t".join(row[name] for name in added_columns)
                if prefix:
                    output.write(prefix + "\t")
                output.write(line)
                if not line.endswith(("\n", "\r")):
                    output.write("\n")

                selected_rows += 1
                total_by_sample[column] += 1
                filter_value = fields[indices["FILTER"]]
                consequence = fields[indices["Variant_Classification"]]
                filters[filter_value] += 1
                consequences[consequence] += 1
                if filter_value == "PASS":
                    pass_by_sample[column] += 1
                    pass_consequences[consequence] += 1
                if build_index is not None:
                    builds[fields[build_index]] += 1

                entrez = positive_entrez(fields[indices["Entrez_Gene_Id"]])
                if entrez is None:
                    invalid_entrez_by_sample[column] += 1
                else:
                    observed_entrez.add(entrez)
                    matching_rows = gene_to_rows.get(entrez)
                    if matching_rows:
                        mapped_by_sample[column] += 1
                        mapped_entrez.add(entrez)
                        presence[matching_rows, column] = 1
                    else:
                        outside_expression_by_sample[column] += 1
                        outside_entrez.add(entrez)
                if time.monotonic() - last_update >= 25:
                    print(json.dumps({"stage": "mutation_extract", "source_rows": source_rows,
                                      "selected_records": selected_rows}), flush=True)
                    last_update = time.monotonic()

    if header is None:
        raise ValueError("Source MAF has no header")
    missing = [samples[i]["mutation_sample_id"] for i in np.flatnonzero(total_by_sample == 0)]
    if missing:
        raise AssertionError(f"{len(missing)} selected DNA barcodes have no source MAF record: {missing[:5]}")
    if int(total_by_sample.sum()) != selected_rows:
        raise AssertionError("Per-sample record counts do not reconcile")
    if not np.array_equal(total_by_sample,
                          mapped_by_sample + invalid_entrez_by_sample + outside_expression_by_sample):
        raise AssertionError("Mapping categories fail to partition all selected records")
    if sum(filters.values()) != selected_rows or sum(consequences.values()) != selected_rows:
        raise AssertionError("Variant annotation counts do not reconcile")
    validated_manifest_counts = []
    for field, observed in [("mutation_record_count", total_by_sample),
                            ("pass_mutation_record_count", pass_by_sample)]:
        if field in samples[0]:
            expected = np.array([int(row[field]) for row in samples], dtype=np.int64)
            mismatches = np.flatnonzero(expected != observed)
            if len(mismatches):
                first = int(mismatches[0])
                raise AssertionError(
                    f"{field} differs from source crosswalk for {len(mismatches)} samples; "
                    f"first {mutation_ids[first]}: expected {expected[first]}, observed {observed[first]}")
            validated_manifest_counts.append(field)

    # Stream rows in the expression-gene order; retain invalid-ID rows as NA.
    with gzip_text_output(matrix_out) as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *sample_ids])
        na_values = ["NA"] * len(samples)
        for index, gene_id in enumerate(gene_ids):
            values = presence[index].tolist() if valid_gene_row[index] else na_values
            writer.writerow([gene_id, *values])

    affected_gene_rows = presence.sum(axis=0, dtype=np.int64)
    sample_summary_path = args.output_dir / "mutation_sample_summary.csv"
    summary_fields = ["sample_id", "patient_id", "mutation_sample_id", "mutation_record_count",
                      "PASS_mutation_record_count", "records_mapped_to_expression_genes",
                      "records_with_invalid_entrez", "records_with_entrez_outside_expression",
                      "expression_gene_rows_with_observed_mutation"]
    with sample_summary_path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=summary_fields)
        writer.writeheader()
        for index, row in enumerate(samples):
            writer.writerow({
                "sample_id": row["sample_id"], "patient_id": row["patient_id"],
                "mutation_sample_id": row["mutation_sample_id"],
                "mutation_record_count": int(total_by_sample[index]),
                "PASS_mutation_record_count": int(pass_by_sample[index]),
                "records_mapped_to_expression_genes": int(mapped_by_sample[index]),
                "records_with_invalid_entrez": int(invalid_entrez_by_sample[index]),
                "records_with_entrez_outside_expression": int(outside_expression_by_sample[index]),
                "expression_gene_rows_with_observed_mutation": int(affected_gene_rows[index]),
            })

    # Confirm the serialized matrix header, rather than only the input list.
    with gzip.open(matrix_out, "rt", encoding="utf-8", newline="") as stream:
        emitted_header = next(csv.reader(stream, delimiter="\t"))
    if emitted_header != ["gene_id", *sample_ids]:
        raise AssertionError("Serialized mutation matrix sample order differs from the manifest")

    summary = {
        "source_MAF": args.source_maf.name,
        "source_records_scanned": source_rows,
        "source_comment_lines_retained": source_comments,
        "source_MAF_column_count": len(header),
        "delivered_MAF_column_count": len(header) + len(added_columns),
        "added_MAF_columns": added_columns,
        "selected_samples": len(samples),
        "selected_unique_patients": len(set(patient_ids)),
        "selected_samples_with_observed_MAF_records": int((total_by_sample > 0).sum()),
        "selected_mutation_records": selected_rows,
        "selected_PASS_mutation_records": int(pass_by_sample.sum()),
        "selected_samples_with_PASS_records": int((pass_by_sample > 0).sum()),
        "selected_records_mapped_to_expression_genes": int(mapped_by_sample.sum()),
        "selected_records_with_invalid_entrez": int(invalid_entrez_by_sample.sum()),
        "selected_records_with_entrez_outside_expression": int(outside_expression_by_sample.sum()),
        "unique_valid_entrez_in_selected_records": len(observed_entrez),
        "unique_entrez_mapped_to_expression": len(mapped_entrez),
        "unique_entrez_outside_expression": len(outside_entrez),
        "matrix_gene_rows": len(genes),
        "matrix_sample_columns": len(samples),
        "matrix_gene_rows_with_valid_entrez": int(valid_gene_row.sum()),
        "matrix_gene_rows_set_to_NA": int((~valid_gene_row).sum()),
        "matrix_unique_valid_expression_entrez": len(gene_to_rows),
        "matrix_cells_equal_to_1": int(presence.sum()),
        "matrix_gene_rows_with_at_least_one_1": int((presence.sum(axis=1) > 0).sum()),
        "FILTER_counts": dict(sorted(filters.items())),
        "Variant_Classification_counts": dict(sorted(consequences.items())),
        "PASS_Variant_Classification_counts": dict(sorted(pass_consequences.items())),
        "NCBI_Build_counts": dict(sorted(builds.items())),
        "mapping_policy": "Exact positive integer Entrez_Gene_Id to expression entrez_gene_id; no symbol fallback",
        "mutation_matrix_definition": {
            "1": "At least one reported MC3 record for this sample maps to this gene",
            "0": "No mapped record observed in this release; not confirmed wild-type or proof of assay coverage",
            "NA": "Expression gene row lacks a valid positive Entrez identifier",
        },
        "variant_filtering": "None: all source FILTER values and all consequence classes retained",
        "record_preservation": "Original selected MAF fields and row order preserved; shared sample_id and patient_id prepended",
        "matrix_order": "Gene order follows gene_mapping.csv; sample order follows selected_samples.csv",
        "per_sample_source_crosswalk_counts_verified": validated_manifest_counts,
        "elapsed_seconds": round(time.monotonic() - start, 2),
    }
    (args.output_dir / "mutation_subset_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "raw_extraction_complete", "selected_samples": len(samples),
                      "selected_records": selected_rows,
                      "mapped_records": int(mapped_by_sample.sum()),
                      "elapsed_seconds": summary["elapsed_seconds"]}), flush=True)
    return summary


def remap_by_exact_symbols(args) -> dict:
    """Rebuild only the matrix/metadata; never rewrite the extracted MAF.

    Ambiguity checks use the selected MAF records: an expression symbol must
    identify exactly one expression row and must not be associated with more
    than one distinct valid ENSG identifier among the selected variant rows.
    """
    start = last_update = time.monotonic()
    samples = read_table(args.samples)
    genes = read_table(args.gene_map)
    if len(samples) != args.expected_samples or len(genes) != args.expected_genes:
        raise ValueError("Unexpected sample/gene manifest size during remapping")
    sample_ids = [row["sample_id"] for row in samples]
    mutation_ids = [row["mutation_sample_id"] for row in samples]
    if len(set(sample_ids)) != len(samples) or len(set(mutation_ids)) != len(samples):
        raise ValueError("Duplicate selected identifiers during remapping")
    dna_to_column = {name: index for index, name in enumerate(mutation_ids)}
    gene_ids = [row["gene_id"] for row in genes]
    original_gene_fields = list(genes[0])
    unusable_symbols = {"", "?", "NA", "NaN", "nan", ".", "-", "0"}

    def usable_symbol(value):
        return bool(value.strip()) and value not in unusable_symbols

    expression_symbol_rows = defaultdict(list)
    for index, row in enumerate(genes):
        if usable_symbol(row["gene_symbol"]):
            expression_symbol_rows[row["gene_symbol"]].append(index)
    unique_expression_symbols = {
        symbol: indices[0] for symbol, indices in expression_symbol_rows.items()
        if len(indices) == 1
    }
    duplicate_expression_symbols = {
        symbol for symbol, indices in expression_symbol_rows.items() if len(indices) > 1
    }

    observed_symbol_ensembl = defaultdict(set)
    observed_symbols = Counter()
    row_sample_counts = defaultdict(Counter)
    total = np.zeros(len(samples), dtype=np.int64)
    passes = np.zeros(len(samples), dtype=np.int64)
    invalid_symbol = np.zeros(len(samples), dtype=np.int64)
    absent_expression_symbol = np.zeros(len(samples), dtype=np.int64)
    duplicate_expression_symbol = np.zeros(len(samples), dtype=np.int64)
    invalid_entrez = np.zeros(len(samples), dtype=np.int64)
    entrez_values = Counter()
    rows_read = 0
    maf_path = args.output_dir / "mutations_2000.maf.gz"
    valid_ensembl = re.compile(r"ENSG[0-9]{11}(?:\.[0-9]+)?\Z")
    with gzip.open(maf_path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader((line for line in stream if not line.startswith("#")), delimiter="\t")
        required = {"sample_id", "patient_id", "Tumor_Sample_Barcode", "Hugo_Symbol",
                    "Entrez_Gene_Id", "Gene", "FILTER"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Extracted MAF lacks identifiers needed for exact-symbol remapping")
        for record in reader:
            column = dna_to_column.get(record["Tumor_Sample_Barcode"])
            if column is None:
                raise AssertionError("Extracted MAF contains a sample outside the manifest")
            if (record["sample_id"] != samples[column]["sample_id"] or
                    record["patient_id"] != samples[column]["patient_id"]):
                raise AssertionError("Extracted MAF common identifiers disagree with manifest")
            rows_read += 1
            total[column] += 1
            passes[column] += record["FILTER"] == "PASS"
            entrez_values[record["Entrez_Gene_Id"]] += 1
            invalid_entrez[column] += positive_entrez(record["Entrez_Gene_Id"]) is None

            symbol = record["Hugo_Symbol"]
            observed_symbols[symbol] += 1
            # Ensembl identifiers only: placeholders never create ambiguity.
            for token in re.split(r"[;,|&]", record["Gene"]):
                token = token.strip()
                if valid_ensembl.fullmatch(token):
                    # Different annotation versions of one stable ENSG ID
                    # refer to the same gene and do not create ambiguity.
                    observed_symbol_ensembl[symbol].add(token.split(".", 1)[0])
            if not usable_symbol(symbol):
                invalid_symbol[column] += 1
            elif symbol in unique_expression_symbols:
                row_sample_counts[unique_expression_symbols[symbol]][column] += 1
            elif symbol in duplicate_expression_symbols:
                duplicate_expression_symbol[column] += 1
            else:
                absent_expression_symbol[column] += 1
            if time.monotonic() - last_update >= 25:
                print(json.dumps({"stage": "exact_symbol_remapping", "records": rows_read}), flush=True)
                last_update = time.monotonic()

    if np.any(total == 0):
        raise AssertionError("Selected sample missing from extracted MAF during remapping")
    for field, observed in [("mutation_record_count", total), ("pass_mutation_record_count", passes)]:
        if field in samples[0]:
            expected = np.array([int(row[field]) for row in samples], dtype=np.int64)
            if not np.array_equal(expected, observed):
                raise AssertionError(f"{field} fails per-sample crosswalk reconciliation during remapping")

    ambiguous_ensembl_symbols = {
        symbol for symbol, identifiers in observed_symbol_ensembl.items()
        if len(identifiers) > 1
    }
    presence = np.zeros((len(genes), len(samples)), dtype=np.uint8)
    mappable = np.zeros(len(genes), dtype=bool)
    mapped = np.zeros(len(samples), dtype=np.int64)
    ambiguous_ensembl = np.zeros(len(samples), dtype=np.int64)
    reason_counts = Counter()
    for row_index, gene in enumerate(genes):
        symbol = gene["gene_symbol"]
        ensembl_ids = sorted(observed_symbol_ensembl.get(symbol, set()))
        if not usable_symbol(symbol):
            reason = "unusable_expression_symbol"
        elif symbol in duplicate_expression_symbols:
            reason = "symbol_identifies_multiple_expression_rows"
        elif symbol in ambiguous_ensembl_symbols:
            reason = "observed_symbol_has_multiple_Ensembl_gene_ids"
            for column, count in row_sample_counts.get(row_index, {}).items():
                ambiguous_ensembl[column] += count
        else:
            mappable[row_index] = True
            if observed_symbols.get(symbol, 0):
                reason = ("exact_unique_symbol_observed_single_Ensembl_id" if ensembl_ids
                          else "exact_unique_symbol_observed_no_valid_Ensembl_id")
            else:
                reason = "not_observed_mapping_unverified"
            for column, count in row_sample_counts.get(row_index, {}).items():
                mapped[column] += count
                presence[row_index, column] = 1
        gene["mutation_matrix_row_mappable"] = "True" if mappable[row_index] else "False"
        gene["mutation_mapping_method"] = "exact_Hugo_Symbol" if mappable[row_index] else "unmappable"
        gene["mutation_mapping_reason"] = reason
        gene["observed_MAF_Ensembl_ids"] = ";".join(ensembl_ids)
        reason_counts[reason] += 1

    if not np.array_equal(total, mapped + invalid_symbol + absent_expression_symbol +
                          duplicate_expression_symbol + ambiguous_ensembl):
        raise AssertionError("Exact-symbol mapping categories do not partition selected records")
    if int(presence.sum()) <= 0:
        raise AssertionError("Mutation matrix is all zero after remapping; do not deliver it")
    if np.any(presence[~mappable]):
        raise AssertionError("Unmappable expression rows contain mutation indicators")

    added_gene_fields = ["mutation_matrix_row_mappable", "mutation_mapping_method",
                         "mutation_mapping_reason", "observed_MAF_Ensembl_ids"]
    gene_fields = original_gene_fields + [x for x in added_gene_fields if x not in original_gene_fields]
    with args.gene_map.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=gene_fields)
        writer.writeheader()
        writer.writerows(genes)

    matrix_path = args.output_dir / "mutation_presence_2000.tsv.gz"
    with gzip_text_output(matrix_path) as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *sample_ids])
        missing = ["NA"] * len(samples)
        for row_index, gene_id in enumerate(gene_ids):
            writer.writerow([gene_id, *(presence[row_index].tolist() if mappable[row_index] else missing)])
    with gzip.open(matrix_path, "rt", encoding="utf-8", newline="") as stream:
        if next(csv.reader(stream, delimiter="\t")) != ["gene_id", *sample_ids]:
            raise AssertionError("Serialized matrix column order is incorrect")

    sample_fields = ["sample_id", "patient_id", "mutation_sample_id", "mutation_record_count",
                     "PASS_mutation_record_count", "records_mapped_to_expression_genes",
                     "records_with_invalid_MAF_symbol", "records_symbol_absent_from_expression",
                     "records_symbol_ambiguous_in_expression", "records_symbol_ambiguous_in_MAF_Ensembl",
                     "records_with_invalid_entrez", "expression_gene_rows_with_observed_mutation"]
    affected_rows = presence.sum(axis=0, dtype=np.int64)
    with (args.output_dir / "mutation_sample_summary.csv").open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=sample_fields)
        writer.writeheader()
        for column, sample in enumerate(samples):
            writer.writerow({
                "sample_id": sample["sample_id"], "patient_id": sample["patient_id"],
                "mutation_sample_id": sample["mutation_sample_id"],
                "mutation_record_count": int(total[column]),
                "PASS_mutation_record_count": int(passes[column]),
                "records_mapped_to_expression_genes": int(mapped[column]),
                "records_with_invalid_MAF_symbol": int(invalid_symbol[column]),
                "records_symbol_absent_from_expression": int(absent_expression_symbol[column]),
                "records_symbol_ambiguous_in_expression": int(duplicate_expression_symbol[column]),
                "records_symbol_ambiguous_in_MAF_Ensembl": int(ambiguous_ensembl[column]),
                "records_with_invalid_entrez": int(invalid_entrez[column]),
                "expression_gene_rows_with_observed_mutation": int(affected_rows[column]),
            })

    summary_path = args.output_dir / "mutation_subset_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    # Remove initial Entrez-only classification fields that no longer describe
    # the delivered symbol-mapped matrix. Preserve the original MAF statistics.
    for key in ["selected_records_with_entrez_outside_expression", "unique_entrez_mapped_to_expression",
                "unique_entrez_outside_expression"]:
        summary.pop(key, None)
    summary.update({
        "selected_mutation_records": rows_read,
        "selected_PASS_mutation_records": int(passes.sum()),
        "selected_records_mapped_to_expression_genes": int(mapped.sum()),
        "selected_records_with_invalid_entrez": int(invalid_entrez.sum()),
        "Entrez_Gene_Id_value_counts": dict(sorted(entrez_values.items())),
        "selected_records_with_invalid_MAF_symbol": int(invalid_symbol.sum()),
        "selected_records_symbol_absent_from_expression": int(absent_expression_symbol.sum()),
        "selected_records_symbol_ambiguous_in_expression": int(duplicate_expression_symbol.sum()),
        "selected_records_symbol_ambiguous_in_MAF_Ensembl": int(ambiguous_ensembl.sum()),
        "matrix_gene_rows_mappable_by_symbol": int(mappable.sum()),
        "matrix_gene_rows_set_to_NA": int((~mappable).sum()),
        "matrix_cells_equal_to_1": int(presence.sum()),
        "matrix_gene_rows_with_at_least_one_1": int((presence.sum(axis=1) > 0).sum()),
        "matrix_mapping_reason_counts": dict(sorted(reason_counts.items())),
        "ambiguous_expression_symbols": sorted(duplicate_expression_symbols),
        "observed_MAF_symbols_with_multiple_Ensembl_ids": {
            symbol: sorted(observed_symbol_ensembl[symbol]) for symbol in sorted(ambiguous_ensembl_symbols)
        },
        "mapping_policy": "Exact case-sensitive Hugo_Symbol to a unique expression gene_symbol; reject symbols with multiple observed valid ENSG IDs; no alias or case normalization",
        "mapping_ambiguity_check_scope": "Selected MAF records only",
        "Ensembl_ambiguity_identifiers": "Valid stable ENSG gene IDs; version suffixes are removed before distinct-ID counting",
        "Entrez_limitation": "All selected source Entrez_Gene_Id values are 0; those fields cannot link genes to the expression matrix",
        "mutation_matrix_definition": {
            "1": "At least one selected source MC3 record maps to this expression row by an unambiguous exact gene symbol",
            "0": "No exact-symbol mapped record observed in this release; not confirmed wild-type, assay coverage, or successful alias resolution",
            "NA": "Expression symbol is unusable, identifies multiple expression rows, or has multiple observed valid Ensembl gene IDs",
        },
        "unobserved_unique_expression_symbols": "Kept as 0 under the observed-record definition, with not_observed_mapping_unverified row metadata",
        "variant_filtering": "None: all source FILTER values and all consequence classes retained in both MAF and observed-record matrix",
        "symbol_remapping_elapsed_seconds": round(time.monotonic() - start, 2),
        "per_sample_source_crosswalk_counts_verified": ["mutation_record_count", "pass_mutation_record_count"],
    })
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "complete", "mapping": "exact_unambiguous_Hugo_Symbol",
                      "selected_samples": len(samples), "selected_records": rows_read,
                      "mapped_records": int(mapped.sum()), "positive_matrix_cells": int(presence.sum()),
                      "NA_gene_rows": int((~mappable).sum()),
                      "remapping_seconds": summary["symbol_remapping_elapsed_seconds"]}), flush=True)
    return summary


def main():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-maf", type=Path,
                        default=base / "source/raw/mc3.v0.2.8.PUBLIC.maf.gz")
    parser.add_argument("--samples", type=Path, default=base / "output/selected_samples.csv")
    parser.add_argument("--gene-map", type=Path, default=base / "output/gene_mapping.csv")
    parser.add_argument("--output-dir", type=Path, default=base / "output")
    parser.add_argument("--expected-samples", type=int, default=2000)
    parser.add_argument("--expected-genes", type=int, default=20531)
    parser.add_argument("--remap-existing", action="store_true",
                        help="Keep an already extracted MAF unchanged and rebuild its symbol-mapped matrix/metadata")
    args = parser.parse_args()
    if not args.remap_existing:
        run(args)
    remap_by_exact_symbols(args)


if __name__ == "__main__":
    main()
