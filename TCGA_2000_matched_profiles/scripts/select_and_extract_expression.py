#!/usr/bin/env python3
"""Select a reproducible matched TCGA subset and extract expression unchanged.

Requires Python's standard library only. Sampling is uniform over unique
patients with exact-vial RNA/DNA matches, not stratified by cancer type.
One pair per patient is chosen by lexicographic full RNA then DNA barcode.
Source QC flags are retained without imposing additional eligibility filters.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import random
import sys
import time

EXPRESSION_NAME = 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'

def open_text(path):
    return gzip.open(path, 'rt', newline='', encoding='utf-8-sig') if path.suffix == '.gz' else path.open('r', newline='', encoding='utf-8-sig')

def select_samples(source, output, n, seed):
    path = source / 'processed/matched_tumour_samples_exact_vial.csv'
    with path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        original_fields = reader.fieldnames
        rows = list(reader)
    if not rows:
        raise ValueError('The exact-vial crosswalk is empty')
    groups = defaultdict(list)
    for row in rows:
        rna, dna, patient = row['expression_sample_id'], row['mutation_sample_id'], row['patient_id']
        if rna[:16] != dna[:16] or rna[:12] != patient or dna[:12] != patient:
            raise ValueError('Crosswalk contains an invalid exact-vial pair')
        groups[patient].append(row)
    patients = sorted(groups)
    if n > len(patients) or n < 1:
        raise ValueError(f'Requested {n} patients from {len(patients)} eligible patients')
    drawn = random.Random(seed).sample(patients, n)
    selected = [min(groups[p], key=lambda r: (r['expression_sample_id'], r['mutation_sample_id'])) for p in drawn]
    selected.sort(key=lambda r: (r['expression_sample_id'][:16], r['expression_sample_id'], r['mutation_sample_id']))
    selected = [{'sample_order_1based': i + 1, 'sample_id': row['expression_sample_id'][:16], **row} for i, row in enumerate(selected)]
    assert len({r['patient_id'] for r in selected}) == n
    assert len({r['sample_id'] for r in selected}) == n
    assert len({r['expression_sample_id'] for r in selected}) == n
    assert len({r['mutation_sample_id'] for r in selected}) == n
    with (output / 'selected_samples.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['sample_order_1based', 'sample_id', *original_fields])
        writer.writeheader()
        writer.writerows(selected)
    eligible_cancers = Counter(groups[p][0]['cancer_type'] for p in patients)
    selected_cancers = Counter(r['cancer_type'] for r in selected)
    with (output / 'cancer_type_counts.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['cancer_type', 'eligible_patients', 'selected_patients'])
        writer.writerows((c, eligible_cancers[c], selected_cancers[c]) for c in sorted(eligible_cancers))
    summary = {
        'seed': seed,
        'sampling_algorithm': 'random.Random(seed).sample(sorted(unique_patient_ids), n)',
        'sampling_unit': 'unique patient',
        'cancer_stratification': False,
        'eligible_exact_vial_pairs': len(rows),
        'eligible_unique_patients': len(patients),
        'selected_samples': n,
        'selected_unique_patients': n,
        'selected_cancer_type_counts': dict(sorted(selected_cancers.items())),
        'selected_cancer_types': len(selected_cancers),
        'representative_pair_rule': 'Minimum full expression_sample_id, then mutation_sample_id, within patient',
        'matrix_column_order': 'Ascending shared 16-character sample_id; same order as selected_samples.csv',
        'source_QC_exclusions_applied': False,
        'selected_RNA_Do_not_use_true': sum(r['qc_Do_not_use_expression'].lower() == 'true' for r in selected),
        'selected_DNA_Do_not_use_true': sum(r['qc_Do_not_use_mutation'].lower() == 'true' for r in selected),
        'selected_RNA_QC_aliquot_found': sum(r['qc_aliquot_found_expression'].lower() == 'true' for r in selected),
        'selected_DNA_QC_aliquot_found': sum(r['qc_aliquot_found_mutation'].lower() == 'true' for r in selected),
        'expected_mutation_records_from_crosswalk': sum(int(r['mutation_record_count']) for r in selected),
        'expected_PASS_mutation_records_from_crosswalk': sum(int(r['pass_mutation_record_count']) for r in selected),
        'selected_samples_with_no_PASS_records': sum(int(r['pass_mutation_record_count']) == 0 for r in selected),
        'python_version': sys.version,
        'expression_transformations': 'Column subsetting and renaming only; original numeric and missing-value strings preserved',
    }
    return selected, summary

def extract_expression(source, output, selected, summary):
    path = source / 'raw' / EXPRESSION_NAME
    if not path.is_file():
        path = path.with_name(path.name + '.gz')
    gene_count = 0
    missing_count = 0
    duplicate_gene_ids = set()
    genes_seen = set()
    entrez_counts = Counter()
    last_progress = time.monotonic()
    with open_text(path) as f, gzip.open(output / 'expression_2000.tsv.gz', 'wt', newline='', compresslevel=3) as expr, (output / 'gene_mapping.csv').open('w', newline='') as gm:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        if header[0] != 'gene_id' or len(header) != len(set(header)):
            raise ValueError('Unexpected expression source header')
        positions = {name: i for i, name in enumerate(header)}
        selected_positions = [positions[r['expression_sample_id']] for r in selected]
        writer = csv.writer(expr, delimiter='\t', lineterminator='\n')
        writer.writerow(['gene_id', *[r['sample_id'] for r in selected]])
        gw = csv.writer(gm)
        gw.writerow(['gene_id', 'gene_symbol', 'entrez_gene_id'])
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f'Incorrect column count in expression row {gene_count + 1}')
            gene_id = row[0]
            if gene_id in genes_seen:
                duplicate_gene_ids.add(gene_id)
            genes_seen.add(gene_id)
            symbol, sep, raw_entrez = gene_id.rpartition('|')
            if not sep:
                symbol, raw_entrez = gene_id, ''
            entrez = str(int(raw_entrez)) if raw_entrez.isdigit() and int(raw_entrez) > 0 else ''
            if entrez:
                entrez_counts[entrez] += 1
            gw.writerow([gene_id, symbol, entrez])
            values = [row[i] for i in selected_positions]
            missing_count += sum(v.strip().lower() in {'na', 'nan', 'null', ''} for v in values)
            writer.writerow([gene_id, *values])
            gene_count += 1
            if time.monotonic() - last_progress > 25:
                print(json.dumps({'stage': 'expression', 'gene_rows': gene_count}), flush=True)
                last_progress = time.monotonic()
    summary.update({
        'expression_gene_rows': gene_count,
        'expression_samples': len(selected),
        'expression_values': gene_count * len(selected),
        'expression_missing_values_preserved': missing_count,
        'expression_duplicate_full_gene_ids': sorted(duplicate_gene_ids),
        'expression_unique_positive_entrez_ids': len(entrez_counts),
        'expression_rows_without_positive_entrez_id': gene_count - sum(entrez_counts.values()),
        'expression_positive_entrez_ids_with_multiple_rows': {k: v for k, v in entrez_counts.items() if v > 1},
    })
    if gene_count != 20531:
        raise ValueError(f'Expected 20531 source expression rows; found {gene_count}')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, default=Path(__file__).resolve().parent / 'source')
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent / 'output')
    parser.add_argument('--n', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected, summary = select_samples(args.source_dir, args.output_dir, args.n, args.seed)
    print(json.dumps({'stage': 'selected', 'patients': len(selected), 'cancer_types': summary['selected_cancer_types'], 'seed': args.seed}), flush=True)
    extract_expression(args.source_dir, args.output_dir, selected, summary)
    (args.output_dir / 'selection_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)

if __name__ == '__main__':
    main()
