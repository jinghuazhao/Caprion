# Caprion MS-based Proteomics Pipeline

(Experimental)

**Location:** `~/Caprion/analysis/pipeline/`

This directory contains scripts for preparing proteomics inputs, configuring and running FragPipe, validating search results, and generating identification and quality-control reports.

The numbered prefixes indicate workflow stages approximately; they do not always establish execution order. Follow the dependencies and instructions below rather than relying on alphabetical sorting alone.

## File list

| File | Purpose | Main output or action |
|---|---|---|
| `01_inventory.sh` | Inventory RAW and mzML files; identify missing matches and duplicates | Inventory report |
| `02_convert_missing.sh` | Convert RAW files without corresponding mzML files | Missing mzML files |
| `03_validate_mzml.sh` | Validate mzML files and report empty, invalid or unmatched files | mzML validation report |
| `04_prepare_database.sh` | Check the target/decoy FASTA database and its integrity | Database validation |
| `05_prepare_workflow.sh` | Prepare a separate FragPipe workflow copy | Workflow file |
| `05b_set_workflow_database.sh` | Set the intended database in the prepared workflow | Updated workflow copy |
| `06_build_manifest.sh` | Build and validate the FragPipe manifest using explicit experiment and fraction assignments | Manifest TSV |
| `07_dry_run_fragpipe.sh` | Check configuration and perform a FragPipe dry run | Dry-run log |
| `08_review_dry_run.sh` | Review dry-run configuration and diagnostic output | Review findings |
| `08b_review_search_parameters.sh` | Review search parameters before production execution | Parameter review findings |
| `09_run_fragpipe_production.sh` | Run the production FragPipe search in an isolated output directory | Search results and production log |
| `monitor_fragpipe_production.sh` | Monitor an existing production run | Monitoring output |
| `09_validate_results.sh` | Validate result files, columns and identification-related outputs | Results validation report |
| `10_summarise_identifications.sh` | Summarise PSMs, peptides, ions and proteins | Identification summary |
| `11_report_sample_qc.sh` | Summarise identification and QC metrics by file | Sample QC report |
| `12_summarise_batch.sh` | Combine available inventory, validation, identification and QC summaries | Batch summary |

## Intended workflow

1. Inventory the input files.
2. Convert missing mzML files and validate the resulting inputs. Revalidate after conversion when necessary.
3. Validate the target/decoy FASTA database.
4. Prepare the workflow and configure its database.
5. Build and validate the manifest.
6. Perform the FragPipe dry run and review the configuration and search parameters.
7. Run the production search using `09_run_fragpipe_production.sh`.
8. Use `monitor_fragpipe_production.sh` as needed while the run is in progress.
9. Validate the completed results.
10. Generate the identification summary, sample QC report and batch summary.

Confirm the dependencies of the workflow-configuration and review scripts before treating this as a strict execution recipe.

## Working and output directories

| Directory | Purpose |
|---|---|
| `dry_run_98_v2/` | Dry-run output and diagnostics |
| `fragpipe_test_home/` | Isolated FragPipe test environment |
| `production_runs/` | Production search outputs, separated by run |
| `reports/` | Validation, identification, QC and batch reports |

These directories contain generated or operational data. They are excluded from version control through `.gitignore` where appropriate. Excluding them from Git does not replace the need for backups and data retention.

## Design principles

- Keep the original scripts and successful six-file pilot unchanged.
- Use a separate output directory for every new search run.
- Never silently replace the validated target/decoy database.
- Do not modify original FragPipe workflows in place.
- Do not infer biological sample identities from filenames.
- Require explicit experiment and fraction assignments before production runs.
- Stop on validation failures rather than continuing with incomplete inputs.
- Preserve production logs and reports for reproducibility.
- Treat modification timestamps as file-history evidence, not definitive execution timestamps.
- Do not change search parameters or regenerate a production database without an explicit, documented decision.

## Reference database and workflow provenance

The production pipeline uses reference assets prepared under:

`~/Caprion/analysis/ZWK/discovery_fresh/`

These assets are separate from the numbered pipeline scripts and should be preserved independently.

### Source protein FASTA

The source FASTA used for target/decoy generation is supplied with FragPipe 22.0:

`/usr/local/Cluster-Apps/ceuadmin/FragPipe/22.0/databases/UP000005640-canonical-proteins-20420.fasta`

This is the source protein database, not the generated production target/decoy FASTA.

### Generated target/decoy database

The generated production FASTA is:

`~/Caprion/analysis/ZWK/discovery_fresh/database_canonical/2026-10-09-decoys-UP000005640-canonical-proteins-20420.fasta.fas`

Philosopher 5.1.1 was used to generate the target/decoy database with the decoy prefix `rev_`.

The validated production database contains:

- Target sequences: 20,420
- Decoy sequences: 20,420
- Total sequences: 40,840

The decoy prefix is used by the downstream identification-filtering configuration.

`04_prepare_database.sh` checks the database as part of the pipeline. The database-generation procedure below documents the separate operation that created the target/decoy FASTA.

Do not replace or regenerate the production database silently. Any new database should be created as a separate version, validated, and documented before use.

### FragPipe workflow

The reference workflow is:

`~/Caprion/analysis/ZWK/discovery_fresh/LFQ-MBR-Sage-comparison.workflow`

The pipeline prepares a separate workflow copy and configures its database through `05_prepare_workflow.sh` and `05b_set_workflow_database.sh`. Keep the original reference workflow unchanged.

### Production provenance

For each production run, record:

- The exact source and generated FASTA paths.
- The source FASTA checksum and generated FASTA checksum.
- The reference workflow and prepared workflow copy.
- The manifest and sample metadata version used.
- The FragPipe and Philosopher versions.
- The relevant search parameters.
- The production output directory and log.
- The database validation results.

Use the production log and run-specific configuration to verify these details. Modification timestamps alone do not establish which assets were used.

For the completed 98-file production run, confirm the manifest version from the production log before recording `production_v2.fp-manifest` as the manifest actually used.

## Re-create the target/decoy database

**Warning:** This is a database-generation operation, not a routine pipeline step. Do not run it over the validated production database without first preserving that database and its checksum. Prefer a separate output directory when creating a new version.

The commands below document the procedure used to generate the target/decoy FASTA from the source FASTA bundled with FragPipe 22.0.

```bash
cd ~/Caprion/analysis/ZWK/discovery_fresh

FP=/usr/local/Cluster-Apps/ceuadmin/FragPipe/22.0
PHIL="$FP/tools/Philosopher/philosopher-v5.1.1"
FASTA="$FP/databases/UP000005640-canonical-proteins-20420.fasta"

DBDIR="$PWD/database_canonical"
mkdir -p "$DBDIR"
cd "$DBDIR"

# Initialise the Philosopher workspace.
"$PHIL" workspace --init

# Generate the target/decoy database.
"$PHIL" database --custom "$FASTA" --prefix rev_

# Inspect the generated FASTA. Confirm its actual path before proceeding.
DB="$PWD/2026-10-09-decoys-UP000005640-canonical-proteins-20420.fasta.fas"

test -s "$DB" || {
    echo "ERROR: Expected database not found: $DB" >&2
    echo "Inspect the generated files and verify the output filename." >&2
    exit 1
}

printf 'Total FASTA records: '
grep -c '^>' "$DB"

printf 'Decoy FASTA records: '
grep -c '^>rev_' "$DB"

printf 'First five FASTA headers:\n'
grep '^>' "$DB" | head -5

# Record the checksum for provenance.
sha256sum "$FASTA" "$DB"
```

For the existing production database, the expected counts are 40,840 total sequences, including 20,420 decoys and 20,420 targets.

If Philosopher produces a different filename or the counts differ, stop and investigate before using the database. Record the actual filename and checksum for any newly generated version.

## Version control

The pipeline scripts and documentation should be version-controlled. Generated search outputs, temporary test environments and reports should be retained according to the project's data-management policy rather than committed indiscriminately.
