# ZWK Independent Reference Analysis

(Experimental with ZWK)

## 1. Objective

Develop an independent reference from the original ZWK pooled mass-spectrometry acquisitions using the original mzML data.

The analysis has four stages:

1. Identify peptides from MS/MS spectra using FragPipe.
2. Audit peptide-ion evidence across acquisitions.
3. Detect MS1 features and generate isotope groups (IGs) independently.
4. Link peptide identifications to independent IGs and compare the resulting reference with the original XLSX assignments.

**Independence requirement:** The original XLSX isotope-group assignments must not be used to construct the independent reference. They are reserved for subsequent comparison.

## 2. Project layout and environment

Project root:

`~/Caprion/analysis/ZWK/reference_ZWK/`

| Directory | Purpose |
|---|---|
| `config/` | FragPipe workflow and input manifest |
| `database/` | Target and decoy protein sequence database |
| `dry_run_98/` | FragPipe dry-run workspace |
| `search_98/` | Completed FragPipe search outputs |
| `reports/` | Audit tables, validation reports and exported results |
| `scripts/` | Numbered workflow and audit scripts |

FragPipe version: **22.0**

Load FragPipe when required:

```bash
module load ceuadmin/FragPipe/22.0
```

The original mzML acquisitions are located at:

`/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/mzML/`

The search manifest contains 98 mzML input paths.

## 3. Database and search configuration

Database:

`database/2026-10-12-decoys-UP000005640-canonical-proteins-20420.fasta.fas`

The FASTA contains 20,420 target entries and 20,420 reverse-decoy entries.

The candidate workflow is:

`config/ZWK-reference-panel-candidate.workflow`

The input manifest is:

`config/ZWK-reference-panel.fp-manifest`

The workflow configuration is the authoritative record of search parameters. The search used MSFragger, Philosopher and IonQuant, with label-free quantification and match-between-runs (MBR) enabled.

The search completed successfully. Its log is:

`reports/fragpipe_search_98.log`

The completed search should not be rerun merely to reproduce the existing outputs.

## 4. Numbered script inventory

All paths below are relative to the project root:

`~/Caprion/analysis/ZWK/reference_ZWK/`

| Script | Responsibility | Status |
|---|---|---|
| `00_check_inputs.py` | Validate required inputs and the search manifest | Completed |
| `01_audit_search_config.py` | Inspect the FragPipe workflow configuration | Completed |
| `02_run_fragpipe_search.sh` | Run FragPipe in headless mode | Completed |
| `03_audit_ion_evidence.py` | Audit the combined peptide-ion evidence | Completed |
| `04_build_ion_evidence.py` | Export ion-by-run evidence in long format and summarise it | Completed |
| `05_build_candidate_panel.py` | Select the provisional candidate ion panel using recurrence of direct MS/MS evidence | Completed |
| `06_validate_spectrum_scan_mapping.py` | Validate FragPipe PSM spectrum names against native mzML MS2 scans | Completed for `xms1`–`xms5` of `szwk021704i19101` |

### Execution order

The numbered scripts describe the current workflow in order. The full FragPipe search has already completed, so avoid rerunning it unnecessarily.

```bash
cd ~/Caprion/analysis/ZWK/reference_ZWK

# Input and configuration checks
python scripts/00_check_inputs.py
python scripts/01_audit_search_config.py

# Run only when a new search is intended
bash scripts/02_run_fragpipe_search.sh

# Evidence audits and exports
python scripts/03_audit_ion_evidence.py
python scripts/04_build_ion_evidence.py
python scripts/05_build_candidate_panel.py

# Scan-mapping validation for the supported acquisitions
python scripts/06_validate_spectrum_scan_mapping.py xms1
python scripts/06_validate_spectrum_scan_mapping.py xms2
python scripts/06_validate_spectrum_scan_mapping.py xms3
python scripts/06_validate_spectrum_scan_mapping.py xms4
python scripts/06_validate_spectrum_scan_mapping.py xms5
```

Check each script's current command-line interface before running it. In particular, script `06` validates spectrum-to-scan mapping; it does not generate isotope groups or establish peptide-to-IG assignments.

### Setup commands and provenance

Commands used to create the candidate workflow, target-decoy database and manifest, and to validate the FragPipe dry run should be retained in a separate setup record, such as `docs/SETUP.md`. They are useful provenance, but are not additional numbered analysis scripts.

### Planned next stage

The next script should be introduced only after the independent MS1 feature-detection and isotope-grouping method has been selected and validated. It should generate new IG identifiers from mzML evidence without using the XLSX assignments as input.

A subsequent stage should link independently generated IGs to peptide identifications using precursor m/z, charge, retention time and scan evidence. These stages are **planned, not yet implemented**.

### Running individual scripts

Examples:

```bash
cd ~/Caprion/analysis/ZWK/reference_ZWK

python scripts/00_check_inputs.py
python scripts/01_audit_search_config.py
```

The FragPipe search script is:

```bash
bash scripts/02_run_fragpipe_search.sh
```

Scripts `03`–`06` should be run according to their command-line interfaces and the intended audit or validation task.

## 5. Completed work and results

### 5.1 FragPipe search

The 98-input search completed successfully in approximately 29.5 minutes. IonQuant exited successfully.

The principal outputs include:

- `search_98/combined_ion.tsv`
- `search_98/combined_modified_peptide.tsv`
- `search_98/combined_peptide.tsv`
- `search_98/combined_protein.tsv`

### 5.2 Ion-evidence audit

The following reports were generated:

- `reports/ion_evidence_long.tsv`
- `reports/ion_evidence_summary.tsv`
- `reports/candidate_reference_ion_panel.tsv`
- `reports/filter_fdr_audit.tsv`

The combined ion table contains 27,984 peptide–charge ion entries.

The long-form evidence table contains 2,630,496 rows, representing ion evidence across acquisitions. The summary table contains 27,984 unique ion entries.

Across run-level evidence cells, the recorded match types were:

| Match type | Count |
|---|---:|
| Direct MS/MS | 169,413 |
| Match between runs (MBR) | 6,360 |
| Unmatched | 2,454,723 |

These counts describe the evidence represented in the combined ion table; they are not counts of unique peptides or biological analytes.

### 5.3 Provisional candidate ion panel

The provisional panel was selected using direct MS/MS evidence in at least two acquisitions.

It contains:

- 23,112 ion entries
- 2,117 protein identifiers
- 2,116 non-missing gene labels

This recurrence threshold is a reproducibility criterion. It is not a separate combined-panel false discovery rate (FDR) estimate.

Direct MS/MS identification must be distinguished from MBR evidence, which transfers quantitative evidence between runs and does not constitute an independent MS/MS identification in the receiving run.

### 5.4 Per-run FDR audit

Per-run Philosopher target-decoy filtering was audited, with final counts available for 94 runs.

The absence of reverse-decoy entries in a combined final output does not, by itself, establish a separately controlled FDR for the combined candidate panel.

### 5.5 PSM-to-native-MS2 scan mapping

For run `szwk021704i19101`, FragPipe PSM spectrum names were matched to native MS2 scan identifiers in five mzML acquisitions.

Report:

`reports/spectrum_scan_mapping_5_acquisitions.tsv`

| Acquisition | PSMs | Scan IDs found | Charge agreement |
|---|---:|---:|---:|
| xms1 | 6,866 | 6,866 | 100% |
| xms2 | 6,974 | 6,974 | 100% |
| xms3 | 6,841 | 6,841 | 100% |
| xms4 | 6,976 | 6,976 | 100% |
| xms5 | 7,108 | 7,108 | 100% |
| **Total** | **34,765** | **34,765** | **100%** |

The maximum absolute retention-time difference reported was approximately 0.00034 seconds.

This validates spectrum-to-scan mapping for these five acquisitions. It does **not** establish an association between each PSM and a particular MS1 isotope group.

## 6. Setup and command provenance

The following commands document important setup operations. They are command history, not additional numbered analysis scripts.

### 6.1 Copying a candidate workflow

The candidate workflow was copied from the discovery workspace, with the source preserved and existing destinations protected from accidental overwrite.

Source:

`~/Caprion/analysis/ZWK/discovery_fresh/LFQ-MBR-Sage-test.workflow`

Destination:

`config/ZWK-reference-panel-candidate.workflow`

### 6.2 Generating the target-decoy database

The database was generated from the canonical UniProt reference FASTA using the Philosopher version supplied with FragPipe 22.0.

The original target FASTA was:

`/usr/local/Cluster-Apps/ceuadmin/FragPipe/22.0/databases/UP000005640-canonical-proteins-20420.fasta`

The resulting target-decoy FASTA is retained under `database/`.

### 6.3 Building the manifest

The manifest was generated from the 98 mzML files. Acquisition labels were derived from file names as neutral reproducible labels; they should not be interpreted as biological sample-group assignments without independent evidence.

### 6.4 Dry-run validation

The FragPipe dry run used:

`config/ZWK-reference-panel-candidate.workflow`

and:

`config/ZWK-reference-panel.fp-manifest`

Its log is:

`reports/fragpipe_dry_run.log`

The dry run and input audit are separate from the completed full search.

## 7. Work remaining

### 7.1 Independent MS1 feature detection and isotope grouping

**This stage has not yet been completed.**

The next step is to inspect the existing IonQuant outputs for suitable MS1 feature-level evidence. If these outputs do not provide sufficient information for independent IG reconstruction, select and document an independent feature-detection and isotope-grouping method.

OpenMS is a candidate open-source framework for MS1 feature detection and isotope-related processing.

Rosetta Elucidator 3.0 / PeakTeller has been mentioned as a possible source of the original IGs, but this provenance remains unconfirmed.

New isotope groups must receive independent identifiers. They must not be assumed to correspond exactly to the original XLSX IG IDs.

### 7.2 Linking peptide identifications to independent IGs

After generating the independent IGs, investigate peptide-to-IG associations using:

- Precursor m/z
- Charge state
- Retention time
- MS/MS scan and precursor evidence
- The relevant MS1 feature and isotope-pattern evidence

Associations should be supported by the available data and should distinguish direct MS/MS identifications from transferred MBR evidence.

### 7.3 Exploratory anchor-stability analysis

Exploratory scripts `07`–`21` have investigated isotope-anchor scoring, mass accuracy, scan-level stability and competing isotope offsets using the existing OpenMS pilot outputs.

The audits show that low individual peak mass errors and apparently complete isotope series do not necessarily identify a unique monoisotopic anchor. Competing offsets can receive similar scores, and ties remain substantial in some candidate groups.

Scripts `20` and `21` and their outputs are retained for diagnostic purposes only. Their results must not be treated as validated monoisotopic assignments:

- Script `20` reports winning-offset rows, which can exceed the number of unique candidates when ties occur.
- Script `21` reuses observed isotope m/z values from an existing audit rather than independently rematching trace peaks for each alternative anchor.

**Current status:** the anchor-scoring method remains exploratory. No final independent isotope-group assignments have been established by this analysis.

Before further scoring or scaling to additional acquisitions, inspect the trace-peak schema and existing matching implementation, verify isotope-index and offset conventions, and validate the rematching procedure on a single candidate using distinct observed trace peaks. Only after this validation should the method be applied across the candidate set.

### 7.4 Comparison with the original XLSX

Only after constructing the independent reference should it be compared with the original workbook:

`~/Caprion/pre_qc_data/pilot/ZWK_EDR_20191002.xlsx`

The workbook is a comparison source, not a construction source for the independent reference.
## 8. Interpretation notes

- The FragPipe search is complete; independent IG generation remains pending.
- A peptide–charge ion entry is not necessarily a unique biological analyte.
- Per-run FDR filtering does not automatically establish a separate FDR for the combined candidate panel.
- MBR evidence is not equivalent to a direct MS/MS identification in the receiving run.
- Successful PSM-to-scan mapping does not establish a PSM-to-IG relationship.
- IonQuant feature definitions are not necessarily equivalent to the original
- Exploratory anchor scores and stability summaries are diagnostic evidence, not validated monoisotopic assignments.
- Competing isotope offsets must be evaluated by independently matching observed trace peaks to theoretical isotope positions; existing summary-level peak matches must not be relabelled as evidence for alternative anchors.
- The independent isotope-group reference remains incomplete until the feature-detection, isotope-grouping and peptide-to-IG association procedures have been validated.
