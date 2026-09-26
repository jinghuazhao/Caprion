# Raw, Normalized Peptides and Protein Abundance

## Overview

The Caprion pilot workbook contains several levels of processed proteomics data. For protein-level abundance analysis, the recommended endpoint is:

**`Protein_All_Peptides`**

The underlying raw and peptide-level data are useful for quality control, troubleshooting, and methodological investigation, but the supplied workbook does not contain enough information to reliably reconstruct the exact Caprion protein-abundance algorithm.

---

## Data hierarchy

The practical data flow appears to be approximately:

```text
Raw IGs
   |
   v
Peptide / feature processing
   |
   v
Normalized Peptides
   |
   v
Protein-level processing / peptide selection or weighting
   |
   v
Protein_All_Peptides
   |
   +--> Protein_DR_Filt_Peptides
```

This represents the observed structure of the supplied data. It should not be interpreted as a documented Caprion algorithm, because the workbook contains exported results rather than the underlying calculation workflow.

---

## 1. Raw IGs

### Sheet

`Raw IGs`

### Contents

This is the lowest-level data available in the workbook. Each row represents an **Isotope Group (IG)** associated with a peptide/protein.

Important fields include:

- `Isotope Group ID`
- `Protein`
- `Modified Peptide Sequence`
- `Monoisotopic m/z`
- `Max Isotope Time Centroid`
- `Charge`
- `ZWK0001` ... `ZWK0200`

The sample columns contain raw intensity-like measurements.

### Example: 1433B_HUMAN

Six isotope groups are present:

- `442662695` — `AKLAEQAERYDDMAAAMK`
- `442590982` — `AVTEQGHELSNEER`
- `442706125` — `AVTEQGHELSNEER`
- `442627470` — `AVTEQGHELSNEERNLLSVAYK`
- `442807041` — modified long peptide
- `442641664` — `YLSEVASGDNK`

There are therefore five unique peptide sequences, because two isotope groups correspond to `AVTEQGHELSNEER`.

### Recommended use

Use Raw IGs when investigating:

- raw MS signal behaviour;
- isotope groups;
- peptide identification;
- retention time and m/z;
- charge;
- missingness;
- raw-data quality control;
- attempts to reproduce the upstream processing pipeline.

For ordinary protein-abundance analysis, Raw IGs are **not the preferred endpoint**.

---

## 2. Normalized Peptides

### Sheet

`Normalized Peptides`

This sheet has essentially the same peptide/IG structure as `Raw IGs`, but the sample measurements have been transformed into normalized peptide-level abundance values.

For example, for `1433B_HUMAN`, the raw value and normalized value for the same IG can be very different in scale:

```text
Raw intensity       Normalized peptide
125873.3281         17.098838
216674.1250         17.779611
198077.9062         17.785887
...
```

### Relationship to raw data

The normalized values are strongly related to `log2(raw intensity)` for several peptides, but they are **not simply**:

```text
Normalized = log2(Raw)
```

nor a single universal linear transformation.

For 1433B_HUMAN, examples of regressions of normalized abundance on `log2(raw)` were:

```text
442590982    slope 0.9547   correlation 0.9865
442706125    slope 0.9710   correlation 0.9939
442627470    slope 0.8517   correlation 0.9071
442662695    slope 0.5878   correlation 0.6968
442641664    slope 1.0069   correlation 0.6992
442807041    slope -0.0297  correlation -0.0172
```

This demonstrates that the normalization cannot be reconstructed as one simple transformation of raw intensity.

### Important finding

Within-IG standardization of `log2(raw)` also did not explain the normalized values. The six IG-specific slopes ranged from approximately `-0.03` to `1.40`.

### Recommended use

`Normalized Peptides` is useful for:

- peptide-level analysis;
- peptide QC;
- investigating peptide behaviour;
- understanding the inputs to protein-level processing.

It is **not necessary to reconstruct these values** if the objective is simply to obtain protein abundances from the supplied processed dataset.

---

## 3. Protein_All_Peptides

### Sheet

`Protein_All_Peptides`

This is the most appropriate dataset for the main **protein-level abundance analysis**.

The structure is:

```text
Protein | ZWK0001 | ZWK0002 | ... | ZWK0200
```

Each row represents a protein and each sample column contains the reported protein abundance.

For example:

```text
1433B_HUMAN
ZWK0001 = 17.3366328529562
ZWK0002 = 17.5882866329177
ZWK0003 = 18.6032096647269
...
```

### What we learned about its construction

We tested a number of simple hypotheses using the underlying normalized peptide values.

The reported protein abundance is **not adequately reproduced by**:

- the arithmetic mean of all normalized peptide/IG values;
- the median of all normalized peptide/IG values;
- the maximum peptide;
- the minimum peptide;
- the geometric mean;
- equal-weight aggregation of unique peptide sequences;
- simply selecting the strongest raw peptide;
- detection rate;
- lower coefficient of variation alone.

Across proteins with exactly two peptide sequences, one peptide was often much more strongly associated with the reported protein abundance than the other. This suggests that peptide evidence is not treated as a simple equal-weight average, although the exact mechanism cannot be established from the exported tables.

For example, among 128 proteins with exactly two unique peptide sequences:

- median correlation of the better peptide with protein abundance: approximately `0.946`;
- median correlation of the other peptide: approximately `0.345`;
- the better peptide had lower regression RMSE in all 128 cases.

This is strong evidence that peptide-level contributions are not simply equal-weighted.

### What cannot be established

The workbook does **not** contain the calculation formulas or a hidden calculation sheet that explains how `Protein_All_Peptides` is generated.

All workbook sheets were visible:

```text
Legend
Samples
Annotations
Raw IGs
Normalized Peptides
Protein_All_Peptides
Protein_DR_Filt_Peptides
```

The numerical values inspected in the calculation sheets are stored values rather than Excel formulas.

Therefore, the exact Caprion algorithm used to produce `Protein_All_Peptides` cannot be reliably reconstructed from this workbook alone.

### Recommendation

**Use `Protein_All_Peptides` as the primary protein abundance dataset.**

There is no need to reverse-engineer the upstream calculation merely to obtain protein abundance values.

---

## 4. Protein_DR_Filt_Peptides

### Sheet

`Protein_DR_Filt_Peptides`

This is another protein-level abundance table, apparently produced after a detection-rate/filtering step.

For `1433B_HUMAN`, the first values illustrate the small but measurable difference:

```text
Protein_All_Peptides       Protein_DR_Filt_Peptides

17.336633                  17.365318
17.588287                  17.630887
18.603210                  18.652622
17.086007                  17.179862
17.689614                  17.773327
...
```

Across the 200 samples for 1433B_HUMAN, the difference between DR-filtered and All-Peptides values was approximately:

```text
Mean difference       +0.0393
SD                     0.0309
Median difference      +0.0392
Minimum                -0.0460
Maximum                +0.1225
```

### Recommended use

Retain this dataset if:

- detection-rate filtering is scientifically relevant;
- you want a sensitivity analysis;
- you need to compare filtered and unfiltered protein abundance.

For the main analysis, use `Protein_All_Peptides` unless there is a specific reason to adopt the DR-filtered version.

---

## 5. Raw ZIP files

The raw/raw.zip material is potentially important for **reproducibility of the upstream proteomics workflow**, but it is considerably less useful for ordinary downstream protein-abundance analysis.

The raw material becomes important if the goal is to:

- reproduce the complete Caprion pipeline;
- investigate peptide identification;
- reprocess the mass-spectrometry data;
- apply different QC thresholds;
- investigate missingness or detection;
- understand how normalized peptides were generated;
- independently calculate protein abundance.

However, without the exact Caprion processing software, parameters, reference data, and workflow documentation, the raw files alone may not be sufficient to reproduce the exact values in `Protein_All_Peptides`.

---

## 6. Practical recommendation

For the current project, use the following hierarchy:

| Dataset | Main purpose | Recommended for protein analysis |
|---|---|---|
| `Protein_All_Peptides` | Reported protein abundance | **Yes — primary** |
| `Protein_DR_Filt_Peptides` | Detection-rate filtered abundance | Yes — sensitivity analysis |
| `Normalized Peptides` | Processed peptide abundance | Optional / QC |
| `Raw IGs` | Raw isotope-group measurements | Mainly QC / reprocessing |
| `raw.zip` | Underlying raw material | Mainly full reprocessing |

### Suggested analysis dataset

The primary matrix should therefore be:

```text
Protein_All_Peptides
```

with:

```text
rows    = proteins
columns = ZWK samples
values  = reported protein abundance
```

The DR-filtered matrix can be retained as a secondary analysis.

---

## 7. Bottom line

The investigation established that:

1. `Protein_All_Peptides` is the appropriate **protein-level endpoint** in the supplied workbook.
2. `Protein_All_Peptides` cannot be reproduced reliably from simple aggregation of `Normalized Peptides`.
3. `Normalized Peptides` is related to raw log2 intensity but is not generated by one simple universal transformation.
4. Raw IG abundance, detection rate, and peptide variability alone do not explain which peptide evidence contributes most strongly to the protein result.
5. The workbook contains no hidden calculation sheets and the relevant values are exported/stored rather than formula-driven.
6. Consequently, the exact Caprion protein-calculation algorithm is not recoverable from the workbook alone.
7. For downstream statistical analysis, it is preferable to use the supplied `Protein_All_Peptides` values rather than attempting to reconstruct them from the raw data.

**Primary recommendation:**

> Use `Protein_All_Peptides` as the main protein abundance dataset, and retain `Protein_DR_Filt_Peptides` for sensitivity or filtered analyses. Treat `Normalized Peptides`, `Raw IGs`, and `raw.zip` as upstream/QC/reprocessing data rather than replacing the supplied protein-level results.
