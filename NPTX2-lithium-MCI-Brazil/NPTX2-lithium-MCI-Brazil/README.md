# NPTX2-lithium-MCI-Brazil

**Code repository for:**
Xiao M-F\*, Kim K\*, et al. *Lithium upregulates NPTX2 in association with cognitive maintenance.* (2026)

\*These authors contributed equally.
Correspondence: Paul F. Worley (pworley@jhmi.edu)

---

## Overview

This repository contains the analysis code supporting the proteomics and biomarker analyses in the manuscript. Data were derived from archived plasma and CSF samples from a randomized, placebo-controlled trial of low-dose lithium carbonate in amnestic mild cognitive impairment (MCI) (Sao Paulo, Brazil trial; Forlenza et al.).

CSF proteins were quantified by three platforms: **ELISA**, **NULISA (CNS Disease Panel 120)**, and **PRM-MS**.

---

## Repository Structure

```
notebooks/
├── 00_data_preprocessing.ipynb
│       Split visit-difference data by treatment arm (Lithium / Placebo)
│       and assay type (ELISA, NULISA, PRM)
│
├── 01_litemia_protein_correlation.ipynb
│       Pearson correlation: plasma lithium (Litemia) vs. CSF protein
│       changes across all analytes; bubble plot visualization
│
├── 02_baseline_proteomics_correlation.ipynb
│       Pearson correlations at baseline between demographic/cognitive
│       variables and CSF proteins (ELISA, NULISA, PRM)
│
├── 03_NULISA_correlation_bubble_plot.ipynb
│       Bubble plot figure: per-protein correlation coefficients across
│       four groups (placebo baseline, lithium baseline, placebo 1y, lithium 1y)
│
└── 04_YWHAZ_NPTX2_cognition_correlation.ipynb
        Pearson correlations between YWHAZ:NPTX2 ratio and cognitive measures

scripts/
└── revision_analysis.py
        Complete, self-contained reanalysis performed for the 2026 revision:
        quality control, metadata assembly, multiple-testing correction,
        between-arm comparisons, interaction tests, network analyses,
        regression-to-the-mean diagnostics, Gene Ontology enrichment,
        quality-control sensitivity analyses, and all revised figures.
```

---

## Analysis performed for the 2026 revision

`scripts/revision_analysis.py` reproduces the revised analyses end to end from
the source data files. It is organized in numbered steps and writes every table
as `.xlsx` and every figure as vector `.pdf`.

Key methodological specifications:

- **NULISA analyte inclusion**: detectability >= 50% of samples above the
  plate-specific limit of detection, following the platform vendor's
  recommendation for differential-expression analysis. Values below the LOD but
  non-zero are retained.
- **Sample inclusion**: >= 70% of targets above LOD (the vendor threshold for CSF).
- **Quantification scale**: NULISA NPQ (log2) and PRM-MS (log2) as provided;
  ELISA log2-transformed so that all three platforms are analyzed on a
  comparable scale.
- **PRM-MS quantification**: absolute concentrations from heavy/light stable
  isotope-labeled peptide ratios, normalized to the peptide-wise median across
  samples, log2-transformed, then averaged across peptides to protein level.
- **Multiple testing**: Benjamini-Hochberg FDR applied within each quantification
  platform (primary) and across the pooled analyte set (sensitivity).
- **Between-arm comparison**: ANCOVA on the one-year change with the baseline
  value as covariate; models with and without age, sex and APOE e4 are reported.
- **Group difference in correlation**: Fisher r-to-z and a formal
  analyte x treatment-arm interaction model.
- **Network analyses**: label-permutation and sample-size-matched down-sampling,
  with a threshold-independent mean absolute Fisher-z measure reported alongside
  edge counts.
- **Regression to the mean**: because baseline CSF NPTX2 differs between arms,
  four diagnostics are reported: the within-arm correlation between the baseline
  value and the subsequent change, the cross-platform reliability of the baseline
  measurement, an ANCOVA adjusted for a composite baseline averaged across the
  three platforms, and a plate-by-arm contingency test.
- **Quality-control sensitivity**: the one NULISA sample falling just below the
  vendor's 70% sample-detectability threshold was retained, because excluding it
  would have broken a longitudinal pair; all NULISA NPTX2 analyses are repeated
  with that participant removed. The plate composition of the longitudinal pairs
  is reported alongside, since a plate-by-arm test alone does not address plate
  effects on within-person change.
- **Enrichment analysis**: Fisher exact tests against Gene Ontology libraries
  using the measured analyte panel as the background set (primary), because the
  panel is composed largely of synaptic and neuronal proteins and a whole-genome
  background would return enrichment reflecting the panel design. The
  whole-genome background is reported alongside it for comparison.

### Output files

Running the script writes the following to `outputs/tables/` and
`outputs/figures/`:

| File | Contents |
| --- | --- |
| `meta_master.xlsx` | Long-format merged metadata, one row per CSF sample |
| `analysis_ready.xlsx` | Participant-level analysis matrix (baseline values and one-year change scores) |
| `Table1_baseline_characteristics.xlsx` | Baseline characteristics by treatment arm |
| `S_Table_QC_and_sample_disposition.xlsx` | Analyte and sample quality control, analysis set sizes, baseline balance |
| `S_Table_association_results.xlsx` | All correlations with FDR, between-arm comparison, interaction tests, within-arm change, leave-one-out diagnostics |
| `S_Table_network_analysis.xlsx` | Correlation-change and network comparisons, and the reconciliation of the originally reported result |
| `S_Table_figure2_values.xlsx` | Plotted values for Figure 2 |
| `S_Table_RTM_diagnostics.xlsx` | Regression-to-the-mean diagnostics |
| `S_Table_GO_enrichment.xlsx` | Gene Ontology enrichment under both background sets |
| `S_Table_lowQC_sensitivity.xlsx` | NULISA low-QC sample sensitivity analysis and plate composition of the longitudinal pairs |

---

## Dependencies

Developed in **Google Colab** (Python 3).

```
pandas
numpy
scipy
statsmodels
matplotlib
openpyxl
gseapy
```

Install via:
```bash
pip install pandas numpy scipy statsmodels matplotlib openpyxl gseapy
```

---

## Data Availability

Raw mass spectrometry data are deposited in **ProteomeXchange / PRIDE**:

- PRM-MS data: `PXD077939`

Peptides were analyzed on an **Orbitrap Fusion Lumos Tribrid** mass spectrometer
coupled to an **Ultimate 3000 RSLCnano** liquid chromatography system
(Thermo Fisher Scientific). Data were processed in Skyline.

Source data files (Excel) referenced in the notebooks are available upon
reasonable request from the corresponding author, pending data sharing
agreements with the Sao Paulo trial team.

> **Note:** Notebook file paths reference Google Drive (`/content/drive/MyDrive/...`).
> Users running locally should update file paths to match their environment.

---

## Citation

> Xiao M-F\*, Kim K\*, Lissit K, Huang Y, Lao Y, Fan M, Talib LL, Huentelman MJ,
> Barnes CA, Forlenza OV, Xu J-C, Na CH, Worley PF.
> *Lithium upregulates NPTX2 in association with cognitive maintenance.* (2026)

---

## License

Released under the [MIT License](LICENSE).

---

## Contact

Kyungdo Kim, Ph.D.
Department of Neurology, Institute for Cell Engineering
Johns Hopkins University School of Medicine
