# NPTX2-lithium-MCI-Brazil

**Code repository for:**  
Xiao M-F\*, Kim K\*, et al. *Lithium upregulates NPTX2 in association with cognitive maintenance.* (2026)

†These authors contributed equally.  
\*Correspondence: Paul F. Worley (pworley@jhmi.edu)

---

## Overview

This repository contains the analysis code supporting the proteomics and biomarker analyses in the manuscript. Data were derived from archived plasma and CSF samples from a randomized, placebo-controlled trial of low-dose lithium carbonate in amnestic mild cognitive impairment (MCI) (São Paulo, Brazil trial; Forlenza et al.).

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
│       Pearson correlation: plasma lithium (Litemia) change vs.
│       CSF protein changes across all analytes; bubble plot visualization
│
├── 02_baseline_proteomics_correlation.ipynb
│       Pearson correlations at baseline (_base groups) between
│       demographic/cognitive variables and CSF proteins (ELISA, NULISA, PRM)
│
├── 03_NULISA_correlation_bubble_plot.ipynb
│       Bubble plot figure: per-protein correlation coefficients across
│       four groups (Base placebo, Baseline lithium, 1y placebo, 1y lithium)
│       formatted for publication (Nature single-column, 1200 dpi)
│
└── 04_YWHAZ_NPTX2_cognition_correlation.ipynb
        Pearson correlations between YWHAZ:NPTX2 ratio (ELISA/NULISA/PRM)
        and cognitive measures (CDR, CDR-SB, ADAS, SLN, TMT-A/B, APOE ε4)
        at baseline, 1-year follow-up, and longitudinal change scores
```

---

## Dependencies

All notebooks were developed in **Google Colab** (Python 3).

```
pandas
numpy
matplotlib
seaborn
scipy
openpyxl
```

Install via:
```bash
pip install pandas numpy matplotlib seaborn scipy openpyxl
```

---

## Data Availability

Raw mass spectrometry data are deposited in **ProteomeXchange / PRIDE**:  
- PRM-MS data: `PXD059266` (Shimadzu LCMS-8060 triple quadrupole)

Source data files (Excel) referenced in notebooks are available upon reasonable request from the corresponding author, pending data sharing agreements with the São Paulo trial team.

> **Note:** Notebook file paths reference Google Drive (`/content/drive/MyDrive/...`). Users running locally should update file paths to match their environment.

---

## Citation

If you use this code, please cite:

> Xiao M-F\*, Kim K\*, Lissit K, Huang Y, Lao Y, Fan M, Talib LL, Huentelman MJ, Barnes CA, Forlenza OV, Xu J-C, Na CH, Worley PF. *Lithium upregulates NPTX2 in association with cognitive maintenance.* (2026)

---

## License

This code is released under the [MIT License](LICENSE).

---

## Contact

Kyungdo Kim, Ph.D.  
Department of Neurology, Institute for Cell Engineering  
Johns Hopkins University School of Medicine
