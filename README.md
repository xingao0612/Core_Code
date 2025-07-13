# Aging Clock Based on Immune Repertoire Features

This repository contains the analysis code and pipelines used in the study titled **"An aging clock based on immune repertoire features: COVID-19 accelerates aging"**.

## 📂 Overview

The code supports data preprocessing, statistical analysis, immune repertoire profiling, and machine learning modeling (XGBoost-based age prediction) across multiple cohorts.

## 🧬 Key Components

- **Data Integration & Cleaning**
  - Merge clinical metadata with TCR/BCR repertoire metrics.
  - Segment-specific CDR3 analysis (nt & aa length).
  - Frequency normalization and transformation.

- **Immune Repertoire Analyses**
  - CDR3 length distribution and aging correlation.
  - Amino acid usage profiling and age correlation.
  - TRBV gene usage and significance screening.
  - Clonality, diversity, and repertoire overlap using `immunarch`.

- **Machine Learning**
  - XGBoost regression for immune age prediction.
  - SHAP-based interpretation.
  - COVID-19 cohort evaluation.

- **Visualization**
  - Publication-ready plots using `ggplot2`, `ComplexHeatmap`, `ggpubr`, `ggsci`, etc.

## 🧪 Requirements

- R >= 4.1
- Packages:
  - `tidymodels`, `xgboost`, `iml`, `fastshap`
  - `ggplot2`, `ggpubr`, `ggsci`, `pheatmap`, `ComplexHeatmap`
  - `immunarch`, `Hmisc`, `readxl`, `dplyr`, `patchwork`

Install missing packages via:

```r
install.packages("your_missing_package")
```

## 📁 File Structure

- `*.csv` – Expression matrices and metadata
- `*.Rdata` – Intermediate analysis objects
- `*.pdf` – Plots for publication
- `*.R` or `.txt` – Core analysis scripts

## ⚙️ Running the Pipeline

1. Prepare input files (`.csv`, `.xlsx`) with correct sample IDs.
2. Source the script or run sections interactively.
3. Outputs will be saved as `.pdf` (plots) and `.csv` (results).

## 🧾 Citation

If you use this code, please cite our paper:

> *An aging clock based on immune repertoire features: COVID-19 accelerates aging*. [Journal Info]

---

**Contact**: [Author Name] – [Institution] – [Email or GitHub]
