# An Aging Clock Based on Immune Repertoire Features: COVID-19 Accelerates Aging

![Workflow](https://img.shields.io/badge/Workflow-TidyModels-blue.svg)
![R Version](https://img.shields.io/badge/R-4.0%2B-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

This repository implements an XGBoost-based aging prediction model using immune repertoire features, demonstrating accelerated immunological aging in COVID-19 patients.

## Requirements
- R >= 4.0
- Packages: `tidymodels`, `xgboost`, `iml`, `fastshap`, etc.

## Code Structure
```r
# Key Steps:
1. Data Loading & Preprocessing
2. XGBoost Model Training with Hyperparameter Tuning
3. Model Evaluation & Prediction
4. Feature Importance & SHAP Analysis
