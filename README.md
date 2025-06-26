# ProteoNexus

[[Website](https://img.shields.io/badge/Website-ProteoNexus.com-blue.svg)](https://www.proteonexus.com/)
[[Data](https://img.shields.io/badge/Data-Science%20Data%20Bank-green.svg)](https://doi.org/10.57760/sciencedb.26857)

**ProteoNexus characterizes genetic architecture, estimates mediation effects, and constructs and evaluates prediction models of plasma proteome.**

This repository contains the official analysis pipeline for ProteoNexus. We provide the source code for pQTL analysis, mediation effect estimation, prediction model training, prediction execution, and data integrity validation tools.

> **IMPORTANT NOTE**
> This repository contains the **backend analysis code and scripts** only. It does not include the source code for the user-facing web server, which is available at [https://www.proteonexus.com/](https://www.proteonexus.com/).

## 🌐 Web Server

For interactive exploration, visualization, and analysis, please visit the ProteoNexus web platform:

**[https://www.proteonexus.com/](https://www.proteonexus.com/)**

The platform is freely available and does not require registration.

## ✨ Key Features

ProteoNexus is a comprehensive platform designed to decode the complex relationships between exposures, the plasma proteome, and disease incidence using data from the UK Biobank (application ID: 144904).

### 🧬 Genetic Architecture Analysis
-   Estimates pQTL summary statistics using **GEMMA**, adjusting for covariates (age, sex, top 18 genetic PCs).
-   Performs fine-mapping analysis with **SuSiE** to identify probable causal SNPs (PIP).
-   Provides interactive visualizations, including Manhattan plots and QQ-plots, on the web server.

### 🔗 Mediation Effect Estimation
-   Systematically investigates the mediating role of 2,919 plasma proteins in three distinct causal pathways:
    1.  **M-P-D**: Measurement → Protein → Disease
    2.  **E-P-D**: Environment → Protein → Disease
    3.  **G-P-D**: SNP → Protein → Disease
-   Identifies thousands of significant mediation pathways, ensuring that protein expression was measured before disease incidence.

### 📈 High-Performance Prediction Models
-   Constructs and evaluates optimized prediction models for 57 incident diseases.
-   Optimizes four machine learning models using a Tree-structured Parzen Estimator (TPE) to maximize AUC:
    -   Penalized Logistic Regression (PLR)
    -   XGBoost
    -   LightGBM
    -   Multi-Layer Perceptron (MLP)
-   The web server provides an interface to predict disease probability using our pre-trained models on user-provided proteomic data.

### 🚻 Sex-Specific Analysis
-   Implements a complete sex-specific pipeline for genetic architecture, mediation analysis, and prediction modeling to uncover sex-differentiated biological mechanisms.

## 💻 Repository Content

This repository includes the core computational scripts for:
-   **pQTL Analysis**: Scripts for running `GEMMA` and `SuSiE`.
-   **Mediation Analysis**: R scripts utilizing the `medflex` package.
-   **Model Training**: Python scripts for hyperparameter optimization of PLR, XGBoost, LightGBM, and MLP models.
-   **Prediction**: Scripts to apply the trained models to new data.
-   **Data Integrity**: Tools for quality control and data validation.

## 📊 Data Availability

The genome-wide pQTL summary statistics generated from 33,325 European participants for 2,919 plasma proteins are publicly available at the Science Data Bank.

-   **Citation**: Kaixin Shao, Peng Huang, Sheng Yang (2025) . Summary statistics for pQTL in UKB-PPP. V1. Science Data Bank.
-   **DOI**: [https://doi.org/10.57760/sciencedb.26857](https://doi.org/10.57760/sciencedb.26857)

## 📜 Citation

If you use the ProteoNexus platform or the code from this repository in your research, please cite:

> Shao K, Luo Z, Huang P, Yang S. ProteoNexus characterizes genetic architecture, estimates mediation effects, and constructs and evaluates prediction models of plasma proteome.

## 🔗 Related Projects from Our Group

-   **[GWAShug](https://github.com/caocenter/GWAShug)**: A comprehensive platform for decoding the shared genetic basis between complex traits.
    > Cao C, Tian M, Li Z, Zhu W, Huang P, Yang S. GWAShug: a comprehensive platform for decoding the shared genetic basis between complex traits based on summary statistics. *Nucleic Acids Res*. 2025 Jan 6;53(D1):D1006-D1015. doi: 10.1093/nar/gkae873.

-   **[PGS-Depot](http://www.pgsdepot.net/)**: A comprehensive resource for polygenic scores constructed by summary statistics-based methods.
    > Cao C, Zhang S, Wang J, Tian M, Ji X, Huang D, Yang S, Gu N. PGS-Depot: a comprehensive resource for polygenic scores constructed by summary statistics based methods. *Nucleic Acids Res*. 2024 Jan 5;52(D1):D963-D971. doi: 10.1093/nar/gkad1029.

-   **[PGSFusion](http://www.pgsfusion.net/)**: A server to streamline polygenic score construction and epidemiological applications.
    > Yang S, Ye X, Ji X, Li Z, Tian M, Huang P, Cao C. PGSFusion streamlines polygenic score construction and epidemiological applications in biobank-scale cohorts. *bioRxiv*.