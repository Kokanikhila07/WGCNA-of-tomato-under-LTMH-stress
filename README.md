# WGCNA-of-tomato-under-LTMH-stress

# WGCNA Analysis of Tomato Pollen Abortion under LTMH Stress (GSE185629)

This repository contains the **Weighted Gene Co-expression Network Analysis (WGCNA)** pipeline used to characterize the molecular-physiological basis of pollen abortion in *Solanum lycopersicum*.

## 📌 Project Overview

This study evaluates the transcriptomic response of developing tomato anthers to **Long-Term Mild Heat (LTMH)** stress. While the original experiment covered multiple stages (tetrad, free microspore, and polarized microspore), this analysis specifically focuses on the **polarized microspore stage**, as it exhibits the most significant damage under prolonged heat stress.

## 📊 Dataset Details

* **Accession:** [GSE185629](https://www.google.com/search?q=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi%3Facc%3DGSE185629)
* **Organism:** *Solanum lycopersicum* (Tomato)
* **Condition:** Control vs. Long-Term Mild Heat (LTMH)
* **Target Stage:** Polarized microspore stage anthers

## 🛠️ Methodology

The analysis was performed using the R package `WGCNA` with the following steps:

### 1. Data Preprocessing

* Expression data (1,468 genes across 8 samples) was transposed to a samples-by-genes format.
* Sample names were standardized, and outlier detection was performed via hierarchical clustering.

### 2. Network Construction

* **Soft Thresholding**: A soft power of 6 was selected for network topology.
* **Blockwise Modules**: Modules were detected using an unsigned TOM type with a minimum module size of 30.

### 3. Module-Trait Relationship

* Module Eigengenes (MEs) were correlated with environmental traits (Control vs LTMH).
* Significant modules were identified using Pearson correlation and Student asymptotic p-values.

## 📈 Key Results

The WGCNA identified critical modules highly correlated with the heat stress condition:

* **Turquoise Module (ME1):** Showed a strong positive correlation of **0.985** (p = 7.77*10^{-6}).
* **Blue Module (ME2):** Showed a strong negative correlation of **-0.926** (p = 9.53*10^{-4}).

### 📁 Repository Structure

Based on your preferred layout, organize your project as follows:

```text
Tomato_Pollen_LTMH_WGCNA/
├── data/
│   ├── expr_matrix_GSE185629.txt    # Expression matrix (Polarized stage)
│   └── metadata.txt                 # Traits: Control vs LTMH
├── scripts/
│   └── WGCNA GSE185629 Script.R     # Your cleaned WGCNA script
├── results/
│   ├── Module_trait_heatmap.png     # Heatmap showing ME1/ME2 correlations
│   ├── CytoscapeInput-edges-turquoise.txt
│   └── CytoscapeInput-nodes-turquoise.txt
└── README.md                        # Documentation

```

## 🚀 Usage

To replicate this analysis:

1. Ensure `WGCNA` and `dplyr` packages are installed in R.
2. Place the expression and metadata files in the `/data` directory.
3. Run `scripts/WGCNA_GSE185629.R`.
