# Sex Differences in Infant EEG Development

Code repository for:

Blanco-Gomez, G. et al. (2025). *Sex differences in early infant brain developmen*

## Overview

This repository contains the analysis code and plotting scripts used to examine biological sex as a primary variable in the early development of EEG-based brain metrics and language outcomes in infants at elevated and typical likelihood of autism spectrum disorder (ASD).


## Repository Structure

```
infant_sex_diff/
    datasets/               Input data files (not shared publicly)
    code/
        sex analyses R/     EEG sex difference models (run this first)
        clustering analyses R/  Language growth and cluster models (run second)
        figures/            Publication Figures notebook (run last)
    figures/
        main/               Main publication figures (output)
        supplementary/      Supplementary figures (output)
    tables/                 Output tables
```

## How to Run

The scripts should be run in the following order:

### Step 1: EEG Sex Differences Analysis

**File:** `code/sex analyses R/Sex differences in EEG development.R`

**What it does:**
This script models whether male and female infants differ in four EEG metrics across the first year of life (6 to 12 months), and whether sex interacts with ASD likelihood or diagnosis. The four EEG metrics are:

- Frontal gamma power
- Auditory network connectivity
- Speech network connectivity
- Power lateralization (gamma)

**Three analysis questions are addressed:**

- Q1: Sex differences in EEG trajectories across the full sample (6 to 12 months)
- Q2: Sex by ASD likelihood interaction (TLA vs ELA, full sample)
- Q3 (Exploratory): Sex by ASD diagnosis interaction (ELA subsample only)


### Step 2: Language Growth and Neurosubtype Analysis

**File:** `code/clustering analyses R/Language_sex_clustering_lme.R`

**What it does:**
This script models whether sex differences in language development (expressive and receptive) from 6 to 36 months interact with HC neurosubtype class. Neurosubtype classes (HC Class 0, 1, 2) were derived using hierarchical clustering applied to five EEG features at 6 months (see Blanco-Gomez et al. 2025b for derivation details). 

**Four analysis questions are addressed:**

- Q1: Expressive language growth by sex and HC class (primary)
- Q2: Expressive language growth by sex, HC class, and ASD likelihood (exploratory, four-way model)
- Q3: Receptive language growth by sex and HC class (primary)
- Q4: Receptive language growth by sex, HC class, and ASD likelihood (exploratory, four-way model)

**Key outputs:**
- Demographics: sex composition by HC class, chi-square tests
- Model summaries and Holm-corrected p-values for primary tests
- Post-hoc pairwise sex contrasts at 6 and 36 months within each cluster
- Summary table saved as: `tables/main/Language_HCclass_sex_corrected_results.docx`
- Supplementary figures saved to `figures/supplementary/`:
  - `S_expressive_language_sex_likelihood.png`: Expressive language trajectories by sex and ASD likelihood
  - `S_expressive_language_sex_diagnosis_ELA.png`: Expressive language by sex and ASD diagnosis (ELA only)
  - `S_expressive_language_sex_likelihood_noASD.png`: Expressive language by sex and likelihood, no-ASD participants only
 

### Step 3: Publication Figures

**File:** `code/figures/Publication Figures.ipynb`

**What it does:**
This Python (Jupyter) notebook generates all main and supplementary publication-quality figures from the pre-processed dataset. 
 

## Dependencies

### R

Tested on R 4.3+.

```r
install.packages(c(
  "tidyverse", "lme4", "lmerTest", "performance",
  "emmeans", "effectsize", "flextable", "officer",
  "lmtest", "car", "sandwich"
))
```

### Python

Tested on Python 3.10+.

```bash
pip install pandas numpy matplotlib seaborn scipy statsmodels pingouin plotly scikit-learn umap-learn
```
 

## Author

Gabriel Blanco-Gomez
