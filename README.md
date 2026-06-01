# Sex Differences in Infant EEG Development

Code repository for:

Blanco-Gomez, G. et al. (2025). *Sex differences in early infant brain development*

## Overview

This repository contains the analysis code and plotting scripts used to examine biological sex as a primary variable in the early development of EEG-based brain metrics and language outcomes in infants at elevated and typical likelihood of autism spectrum disorder (ASD).


## Repository Structure

```
infant_sex_diff/
    datasets/               Input data files (not shared publicly)
    code/
        eeg/                EEG extraction tools
        sex analyses R/     EEG sex difference models 
        clustering analyses R/  Language growth and cluster models 
        figures/            Publication Figures notebook 
    figures/
        main/               Main publication figures (output)
        supplementary/      Supplementary figures (output)
    tables/                 Output tables
```

## How to Run

The scripts should be run in the following order:

### Step 1: EEG Variable extraction

**File:** `code/eeg/EEG_feature_extraction.ipynb

**What it does:**
This script extracts various EEG metrics across the first year of life (6 to 12 months) including:

- Frontal gamma power
- Auditory network connectivity
- Speech network connectivity
- Power lateralization (gamma)

### Step 2: EEG Sex Differences Analysis

**File:** `code/sex analyses R/Sex differences in EEG development.R`

**What it does:**
This script models whether male and female infants differ in four EEG metrics across the first year of life (6 to 12 months), and whether sex interacts with ASD likelihood or diagnosis

**Three analysis questions are addressed:**

- Q1: Sex differences in EEG trajectories across the full sample (6 to 12 months)
- Q2: Sex by ASD likelihood interaction (TLA vs ELA, full sample)
- Q3 (Exploratory): Sex by ASD diagnosis interaction (ELA subsample only)


### Step 3: Language Growth and Neurosubtype Analysis

**File:** `code/clustering analyses R/Language_sex_clustering_lme.R`

**What it does:**
This script models whether sex differences in language development (expressive and receptive) from 6 to 36 months interact with HC neurosubtype class. Neurosubtype classes (HC Class 0, 1, 2) were derived using hierarchical clustering applied to five EEG features at 6 months (see  Blanco-Gomez, G., Wright, N., O'Reilly, C. et al. (2025). EEG neurosubtyping of infants predicts language trajectories. Journal of Neural Transmission. https://doi.org/10.1007/s00702-025-03063-2 derivation details). 

**Four analysis questions are addressed:**

Q1: Do EEG measures primarily associated with ASD including frontal gamma power, functional connectivity across networks (i.e., auditory network, speech network and language network), gamma lateralisation, peak alpha frequency and theta phase consistency differ by sex? 
Q1.2: Are these differences modulated by ASD sex likelihood or ASD diagnosis? (Exploratory)
Q2: Where differences in EEG measures are identified (speech connectivity), does biological sex moderate the relationship between continuous EEG measures and language development from 6 to 36 months? 
Q3: Given that individual EEG metrics may not capture the full complexity of early brain organization, can a re-analysis of previous multivariate clustering results shed light on how a combination of EEG brain metrics relates to sex-differentiated language trajectories?


### Step 4: Publication Figures

**File:** `code/figures/Publication Figures.ipynb`

**What it does:**
This Python (Jupyter) notebook generates all main and supplementary publication-quality figures from the pre-processed dataset. 
 

### Dependencies

#### R

Tested on R 4.3+.

```r
install.packages(c(
  "tidyverse", "lme4", "lmerTest", "performance",
  "emmeans", "effectsize", "flextable", "officer",
  "lmtest", "car", "sandwich"
))
```

#### Python

Tested on Python 3.10+.

```bash
pip install pandas numpy matplotlib seaborn scipy statsmodels pingouin plotly scikit-learn umap-learn
```
 

## Author

Gabriel Blanco-Gomez
