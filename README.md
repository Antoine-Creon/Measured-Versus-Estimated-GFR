# Measured Versus Estimated GFR

<p align="center">
  <img src="docs/KI_logo.jpg" alt="Project Logo" width="250">
</p>

## Project Overview  
This repository contains the analysis code for the project *[PROJECT NAME]*, led by the [Juan Jesus Carrero group at Karolinska Institutet](https://ki.se/en/research/research-areas-centres-and-networks/research-groups/cardio-renal-epidemiology-juan-jesus-carreros-research-group
) and associated with the manuscript submitted to **[JOURNAL NAME]**.  

### Team 
[Edouard L. Fu](https://edouard-fu.github.io/)<sup>1,2*</sup>, [Antoine Créon](https://antoine-creon.github.io/)<sup>1\*</sup>, Josef Coresh<sup>3,4</sup>, Morgan E. Grams<sup>5</sup>, Michael G. Shlipak<sup>6</sup>, Lesley A. Inker<sup>7</sup>, Andrew S. Levey<sup>7</sup>, [Juan-Jesus Carrero](https://ki.se/en/people/juan-jesus-carrero)<sup>1,8</sup>

- <sup>1</sup> Department of Medical Epidemiology and Biostatistics, Karolinska Institute, Stockholm, Sweden
- <sup>2</sup> Department of Clinical Epidemiology, Leiden University Medical Center, Leiden, the Netherlands
- <sup>3</sup> Optimal Aging Institute and Division of Epidemiology, Department of Population Health, New York University Grossman School of Medicine, New York, New York
- <sup>4</sup> Department of Epidemiology, Johns Hopkins University Bloomberg School of Public Health, Baltimore, Maryland
- <sup>5</sup> Division of Precision Medicine, Department of Medicine, New York University Grossman School of Medicine, New York, New York
- <sup>6</sup> Kidney Health Research Collaborative, Department of Medicine, San Francisco Veterans Affairs Health Care System and University of California San Francisco, San Francisco, California, USA
- <sup>7</sup> Division of Nephrology, Department of Internal Medicine, Tufts Medical Center, Boston, Massachusetts
- <sup>8</sup> Division of Nephrology, Department of Clinical Sciences, Karolinska Institute, Danderyd Hospital, Stockholm, Sweden

<sup>\*</sup> These authors contributed equally to this work as co-first authors

Correspondence: [Dr. Edouard L. Fu](e.l.fu@lumc.n)


### Objectives  
To quantify associations between mGFR and a range of adverse clinical outcomes, and to investigate which eGFR most closely resembled these associations.

### Design
Cohort study using routinely collected healthcare data

### Setting
Stockholm, Sweden, January 2011 to December 2021

### Participants
6,174 participants aged 18 years or older with plasma iohexol testing and concurrent measurements of creatinine and cystatin C.

### Exposures

mGFR, eGFR<sub>cr</sub>, eGFR<sub>cys</sub> and eGFR<sub>cr-cys</sub>.

## Installation & Dependencies  
### Operating system 
Windows 10 Enterprise Version	22H2

### Software
#### Programming language
All analyses were performed using R, version 4.4.2 (2024-10-31 ucrt) *Pile of Leaves*.

#### IDE
IDE RStudio 2024.12.1+563 *Kousa Dogwood* for windows

### R (meta)packages

- `qs` 0.27.2
- `here` 1.0.1
- `tidyverse` 2.0.0
- `labelled` 2.13.0
- `glue` 1.8.0
- `mice` 3.16.0
- `micemd` 1.10.0
- `mitools` 2.4
- `rms` 6.8.2
- `Hmisc` 5.2.0
- `survival` 3.7.0
- `ragg` 1.3.3
- `ggokabeito` 0.1.0
- `gt` 0.11.1

## Description of scripts

| Script | Description |
|--------|-------------|
| `01_missing-data-imputation.R` | Performs multiple imputation using MICE (Multivariate Imputation by Chained Equations) for missing UACR and BMI values. Includes cubic spline transformations and Nelson-Aalen cumulative hazard estimates as auxiliary variables. |
| `02_helper-functions-survival.R` | Defines helper functions for eGFR calculations (CKD-EPI 2009/2021, EKFC), Cox proportional hazards model fitting, term plot extraction, Rubin's rules pooling for multiply imputed data, and visualization of hazard ratios. |
| `03_analysis-preparation.R` | Prepares the analysis environment by loading packages, defining predictors, covariates, and outcome events (death, KFRT, AKI, MACE, heart failure). |
| `04_main-analysis_CKDEPI2021.R` | Conducts the main survival analysis using CKD-EPI 2021 equations. Fits Cox models on multiply imputed data, generates hazard ratio plots and summary tables for all outcomes. Excludes patients with prior events for incident outcome analyses. |
| `05_supporting-analyses_EKFC_CKDEPI2009-12.R` | Performs supporting analyses using alternative eGFR equations (EKFC and CKD-EPI 2009/2012). Generates distribution plots and hazard ratio estimates for comparison with CKD-EPI 2021 results. |
| `06_sensitivity-analysis_without-UACR-enrichment.R` | Sensitivity analysis excluding the conversion of protein-to-creatinine ratio (PCR) and dipstick results into UACR values. Tests robustness of findings to UACR enrichment assumptions. |
| `07_sensitivity-analysis_without-KTR.R` | Sensitivity analysis excluding kidney transplant recipients (KTR) from the study population. Evaluates whether results are driven by this subgroup with unique clinical characteristics. |
| `08_conditional_incidence-rates.R` | Computes conditional incidence rates by eGFR category using Poisson regression with pooling across multiply imputed datasets. Estimates adjusted incidence rates at median covariate values. |


## Data sharing statement
The repository does not contain the raw data used in the analysis. Data will be available for collaborative research under reasonable request and fulfillment of GDPR regulations. For inquiries, please send your proposal to [the Steering Committee of the SCREAM project](juan.jesus.carrero@ki.se).


