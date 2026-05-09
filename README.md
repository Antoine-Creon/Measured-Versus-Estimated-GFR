# Measured Versus Estimated GFR

<p align="center">
  <img src="docs/KI_logo.jpg" alt="Project Logo" width="250">
</p>

## Project Overview  
This repository contains the analysis code for the project _Measured and estimated glomerular filtration rates and risk of adverse health outcomes_, published in **[JOURNAL NAME]**.  

## Team 

This study was led by [Edouard L. Fu](https://edouard-fu.github.io/)<sup>1,2*</sup>, [Antoine Créon](https://antoine-creon.github.io/)<sup>2\*</sup> and [Juan Jesus Carrero](https://ki.se/en/research/research-areas-centres-and-networks/research-groups/cardio-renal-epidemiology-juan-jesus-carreros-research-group)<sup>2,10</sup> as part of an international collaboration including Morgan E. Grams<sup>3</sup>, Josef Coresh<sup>4,5</sup>, Arvid Sjölander<sup>2</sup>, Anne-Laure Faucon<sup>2</sup>, Michelle M. Estrella<sup>6,7,8</sup>, Friedo W. Dekker<sup>1</sup>, Michael G. Shlipak<sup>8</sup>, Lesley A. Inker<sup>9</sup> and Andrew S. Levey<sup>9</sup>.

- <sup>1</sup> Department of Clinical Epidemiology, Leiden University Medical Center, Leiden, the Netherlands
- <sup>2</sup> Department of Medical Epidemiology and Biostatistics, Karolinska Institute, Stockholm, Sweden
- <sup>3</sup> Division of Precision Medicine, Department of Medicine, New York University Grossman School of Medicine, New York, New York
- <sup>4</sup> Optimal Aging Institute and Division of Epidemiology, Department of Population Health, New York University Grossman School of Medicine, New York, New York
- <sup>5</sup> Department of Epidemiology, Johns Hopkins University Bloomberg School of Public Health, Baltimore, Maryland
- <sup>6</sup> Division of Nephrology, Department of Medicine, University of California, San Francisco, San Francisco, California, USA
- <sup>7</sup> Renal Section, Medical Service, San Francisco VA Health Care System, San Francisco, California, USA
- <sup>8</sup> Kidney Health Research Collaborative, Department of Medicine, San Francisco Veterans Affairs Health Care System and University of California San Francisco, San Francisco, California, USA
- <sup>9</sup> Division of Nephrology, Department of Medicine, Tufts Medical Center, Boston, Massachusetts
- <sup>10</sup> Division of Nephrology, Department of Clinical Sciences, Karolinska Institute, Danderyd Hospital, Stockholm, Sweden

<sup>\*</sup> These authors contributed equally to this work as co-first authors

Corresponding author: [Dr. Edouard L. Fu](mailto:e.l.fu@lumc.nl)

## Installation & Dependencies  
### Operating system 
Windows 10 Enterprise Version	22H2

### Software
Analyses were run under R 4.4.2 and repeated under R 4.5.2 with consistent results.

#### R (meta)packages

- `tidyverse` 2.0.0
- `here` 1.0.2
- `qs` 0.27.3
- `Hmisc` 5.2-4
- `survival` 3.8-3
- `rms` 8.1-0
- `mice` 3.19.0
- `mitools` 2.4
- `ggokabeito` 0.1.0
- `glue` 1.8.0
- `scales` 1.4.0
- `rlang` 1.1.6
- `splines` 4.5.2
- `parallelly` 1.46.0
- `boot` 1.3-32
- `conflicted` 1.2.0
- `pacman` 0.5.1

## Description of scripts

| Script | Description |
|--------|-------------|
| `01_helper-functions.R` | Helper functions for eGFR equations, Cox models, termplot extraction, pooling across imputations, plotting, and conditional incidence rates. |
| `02_analysis-preparation.R` | Loads packages and defines predictors, covariates, and outcome events. |
| `03_missing-data-imputation.R` | Multiple imputation for missing UACR and BMI using MICE, including splines and Nelson-Aalen estimates. |
| `04_main-analysis_CKDEPI2021.R` | Main survival analysis using CKD-EPI 2021; plots and tables for all outcomes; incident analyses exclude prior events. |
| `05_supporting-analyses_EKFC_CKDEPI2009-12.R` | Supporting analyses with EKFC and CKD-EPI 2009/2012 equations; plots and tables. |
| `06_sensitivity-analysis_without-UACR-enrichment.R` | Sensitivity analysis without UACR enrichment (no PCR/dipstick conversion). |
| `07_sensitivity-analysis_without-KTR.R` | Sensitivity analysis excluding kidney transplant recipients. |
| `08_sensitivity-analysis_mortality-in-full-creat-pop.R` | Mortality analysis in full creatinine population vs study population; conditional IRs; hazard ratio contrasts. |
| `09_conditional_incidence-rates.R` | Conditional incidence rates for continuous GFR across outcomes; includes CKD-EPI 2009, EKFC, and sensitivity variants. |
| `10_ratios-of-hazard-ratios.R` | Bootstrap + multiple imputation to estimate ratios of hazard ratios (BOOT-MI). |

## Data sharing statement
The repository does not contain the raw data used in the analysis. Data will be available for collaborative research under reasonable request and fulfillment of GDPR regulations. For inquiries, please send your proposal to [the Steering Committee of the SCREAM project](juan.jesus.carrero@ki.se).


