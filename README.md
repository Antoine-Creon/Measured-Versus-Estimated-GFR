# Measured Versus Estimated GFR

<p align="center">
  <img src="docs/KI_logo.jpg" alt="Project Logo" width="250">
</p>

## Project Overview  
This repository contains the analysis code for the project *[PROJECT NAME]*, led by the [Juan Jesus Carrero group at Karolinska Institutet](https://ki.se/en/research/research-areas-centres-and-networks/research-groups/cardio-renal-epidemiology-juan-jesus-carreros-research-group
) and associated with the manuscript submitted to **[JOURNAL NAME]**.  

### Team 
[Edouard L. Fu](https://edouard-fu.github.io/)<sup>1,2*</sup>, [Antoine Créon](https://ki.se/en/people/antoine-creon)<sup>1\*</sup>, Josef Coresh<sup>3,4</sup>, Morgan E. Grams<sup>5</sup>, Michael G. Shlipak<sup>6</sup>, Lesley A. Inker<sup>7</sup>, Andrew S. Levey<sup>7</sup>, [Juan-Jesus Carrero](https://ki.se/en/people/juan-jesus-carrero)<sup>1,8</sup>

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
6,059 participants aged 18 years or older with plasma iohexol testing and concurrent measurements of creatinine and cystatin C.

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
- `smcfcs` 1.9.0
- `rms` 6.8.2
- `Hmisc` 5.2.0
- `survival` 3.7.0
- `ragg` 1.3.3
- `ggokabeito` 0.1.0
- `gt` 0.11.1

## Description of scripts (which script does what).

- Data preparation and cleaning
  - `01_covariates-naming.R`: labeling of the variables. Definitions of the variables are provided in the supplementary material.
  - `02_missing-data-imputation.R`: missing data description, and imputation using MICE or SMC-FCS.
- Main analysis
  - `03_helper-functions-survival.R`: helper functions used in the survival analyses
  - `04_main-analysis_CKDEPI2021.R`: mGFR and eGFR versus health outcomes after imputation by MICE
- Supporting analyses:
  - `05_supporting-analyses_EKFC_CKDEPI2009-12.R`: main analysis but using the EKFC and CKD-EPI 2009-2012 equations
- Sensitivity analyses:
  - `06_sensitivity-analysis_without-UACR-enrichment.R`: analysis without converting PCR and dipstick results to UACR.
  - `07_sensitivity-analysis_without-KTR.R`: analysis after exclusion of kidney transplant recipients
  - `08_sensitivity-analysis_prevalent-HF_incident-MACE-AKI.R`: analysis after 
    - exclusion of individuals with history of AKI for the outcome AKI
    - exclusion of individuals with history of MACE for the outcome MACE
    - inclusion of individuals with history of heart failure for the outcome hospitalization with heart failure
  - `09_sensitivity-analysis_SMC-FCS.R`: main analysis after multiple imputation by SMC-FCS
  - `10_sensitivity-analyses_complete-cases.R`: complete cases analysis


## Data sharing statement
The repository does not contain the raw data used in the analysis. Data will be available for collaborative research under reasonable request and fulfillment of GDPR regulations. For inquiries, please send your proposal to [the Steering Committee of the SCREAM project](juan.jesus.carrero@ki.se).


