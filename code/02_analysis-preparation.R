################################################################################
## mGFR and outcomes in SCREAM ##
# Code for main analysis GFR versus outcomes
# Written by: Antoine Creon
# Date: 2024-11-21
################################################################################

################################################################################
# LOAD FUNCTIONS ###############################################################
################################################################################

pacman::p_load(
  tidyverse,
  here,
  ggokabeito,
  qs,
  Hmisc,
  survival,
  rms,
  mice,
  mitools
)

# Helper functions for the survival analysis
source(here::here("code", "01_helper-functions.R"))


################################################################################
# PREPARE THE MAIN ANALYSIS ####################################################
################################################################################

## VECTOR OF PREDICTORS AND COVARIATES -----------------------------------------

# List the predictors
predictors <- list(
  "mgfr",
  "ckd_epi_2021_cr",
  "ckd_epi_2012_cys",
  "ckd_epi_2021_cr_cys"
)

# Create a vector of the confounders
## UACR not imputed
covariates <- c(
  "age",
  "female",
  "bmi",
  "history_mi",
  "hyperten",
  "hf",
  "history_stroke",
  "history_AF",
  "pvd",
  "dm",
  "copd",
  "cancer",
  "liver",
  "transplant",
  "antihypertensives",
  "lipid",
  "nsaid",
  "log_uacr_enriched",
  "glucocorticosteroids"
)

# List the events and give them a proper name
events <- list(
  "All-cause Death" = "death",
  "KFRT" = "rrt",
  "AKI" = "aki",
  "MACE" = "MACE_wo_HF",
  "Heart Failure" = "heart_failure"
)
