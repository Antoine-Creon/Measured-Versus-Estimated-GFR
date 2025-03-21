################################################################################
# Project: mGFR and outcomes in SCREAM 
# Purpose: Code for naming of covariates
# Written by: Antoine Creon
# Date: 2024-11-02
################################################################################

################################################################################
###                         LOAD DATA AND PACKAGES                           ###
################################################################################

pacman::p_load(tidyverse,
               here,
               labelled,
               qs)

# Dataset with outcomes not censored at KFRT
data_uncens <- qread(file = here::here("data","cleaned","data_not_cens_kfrt.qs"))

# Define antihypertensives and factor mGFR and BMI
data_uncens <- data_uncens %>%
   dplyr::mutate(bmi_cat = factor(bmi_cat, 
                                  levels = 0:4, 
                                  labels = c("Missing", "<20", "20 to <25", "25 to <30", ">=30")),
                 mgfr_category = factor(mgfr_category),
                 antihypertensives = if_else(bblock == 1 | ccb == 1 | diur == 1 | rasi == 1, 1, 0))

                                                
################################################################################
###           NAME VARIABLES TO BE DESCRIBED, REMOVE UNUSED ONES             ###
################################################################################

var_names <- list(age = "Age",
                  age_cat = "Age group",
                  female = "Female",
                  bmi = "Body Mass Index",
                  bmi_cat = "Body Mass Index Group",
                  creat = "Creatinine",
                  cystatin = "Cystatin C",
                  ckd_epi_2021_cr = "eGFRcr, ml/min/1.73m2",
                  ckd_epi_2012_cys = "eGFRcys, ml/min/1.73m2", 
                  ckd_epi_2021_cr_cys = "eGFRcr-cys, ml/min/1.73m2",
                  mgfr = "mGFR, ml/min/1.73m2",
                  mgfr_category = "mGFR category",
                  uacr_enriched = "Urine Albumin-to-Creatinine Ratio (mg/mmol)",
                  albu_categ = "Urine Albumin-to-Creatinine Ratio category",
                  uacr_enriched_missing_indicator = "UACR missing",
                  hyperten = "Hypertension",
                  history_mi = "Myocardial infarction",
                  ihd = "Other ischemic heart disease",
                  hf= "Heart failure",
                  history_stroke = "Stroke",
                  cevd = "Other cerebrovascular disease",
                  arrh = "Arrhythmia (including atrial fibrillation)",
                  history_AF = "Atrial fibrillation",
                  pvd = "Peripheral vascular disease",
                  dm = "Diabetes mellitus",
                  cancer = "Cancer in previous year",
                  copd = "Chronic obstructive pulmonary disease",
                  liver = "Liver disease",
                  history_aki = "Acute kidney injury",
                  transplant = "Kidney Transplant Recipient",
                  kidney_donor = "Kidney donor",
                  glucocorticosteroids = "Glucocorticosteroids",
                  bblock = "Beta blocker",
                  ccb = "Calcium channel blocker",
                  diur = "Diuretic",
                  rasi = "ACEi/ARB",
                  lipid = "Lipid lowering drug",
                  nsaid = "NSAID",
                  antihypertensives = "Antihypertensive drugs",
                  calendar_year = "Calendar Year"
                  )

# Save variable names
write_rds(var_names, file = here::here("data","cleaned","variable_names.rds"))

# Covariates to include in models
covariates <-   c("age", "female", "bmi", "history_mi",        "hyperten", "hf", "history_stroke", "history_AF",   "pvd", "dm", "copd", "cancer", "liver",                "transplant", "antihypertensives",             "lipid", "nsaid", "log_uacr_enriched", "glucocorticosteroids") 
write_rds(covariates, file = here::here("data","cleaned","covariates_names.rds"))

# Events
events <- c("death", "rrt", "aki", "heart_failure", "MACE_wo_HF")

# predictors of interest
predictors <- list("mgfr", "ckd_epi_2021_cr", "ckd_epi_2012_cys", "ckd_epi_2021_cr_cys")

# Name variables and keep those of interest
data_uncens <- data_uncens %>% 
   select(names(var_names), 
          starts_with("time_to"),
          starts_with("history_"),
          all_of(events),
          starts_with("log_")) %>%
   set_variable_labels(.labels = var_names)

write_rds(data_uncens, file = here::here("data","cleaned","data_not_cens_kfrt_named.rds"))

