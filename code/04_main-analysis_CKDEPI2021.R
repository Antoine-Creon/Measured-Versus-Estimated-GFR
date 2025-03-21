################################################################################
# Project: mGFR and outcomes in SCREAM 
# Purpose: code for main analysis of m/eGFR versus outcomes, after MICE
# Written by: Antoine Creon
# Date: 2024-11-13
################################################################################

# steps to follow for each GFR estimating equation or mGFR:
# - compute coxph fit in each imputed dataset
# - extract the termplot in each imputed dataset
# - pool the coefficients (beta and their se) using Rubin's rules 
# - Merge results for all GFR methods to plot/tabulate them
# Perform these steps for each studied outcome

################################################################################
###                 Load imputed dataset and helper functions                ###
################################################################################

source(here::here("code", "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"))

data_imp <- read_rds(file = here::here("data","cleaned","not_cens_kfrt_imp_cs.rds"))

# FILTER patients with history of heart failure
data_imp_wo_hf <- data_imp |>
  complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
  filter(hf == 0) |> # Apply filter to each dataset
  mice::as.mids() # reconstruct into a mids object if needed


################################################################################
###                             CUBIC SPLINES                                ###
################################################################################

#  ------------------------------ Prevalent patients ---------------------------

plots_ref90_imp <- events[!events %in% c("heart_failure")] %>%
  imap(~ plot_imputed_outcomes(.predictor = predictors, 
                               outcome_element = .x, 
                               outcome_name = .y, 
                               .covariates = covariates, 
                               .imp_data = data_imp, 
                               .ref = 90, 
                               .trunc = 120)
  )

#  --------------------------- Incident patients -------------------------------

hf_imp <- plot_imputed_outcomes(.predictor = predictors, 
                                outcome_element = "heart_failure", 
                                outcome_name = "Heart failure", 
                                .covariates = covariates[!covariates %in% c("hf")], 
                                .imp_data = data_imp_wo_hf, 
                                .ref = 90, 
                                .trunc = 120) 


#  ---------------- Make a single list of all plots and save -------------------

MICE_withKTR_cs <- append(x = plots_ref90_imp, list("HF" = hf_imp))

save(MICE_withKTR_cs, file = here::here("output","r_objects","gfr-vs-outcomes_MICE_withKTR_cs.rda"))


#  -------------------- Make a table with the main results ---------------------
# a column for GFR level and each outcome
# a row for 15-30-45-45-60 ml/min 
# 4 groups of rows cr, cys, cys-cr and mGFR

# A. extract the summarized HRs for each outcome and predictor 
tbl_ref90_imp <- events[!events %in% c("heart_failure")] %>%
  map2(.x = .,
       .y = names(.),
       ~ summarize_HR_MICE(.predictor = predictors, 
                           outcome_element = .x, 
                           outcome_name = {{.y}}, 
                           .covariates = covariates, 
                           .imp_data = data_imp, 
                           .ref = 90, 
                           .trunc = 120)
  )

tbl_hf_imp <- summarize_HR_MICE(.predictor = predictors, 
                                outcome_element = "heart_failure", 
                                outcome_name = "Heart Failure", 
                                .covariates = covariates[!covariates %in% c("hf")], 
                                .imp_data = data_imp_wo_hf, 
                                .ref = 90, 
                                .trunc = 120)


# B. merge all summarized tables
MICE_withKTR_cs_tbl <- purrr::reduce(.x = append(x = tbl_ref90_imp, list("HF" = tbl_hf_imp)),
                                     .f = function(x, y) {left_join(x, y, by = c("x1", "predictor"))}) %>%
  relocate(predictor, GFR = x1, `All-cause Death`,  KFRT, AKI, `Heart Failure`, MACE)

# C. Save  as R object to be called and modified in analysis reports
save(MICE_withKTR_cs_tbl, file = here::here("output","r_objects","gfr-vs-outcomes_MICE_withKTR_cs_tbl.rda"))


