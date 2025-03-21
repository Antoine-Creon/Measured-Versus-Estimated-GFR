################################################################################
# Project: mGFR and outcomes in SCREAM 
# Purpose: code for complete cases analysis of m/eGFR versus outcomes
# Written by: Antoine Creon
# Date: 2024-11-13
################################################################################

################################################################################
###                        LOAD DATA AND FUNCTIONS                           ###
################################################################################

source(here::here("code", "05_helper-functions-survival.R"))
source(here::here("code", "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"))

data <- read_rds(file = here::here("data","cleaned","data_not_cens_kfrt_named.rds"))

# Dataset without prevalent HF patients
data_wo_hf <- data %>%
   filter(hf == 0)


################################################################################
###                                  ANALYSIS                                ###
################################################################################

#  --------------------------- PLOTS CUBIC SPLINES -----------------------------

# Outcomes with ref 90 and truncation 120 ml/min
plots_ref90 <- events[!events %in% c( "heart_failure")] %>%
   imap(~ plot_outcomes(outcome_element = .x, outcome_clean = .y, 
                        .predictor = unlist(predictors),
                        .covariates = covariates,
                        .data = data, .ref = 90, .trunc = 120))
   
hf <- predictors %>%
   map(~ compute_fits("heart_failure", 
                      .x, 
                      covariates[!covariates %in% c("hf")],
                      data_wo_hf)) %>%
   map(~ extract_termplot(.x, ref = 90, 120)) %>%
   list_rbind() %>%
   plot_HRs(., .outcome = "HF (incident patients)")


#  ---------------- Make a single list of all plots and save -------------------

# Append lists
CCA_uncens_kfrt <- append(x = plots_ref90, list("HF" = hf))

# Save as objects to call them in analysis reports
save(CCA_uncens_kfrt, file = here::here("output","r_objects","gfr-vs-outcomes_CCA_withKTR_cs_plots.rda"))


#  -------------------- Make a table with the main results ---------------------
# a column for GFR level and each outcome
# a row for 15-30-45-45-60 ml/min 
# 4 groups of rows cr, cys, cys-cr and mGFR

# A. extract the summarized HRs for each outcome and predictor 
tbl_ref90 <- events[!events %in% c("heart_failure")] %>%
   imap(~ summarize_HR_CCA(outcome_element = .x, outcome_name = {{.y}}, 
                           .predictor = unlist(predictors),
                           .covariates = covariates,
                           .data = data, .ref = 90, .trunc = 120))

tbl_hf <- predictors %>%
   map(~ compute_fits("heart_failure", 
                      .x, 
                      covariates[!covariates %in% c("hf")],
                      data_wo_hf)) %>%
   map(~ extract_termplot(.x, ref = 90, 120)) %>%
   list_rbind() %>%
   format_table(., .outcome = "Heart failure")


# B. merge all summarized tables
CCA_HRs <- purrr::reduce(.x = append(x = tbl_ref90, list("Heart failure" = tbl_hf)),
                          .f = function(x, y) {left_join(x, y, by = c("x1", "predictor"))}) %>%
   relocate(predictor, GFR = x1, `All-cause Death`,  KFRT, AKI, MACE, `Heart failure`)
            

# C. Save  as R object to be called and modified in analysis reports
save(CCA_HRs, file = here::here("output","r_objects","gfr-vs-outcomes_CCA_withKTR_cs_tbl.rda"))