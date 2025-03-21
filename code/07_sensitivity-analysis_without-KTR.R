################################################################################
# Project: mGFR and outcomes in SCREAM 
# Purpose: code for sensitivity analysis excluding kidney transplant recipients
# Written by: Antoine Creon
# Date: 2024-11-13
################################################################################

################################################################################
###                 Load imputed dataset and helper functions                ###
################################################################################

source(here::here("code", "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"))

data_imp <- read_rds(file = here::here("data","cleaned","not_cens_kfrt_imp_cs.rds"))

# FILTER KTR
data_imp_wo_ktr <- data_imp |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(transplant == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed

# FILTER HF and KTR 
data_imp_wo_hf_ktr <- data_imp_wo_ktr |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(hf == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed


################################################################################
###                              CUBIC SPLINES                               ###
################################################################################

plots_ref90_imp_wo_ktr <- events[!events %in% c("heart_failure")] %>%
   imap(~ plot_imputed_outcomes(.predictor = predictors, 
                                outcome_element = .x, 
                                outcome_name = .y, 
                                .covariates = covariates[!covariates %in% c("transplant")], 
                                .imp_data = data_imp_wo_ktr, 
                                .ref = 90, 
                                .trunc = 120)
   )

hf_imp_wo_ktr <- plot_imputed_outcomes(.predictor = predictors, 
                                       outcome_element = "heart_failure", 
                                       outcome_name = "Heart failure", 
                                       .covariates =  covariates[!covariates %in% c("transplant", "hf")], 
                                       .imp_data = data_imp_wo_hf_ktr, 
                                       .ref = 90, 
                                       .trunc = 120) 


MICE_withoutKTR_cs <- append(x = plots_ref90_imp_wo_ktr, list("HF" = hf_imp_wo_ktr ))
save(MICE_withoutKTR_cs, file = here::here("output","r_objects","gfr-vs-outcomes_MICE_withoutKTR_cs.rda"))


#  -------------------- Make a table with the main results ---------------------

tbl_ref90_imp_wo_ktr <- events[!events %in% c("heart_failure")] %>%
   map2(.x = .,
        .y = names(.),
        ~ summarize_HR_MICE(.predictor = predictors, 
                            outcome_element = .x, 
                            outcome_name = {{.y}}, 
                            .covariates = covariates[!covariates %in% c("transplant")], 
                            .imp_data = data_imp_wo_ktr, 
                            .ref = 90, 
                            .trunc = 120)
   )

tbl_hf_imp_wo_ktr <- summarize_HR_MICE(.predictor = predictors, 
                                       outcome_element = "heart_failure", 
                                       outcome_name = "Heart Failure", 
                                       .covariates = covariates[!covariates %in% c("hf", "transplant")], 
                                       .imp_data = data_imp_wo_hf_ktr, 
                                       .ref = 90, 
                                       .trunc = 120)

MICE_withoutKTR_cs_tbl <- purrr::reduce(.x = append(x = tbl_ref90_imp_wo_ktr, list("HF" = tbl_hf_imp_wo_ktr)),
                                        .f = function(x, y) {left_join(x, y, by = c("x1", "predictor"))}) %>%
   relocate(predictor, GFR = x1, `All-cause Death`,  KFRT, AKI, `Heart Failure`, MACE)

save(MICE_withoutKTR_cs_tbl, file = here::here("output","r_objects","gfr-vs-outcomes_MICE_withoutKTR_cs_tbl.rda"))


