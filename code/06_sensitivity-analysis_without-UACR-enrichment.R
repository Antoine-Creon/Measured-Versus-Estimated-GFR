################################################################################
# Project: mGFR and outcomes in SCREAM 
# Purpose: code for sensitivity analysis, without converting PCR/disptick results
#          into UACR
# Written by: Antoine Creon
# Date: 2025-03-13
################################################################################

################################################################################
###                     Load data and helper functions                       ###
################################################################################

source(here::here("code", "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"))

wo_uacr <- read_rds(file = here::here("data","cleaned","MICE_UACR-not-enriched.rds")) # imputation performed in dedicated script

# FILTER patients with history of heart failure
wo_uacr_wo_hf <- wo_uacr |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(hf == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed



################################################################################
###                             Perform analysis                             ###
################################################################################


wo_converted_UACR <- events[!events %in% c("heart_failure")] %>%
   imap(~ plot_imputed_outcomes(.predictor = predictors, 
                                outcome_element = .x, 
                                outcome_name = .y, 
                                .covariates = covariates, 
                                .imp_data = wo_uacr, 
                                .ref = 90, 
                                .trunc = 120)
   )

wo_converted_UACR_HF <- plot_imputed_outcomes(.predictor = predictors, 
                                outcome_element = "heart_failure", 
                                outcome_name = "Heart failure", 
                                .covariates = covariates[!covariates %in% c("hf")], 
                                .imp_data = wo_uacr_wo_hf, 
                                .ref = 90, 
                                .trunc = 120) 

wo_UACR_plots <- append(x = wo_converted_UACR, list("HF" = wo_converted_UACR_HF))
save(wo_UACR_plots, file = here::here("output","r_objects","MICE_wo-UACR-conversion.rda"))

# Tables with HRs
wo_converted_UACR_tbl <- events[!events %in% c("heart_failure")] %>%
   map2(.x = .,
        .y = names(.),
        ~ summarize_HR_MICE(.predictor = predictors, 
                            outcome_element = .x, 
                            outcome_name = {{.y}}, 
                            .covariates = covariates, 
                            .imp_data = wo_uacr, 
                            .ref = 90, 
                            .trunc = 120)
   )

wo_converted_UACR_HF_tbl <- summarize_HR_MICE(.predictor = predictors, 
                                outcome_element = "heart_failure", 
                                outcome_name = "Heart Failure", 
                                .covariates = covariates[!covariates %in% c("hf")], 
                                .imp_data = wo_uacr_wo_hf, 
                                .ref = 90, 
                                .trunc = 120)

wo_UACR_tbl <- purrr::reduce(.x = append(x = wo_converted_UACR_tbl, list("HF" = wo_converted_UACR_HF_tbl)),
                                     .f = function(x, y) {left_join(x, y, by = c("x1", "predictor"))}) %>%
   relocate(predictor, GFR = x1, `All-cause Death`,  KFRT, AKI, `Heart Failure`, MACE)

save(wo_UACR_tbl, file = here::here("output","r_objects","MICE_wo-UACR-conversion_tbl.rda"))
