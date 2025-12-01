################################################################################
# Project: mGFR and outcomes in SCREAM
# Purpose: code for sensitivity analysis excluding kidney transplant recipients
# Written by: Antoine Creon
# Date: 2024-11-13
################################################################################

################################################################################
###                              LOAD DATA                                   ###
################################################################################

# Load packages, helper functions and covariate/predictor lists
source(here::here(
   "code",
   "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"
))

# MICE-imputed data
data_imp <- read_rds(
   file = here::here("data", "cleaned", "not_cens_kfrt_imp_cs.rds")
)

# EXCLUDE patients with history of kidney transplantation
data_imp_wo_ktr <- data_imp |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(transplant == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed

# EXCLUDE patients with history of MACE (MI and stroke)
data_imp_wo_mace_ktr <- data_imp_wo_ktr |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(history_mi == 0, history_stroke == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed

# EXCLUDE patients with history of heart failure
data_imp_wo_hf_ktr <- data_imp_wo_ktr |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(hf == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed

# EXCLUDE patients with history of AKI
data_imp_wo_aki_ktr <- data_imp_wo_ktr |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(history_aki == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed


################################################################################
###                           PERFORM THE ANALYSIS                           ###
################################################################################

#  -------------------------------- PLOTS --------------------------------------

# Death and KFRT
death_kfrt_woKTR_plots <- c(`All-cause Death` = "death", "KFRT" = "rrt") %>%
   imap(
      ~ plot_imputed_outcomes(
         .predictor = predictors,
         outcome_element = .x,
         outcome_name = .y,
         .covariates = covariates,
         .imp_data = data_imp_wo_ktr,
         .ref = 90,
         .trunc = 120
      )
   )


# Heart failure excluding patients with history of HF
wo_hf_woKTR_ <- plot_imputed_outcomes(
   .predictor = predictors,
   outcome_element = "heart_failure",
   outcome_name = "Heart Failure",
   .covariates = covariates[!covariates %in% c("hf")],
   .imp_data = data_imp_wo_hf_ktr,
   .ref = 90,
   .trunc = 120
)

# MACE excluding patients with history of stroke or MI
wo_MACE_woKTR_ <- plot_imputed_outcomes(
   .predictor = predictors,
   outcome_element = "MACE_wo_HF",
   outcome_name = "MACE",
   .covariates = covariates[!covariates %in% c("history_stroke", "history_mi")],
   .imp_data = data_imp_wo_mace_ktr,
   .ref = 90,
   .trunc = 120
)

# AKI excluding patients with history of AKI
wo_AKI_woKTR_ <- plot_imputed_outcomes(
   .predictor = predictors,
   outcome_element = "aki",
   outcome_name = "AKI",
   .covariates = covariates[!covariates %in% c("history_aki")],
   .imp_data = data_imp_wo_aki_ktr,
   .ref = 90,
   .trunc = 120
)


# Append all lists
woKTR_plots <- append(
   x = death_kfrt_woKTR_plots,
   list("HF" = wo_hf_woKTR_, "MACE" = wo_MACE_woKTR_, "AKI" = wo_AKI_woKTR_)
)

# Save them as R objects to be called and modified in analysis reports
save(woKTR_plots, file = here::here("output", "r_objects", "woKTR_plots.rda"))


#  ------------------------------- TABLES --------------------------------------
## Death and KFRT
death_kfrt_woKTR_tbl <- c(`All-cause Death` = "death", "KFRT" = "rrt") %>%
   map2(
      .x = .,
      .y = names(.),
      ~ summarize_HR_MICE(
         .predictor = predictors,
         outcome_element = .x,
         outcome_name = {{ .y }},
         .covariates = covariates,
         .imp_data = data_imp_wo_ktr,
         .ref = 90,
         .trunc = 120
      )
   )

## Heart failure excluding patients with history of HF
wo_hf_woKTR_tbl <- summarize_HR_MICE(
   .predictor = predictors,
   outcome_element = "heart_failure",
   outcome_name = "Heart Failure",
   .covariates = covariates[!covariates %in% c("hf")],
   .imp_data = data_imp_wo_hf_ktr,
   .ref = 90,
   .trunc = 120
)

## MACE excluding patients with history of stroke or MI
wo_MACE_woKTR_tbl <- summarize_HR_MICE(
   .predictor = predictors,
   outcome_element = "MACE_wo_HF",
   outcome_name = "MACE",
   .covariates = covariates[!covariates %in% c("history_stroke", "history_mi")],
   .imp_data = data_imp_wo_mace_ktr,
   .ref = 90,
   .trunc = 120
)

## AKI excluding patients with history of AKI
wo_AKI_woKTR_tbl <- summarize_HR_MICE(
   .predictor = predictors,
   outcome_element = "aki",
   outcome_name = "AKI",
   .covariates = covariates[!covariates %in% c("history_aki")],
   .imp_data = data_imp_wo_aki_ktr,
   .ref = 90,
   .trunc = 120
)

# merge all summarized tables
woKTR_tbl <- purrr::reduce(
   .x = append(
      x = death_kfrt_woKTR_tbl,
      list(
         "Heart Failure" = wo_hf_woKTR_tbl,
         "MACE" = wo_MACE_woKTR_tbl,
         "AKI" = wo_AKI_woKTR_tbl
      )
   ),
   .f = function(x, y) {
      left_join(x, y, by = c("x1", "predictor"))
   }
) %>%
   relocate(
      predictor,
      GFR = x1,
      `All-cause Death`,
      KFRT,
      AKI,
      MACE,
      `Heart Failure`
   )


# Save  as R object to be called and modified in analysis reports
save(woKTR_tbl, file = here::here("output", "r_objects", "woKTR_tbl.rda"))
