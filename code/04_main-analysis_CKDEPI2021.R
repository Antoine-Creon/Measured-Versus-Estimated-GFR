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

source(here::here(
  "code",
  "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"
))

data_imp <- read_rds(
  file = here::here("data", "cleaned", "not_cens_kfrt_imp_cs.rds")
)

# EXCLUDE patients with history of MACE (MI and stroke)
data_imp_wo_mace <- data_imp |>
  complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
  filter(history_mi == 0, history_stroke == 0) |> # Apply filter to each dataset
  mice::as.mids() # reconstruct into a mids object if needed

# EXCLUDE patients with history of heart failure
data_imp_wo_hf <- data_imp |>
  complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
  filter(hf == 0) |> # Apply filter to each dataset
  mice::as.mids() # reconstruct into a mids object if needed

# EXCLUDE patients with history of AKI
data_imp_wo_aki <- data_imp |>
  complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
  filter(history_aki == 0) |> # Apply filter to each dataset
  mice::as.mids() # reconstruct into a mids object if needed


################################################################################
###                              PERFORM ANALYSES                            ###
################################################################################

#  --------------------------------- Plots -------------------------------------

# Prevalent patients (for all-cause mortality and KFRT)
plots_ref90_imp <- events %>%
  imap(
    ~ plot_imputed_outcomes(
      .predictor = predictors,
      outcome_element = .x,
      outcome_name = .y,
      .covariates = covariates,
      .imp_data = data_imp,
      .ref = 90,
      .trunc = 120
    )
  )

# Heart failure excluding patients with history of HF
wo_hf <- plot_imputed_outcomes(
  .predictor = predictors,
  outcome_element = "heart_failure",
  outcome_name = "Heart failure",
  .covariates = covariates[!covariates %in% c("hf")],
  .imp_data = data_imp_wo_hf,
  .ref = 90,
  .trunc = 120
)

# MACE excluding patients with history of stroke or MI
wo_MACE <- plot_imputed_outcomes(
  .predictor = predictors,
  outcome_element = "MACE_wo_HF",
  outcome_name = "MACE",
  .covariates = covariates[!covariates %in% c("history_stroke", "history_mi")],
  .imp_data = data_imp_wo_mace,
  .ref = 90,
  .trunc = 120
)

# AKI excluding patients with history of AKI
wo_AKI <- plot_imputed_outcomes(
  .predictor = predictors,
  outcome_element = "aki",
  outcome_name = "AKI",
  .covariates = covariates[!covariates %in% c("history_aki")],
  .imp_data = data_imp_wo_aki,
  .ref = 90,
  .trunc = 120
)


#  -------------------------------- Tables -------------------------------------

# Prevalent patients (for all-cause mortality and KFRT)
tbl_ref90_imp <- events %>%
  map2(
    .x = .,
    .y = names(.),
    ~ summarize_HR_MICE(
      .predictor = predictors,
      outcome_element = .x,
      outcome_name = {{ .y }},
      .covariates = covariates,
      .imp_data = data_imp,
      .ref = 90,
      .trunc = 120
    )
  )

# Heart failure excluding patients with history of HF
wo_hf_tbl <- summarize_HR_MICE(
  .predictor = predictors,
  outcome_element = "heart_failure",
  outcome_name = "Heart failure",
  .covariates = covariates[!covariates %in% c("hf")],
  .imp_data = data_imp_wo_hf,
  .ref = 90,
  .trunc = 120
)

# MACE excluding patients with history of stroke or MI
wo_MACE_tbl <- summarize_HR_MICE(
  .predictor = predictors,
  outcome_element = "MACE_wo_HF",
  outcome_name = "MACE",
  .covariates = covariates[!covariates %in% c("history_stroke", "history_mi")],
  .imp_data = data_imp_wo_mace,
  .ref = 90,
  .trunc = 120
)

# AKI excluding patients with history of AKI
wo_AKI_tbl <- summarize_HR_MICE(
  .predictor = predictors,
  outcome_element = "aki",
  outcome_name = "AKI",
  .covariates = covariates[!covariates %in% c("history_aki")],
  .imp_data = data_imp_wo_aki,
  .ref = 90,
  .trunc = 120
)
