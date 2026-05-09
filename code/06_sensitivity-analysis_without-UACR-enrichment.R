################################################################################
# Project: mGFR and outcomes in SCREAM
# Purpose: code for sensitivity analysis, without converting PCR/disptick results
#          into UACR
# Written by: Antoine Creon
# Date: 2025-03-13
################################################################################

################################################################################
# Load data and helper functions ###############################################
################################################################################

source(here::here("code", "02_analysis-preparation.R"))

wo_uacr <- read_rds(
   file = here::here("data", "cleaned", "MICE_UACR-not-enriched.rds")
)

## Prepare datasets ------------------------------------------------------------

# EXCLUDE patients with history of MACE (MI and stroke)
data_imp_wo_mace <- wo_uacr |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(history_mi == 0, history_stroke == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed

# EXCLUDE patients with history of heart failure
data_imp_wo_hf <- wo_uacr |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(hf == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed

# EXCLUDE patients with history of AKI
data_imp_wo_aki <- wo_uacr |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(history_aki == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed

################################################################################
# Perform analysis #############################################################
################################################################################

## Plots -----------------------------------------------------------------------

death_kfrt_wo_uacr <- c(`All-cause Death` = "death", "KFRT" = "rrt") %>%
   imap(
      ~ plot_imputed_outcomes(
         .predictor = predictors,
         outcome_element = .x,
         outcome_name = .y,
         .covariates = covariates,
         .imp_data = wo_uacr,
         .ref = 90,
         .trunc = 120
      )
   )


# Heart failure excluding patients with history of HF
wo_hf_wo_uacr <- plot_imputed_outcomes(
   .predictor = predictors,
   outcome_element = "heart_failure",
   outcome_name = "Heart failure",
   .covariates = covariates[!covariates %in% c("hf")],
   .imp_data = data_imp_wo_hf,
   .ref = 90,
   .trunc = 120
)

# MACE excluding patients with history of stroke or MI
wo_MACE_wo_uacr <- plot_imputed_outcomes(
   .predictor = predictors,
   outcome_element = "MACE_wo_HF",
   outcome_name = "MACE",
   .covariates = covariates[!covariates %in% c("history_stroke", "history_mi")],
   .imp_data = data_imp_wo_mace,
   .ref = 90,
   .trunc = 120
)

# AKI excluding patients with history of AKI
wo_AKI_wo_uacr <- plot_imputed_outcomes(
   .predictor = predictors,
   outcome_element = "aki",
   outcome_name = "AKI",
   .covariates = covariates,
   .imp_data = data_imp_wo_aki,
   .ref = 90,
   .trunc = 120
)

wo_converted_UACR <- append(
   death_kfrt_wo_uacr,
   list("HF" = wo_hf_wo_uacr, "MACE" = wo_MACE_wo_uacr, "AKI" = wo_AKI_wo_uacr)
)

# wo_UACR_plots <- append(x = wo_converted_UACR, list("HF" = wo_converted_UACR_HF))
save(
   wo_converted_UACR,
   file = here::here("output", "r_objects", "MICE_wo-UACR-conversion.rda")
)


## Tables ----------------------------------------------------------------------

# Tables with HRs
wo_converted_UACR_tbl <- c(`All-cause Death` = "death", "KFRT" = "rrt") %>%
   map2(
      .x = .,
      .y = names(.),
      ~ summarize_HR_MICE(
         .predictor = predictors,
         outcome_element = .x,
         outcome_name = {{ .y }},
         .covariates = covariates,
         .imp_data = wo_uacr,
         .ref = 90,
         .trunc = 120
      )
   )

# Heart failure excluding patients with history of HF
wo_hf_wo_uacr_tbl <- summarize_HR_MICE(
   .predictor = predictors,
   outcome_element = "heart_failure",
   outcome_name = "Heart Failure",
   .covariates = covariates[!covariates %in% c("hf")],
   .imp_data = data_imp_wo_hf,
   .ref = 90,
   .trunc = 120
)

# MACE excluding patients with history of stroke or MI
wo_MACE_wo_uacr_tbl <- summarize_HR_MICE(
   .predictor = predictors,
   outcome_element = "MACE_wo_HF",
   outcome_name = "MACE",
   .covariates = covariates[!covariates %in% c("history_stroke", "history_mi")],
   .imp_data = data_imp_wo_mace,
   .ref = 90,
   .trunc = 120
)

# AKI excluding patients with history of AKI
wo_AKI_wo_uacr_tbl <- summarize_HR_MICE(
   .predictor = predictors,
   outcome_element = "aki",
   outcome_name = "AKI",
   .covariates = covariates,
   .imp_data = data_imp_wo_aki,
   .ref = 90,
   .trunc = 120
)

wo_UACR_listoftbl <- append(
   wo_converted_UACR_tbl,
   list(
      "Heart Failure" = wo_hf_wo_uacr_tbl,
      "MACE" = wo_MACE_wo_uacr_tbl,
      "AKI" = wo_AKI_wo_uacr_tbl
   )
)

wo_UACR_tbl <- purrr::reduce(.x = wo_UACR_listoftbl, .f = function(x, y) {
   left_join(x, y, by = c("x1", "predictor"))
}) %>%
   relocate(
      predictor,
      GFR = x1,
      `All-cause Death`,
      KFRT,
      AKI,
      `Heart Failure`,
      MACE
   )

save(
   wo_UACR_tbl,
   file = here::here("output", "r_objects", "MICE_wo-UACR-conversion_tbl.rda")
)
