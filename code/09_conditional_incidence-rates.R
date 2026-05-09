################################################################################
# Project: mGFR and outcomes in SCREAM
# Purpose: Conditional incidence rates, with continuous GFR
# Written by: Antoine Creon
# Date: 2026-04-29
################################################################################

################################################################################
# LOAD DATA AND PACKAGES #######################################################
################################################################################

source(here::here("code", "02_analysis-preparation.R"))

## Load data -------------------------------------------------------------------

# Data not imputed
data <- read_rds(
  file = here::here("data", "cleaned", "data_not_cens_kfrt_named.rds")
)

# Data imputed
data_imp <- read_rds(
  file = here::here("data", "cleaned", "not_cens_kfrt_imp_cs.rds")
)

# Add EKFC and CKDEPI2009-2012 eGFR to original dataset
data <- data |>
  mutate(
    ekfc_cr = ekfc_cr(creat, age, female),
    ekfc_cys = ekfc_cys(cystatin, age),
    ckd_epi_2009_cr = ckd_epi_2009_cr(creat, age, female),
    ckd_epi_2012_cr_cys = ckd_epi_2012_cr_cys(creat, cystatin, age, female)
  ) |>
  rowwise() |>
  mutate(ekfc_combined = mean(c(ekfc_cr, ekfc_cys))) |>
  ungroup() |>
  relocate(ekfc_cr, ekfc_combined, ekfc_cys)

# Add EKFC eGFR with cubic splines to the imputed dataset
data_imp <- data_imp |>
  complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
  mutate(
    ekfc_cr_s1 = rep(data$ekfc_cr, 51),
    ekfc_cys_s1 = rep(data$ekfc_cys, 51),
    ekfc_combined_s1 = rep(data$ekfc_combined, 51),
    ckd_epi_2009_cr_s1 = rep(data$ckd_epi_2009_cr, 51),
    ckd_epi_2012_cr_cys_s1 = rep(data$ckd_epi_2012_cr_cys, 51)
  ) |>
  mice::as.mids() # reconstruct into a mids object if needed

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


# ##############################################################################
# MAIN ANALYSIS  ###############################################################
# ##############################################################################

## Define conditions for each analysis -----------------------------------------

conditions <- crossing(
  outcome = c("death", "rrt", "aki", "MACE_wo_HF", "heart_failure"),
  predictor = c(
    "mgfr",
    "ckd_epi_2021_cr",
    "ckd_epi_2012_cys",
    "ckd_epi_2021_cr_cys"
  )
)

conditions <- conditions |>
  mutate(
    dataset = case_when(
      outcome %in% c("death", "rrt") ~ list(data_imp),
      outcome == "heart_failure" ~ list(data_imp_wo_hf),
      outcome == "aki" ~ list(data_imp_wo_aki),
      outcome == "MACE_wo_HF" ~ list(data_imp_wo_mace)
    ),
    covariates = case_when(
      outcome %in% c("death", "rrt") ~ list(covariates),
      outcome == "heart_failure" ~ list(covariates[
        -which(covariates %in% c("hf"))
      ]),
      outcome == "aki" ~ list(covariates),
      outcome == "MACE_wo_HF" ~ list(covariates[
        !covariates %in% c("history_stroke", "history_mi")
      ])
    )
  )

## Compute conditional incidence rates for each combination --------------------

IR_main <- conditions |>
  pmap(
    \(outcome, predictor, dataset, covariates) {
      cond_IR(
        imputed_dataset = dataset,
        .outcome = outcome,
        .predictor = predictor,
        .covariates = covariates
      )
    },
    .progress = TRUE
  ) |>
  bind_rows()

save(
  IR_main,
  file = here::here(
    "output",
    "r_objects",
    "conditional_incidence_rates_continuous_GFR.rda"
  )
)


# ##############################################################################
# CKD-EPI 2009-2012  ###########################################################
# ##############################################################################

## Define conditions for each analysis -----------------------------------------

conditions_ckdepi2009 <- conditions |>
  mutate(
    predictor = recode(
      predictor,
      "ckd_epi_2021_cr" = "ckd_epi_2009_cr",
      "ckd_epi_2021_cr_cys" = "ckd_epi_2012_cr_cys"
    )
  )


## Compute conditional incidence rates for each combination --------------------

IR_ckdepi2009 <- conditions_ckdepi2009 |>
  pmap(
    \(outcome, predictor, dataset, covariates) {
      cond_IR(
        imputed_dataset = dataset,
        .outcome = outcome,
        .predictor = predictor,
        .covariates = covariates
      )
    },
    .progress = TRUE
  ) |>
  bind_rows()

save(
  IR_ckdepi2009,
  file = here::here(
    "output",
    "r_objects",
    "conditional_incidence_rates_continuous_GFR_CKDEPI2009.rda"
  )
)


# ##############################################################################
# EKFC  ########################################################################
# ##############################################################################

## Define conditions for each analysis -----------------------------------------

conditions_ekfc <- conditions |>
  mutate(
    predictor = recode(
      predictor,
      "ckd_epi_2021_cr" = "ekfc_cr",
      "ckd_epi_2012_cys" = "ekfc_cys",
      "ckd_epi_2021_cr_cys" = "ekfc_combined"
    )
  )

## Compute conditional incidence rates for each combination --------------------

IR_ekfc <- conditions_ekfc |>
  pmap(
    \(outcome, predictor, dataset, covariates) {
      cond_IR(
        imputed_dataset = dataset,
        .outcome = outcome,
        .predictor = predictor,
        .covariates = covariates
      )
    },
    .progress = TRUE
  ) |>
  bind_rows()

save(
  IR_ekfc,
  file = here::here(
    "output",
    "r_objects",
    "conditional_incidence_rates_continuous_GFR_EKFC.rda"
  )
)


# ##############################################################################
# WITHOUT UACR CONVERSION  #####################################################
# ##############################################################################

## Load data -------------------------------------------------------------------

wo_uacr <- read_rds(
  file = here::here("data", "cleaned", "MICE_UACR-not-enriched.rds")
)

# EXCLUDE patients with history of MACE (MI and stroke)
imp_wo_uacr_wo_mace <- wo_uacr |>
  complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
  filter(history_mi == 0, history_stroke == 0) |> # Apply filter to each dataset
  mice::as.mids() # reconstruct into a mids object if needed

# EXCLUDE patients with history of heart failure
imp_wo_uacr_wo_hf <- wo_uacr |>
  complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
  filter(hf == 0) |> # Apply filter to each dataset
  mice::as.mids() # reconstruct into a mids object if needed

# EXCLUDE patients with history of AKI
imp_wo_uacr_wo_aki <- wo_uacr |>
  complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
  filter(history_aki == 0) |> # Apply filter to each dataset
  mice::as.mids() # reconstruct into a mids object if needed

## Define conditions for each analysis -----------------------------------------

conditions_wo_uacr <- conditions |>
  mutate(
    dataset = case_when(
      outcome %in% c("death", "rrt") ~ list(wo_uacr),
      outcome == "heart_failure" ~ list(imp_wo_uacr_wo_hf),
      outcome == "aki" ~ list(imp_wo_uacr_wo_aki),
      outcome == "MACE_wo_HF" ~ list(imp_wo_uacr_wo_mace)
    )
  )

## Compute conditional incidence rates for each combination --------------------

IR_wo_uacr <- conditions_wo_uacr |>
  pmap(
    \(outcome, predictor, dataset, covariates) {
      cond_IR(
        imputed_dataset = dataset,
        .outcome = outcome,
        .predictor = predictor,
        .covariates = covariates
      )
    },
    .progress = TRUE
  ) |>
  bind_rows()

save(
  IR_wo_uacr,
  file = here::here(
    "output",
    "r_objects",
    "conditional_incidence_rates_continuous_GFR_WO_UACR.rda"
  )
)


# ##############################################################################
# WITHOUT KTR  #################################################################
# ##############################################################################

## Filter patients  ------------------------------------------------------------

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

## Define conditions for each analysis -----------------------------------------

conditions_wo_ktr <- conditions |>
  mutate(
    dataset = case_when(
      outcome %in% c("death", "rrt") ~ list(data_imp_wo_ktr),
      outcome == "heart_failure" ~ list(data_imp_wo_hf_ktr),
      outcome == "aki" ~ list(data_imp_wo_aki_ktr),
      outcome == "MACE_wo_HF" ~ list(data_imp_wo_mace_ktr)
    ),
    covariates = case_when(
      outcome %in% c("death", "rrt") ~ list(.env$covariates[
        !.env$covariates %in% c("transplant")
      ]),
      outcome == "heart_failure" ~ list(.env$covariates[
        !.env$covariates %in% c("transplant", "hf")
      ]),
      outcome == "aki" ~ list(.env$covariates[
        !.env$covariates %in% c("transplant")
      ]),
      outcome == "MACE_wo_HF" ~ list(.env$covariates[
        !.env$covariates %in% c("transplant", "history_stroke", "history_mi")
      ])
    )
  )

## Compute conditional incidence rates for each combination --------------------

IR_wo_ktr <- conditions_wo_ktr |>
  pmap(
    \(outcome, predictor, dataset, covariates) {
      cond_IR(
        imputed_dataset = dataset,
        .outcome = outcome,
        .predictor = predictor,
        .covariates = covariates
      )
    },
    .progress = TRUE
  ) |>
  bind_rows()

save(
  IR_wo_ktr,
  file = here::here(
    "output",
    "r_objects",
    "conditional_incidence_rates_continuous_GFR_WO_KTR.rda"
  )
)

# ##############################################################################
# NO RESTRICTION TO INCIDENT EVENTS  ###########################################
# ##############################################################################

## Define conditions for each analysis -----------------------------------------

conditions <- crossing(
  outcome = c("aki", "MACE_wo_HF", "heart_failure"),
  predictor = c(
    "mgfr",
    "ckd_epi_2021_cr",
    "ckd_epi_2012_cys",
    "ckd_epi_2021_cr_cys"
  )
)

conditions_prev_comorb <- conditions |>
  mutate(
    dataset = case_when(
      outcome == "heart_failure" ~ list(data_imp),
      outcome == "aki" ~ list(data_imp),
      outcome == "MACE_wo_HF" ~ list(data_imp)
    ),
    covariates = case_when(
      outcome %in% c("heart_failure", "aki", "MACE_wo_HF") ~ list(covariates)
    )
  )

## Compute conditional incidence rates for each combination --------------------

IR_prevalent_comorb <- conditions_prev_comorb |>
  pmap(
    \(outcome, predictor, dataset, covariates) {
      cond_IR(
        imputed_dataset = dataset,
        .outcome = outcome,
        .predictor = predictor,
        .covariates = covariates
      )
    },
    .progress = TRUE
  ) |>
  bind_rows()

save(
  IR_prevalent_comorb,
  file = here::here(
    "output",
    "r_objects",
    "conditional_incidence_rates_continuous_GFR_prev_comorb.rda"
  )
)
