################################################################################
# Project: mGFR and outcomes in SCREAM
# Purpose: Bootstrap-MI (BOOT-MI) contrasts of log hazard ratios
# Author: Antoine Creon
# Date: 2026-04-19
# Note: This script performs 500 bootstrap replicates with nested MICE imputation
#       Order: Bootstrap sampling -> MICE -> Analysis -> Pool estimates
################################################################################

# Part of this script runs on a high-performance computing cluster and
# depends on external variables storing exact paths
# (project_dir_HPC, project_dir_P, library_path_HPC)

################################################################################
# PREPARE THE ANALYSIS #########################################################
################################################################################

## Paths -----------------------------------------------------------------------
# On the computer cluster
setwd(project_dir_HPC)

## Results directory
results_dir_HPC <- paste(
  "inc-events",
  "MICE",
  "logHR",
  sep = "_"
)

results_path_HPC <- file.path(project_dir_HPC, "output", results_dir_HPC)

if (!dir.exists(results_path_HPC)) {
  dir.create(results_path_HPC, recursive = TRUE)
}

# On the P: drive (local)
results_path_P <- file.path(
  project_dir_P,
  "output",
  results_dir_HPC
)

if (!dir.exists(results_path_P)) {
  dir.create(results_path_P, recursive = TRUE)
}

# Library path for cluster workers
.libPaths(library_path_HPC)

## Packages --------------------------------------------------------------------
packages <- c(
  "dplyr",
  "readr",
  "tidyr",
  "purrr",
  "here",
  "conflicted",
  "survival",
  "rms",
  "Hmisc",
  "mice",
  "mitools",
  "boot"
)

for (pkg in packages) {
  library(pkg, character.only = TRUE)
}

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

## Data ------------------------------------------------------------------------
# Read the raw dataset (pre-imputation)
data_raw <- read_rds(
  file = file.path(
    project_dir_P,
    "data",
    "cleaned",
    "data_not_cens_kfrt_named.rds"
  )
)

# Read covariate names
covariates <- read_rds(
  file = file.path(project_dir_P, "data", "cleaned", "covariates_names.rds")
)

# List the predictors
predictors <- list(
  "mgfr",
  "ckd_epi_2021_cr",
  "ckd_epi_2012_cys",
  "ckd_epi_2021_cr_cys"
)

# List the events
events <- list(
  "All-cause Death" = "death",
  "KFRT" = "rrt",
  "AKI" = "aki",
  "MACE" = "MACE_wo_HF",
  "Heart Failure" = "heart_failure"
)

## Seed for reproducibility
set.seed(19920903, "L'Ecuyer-CMRG")


################################################################################
# HELPER FUNCTIONS #############################################################
################################################################################

## Prepare data with cubic splines (for MICE imputation) -----------------------
prepare_data_for_imputation <- function(.data) {
  # Helper: compute cubic splines
  compute_cubic_splines <- function(.data, .variable) {
    rcspline.eval(
      x = .data[[.variable]],
      knots = quantile(.data[[.variable]], probs = c(0.05, 0.35, 0.65, 0.95)),
      nk = 4,
      norm = 2,
      pc = FALSE,
      inclx = TRUE
    ) |>
      as.data.frame() |>
      rename_with(~ paste0(.variable, "_s", seq_along(.)), everything())
  }

  # Helper: compute Nelson-Aalen
  create_nelsonaalen_column <- function(.data, event) {
    new_col_name <- paste0("H0_", event)
    time_col <- paste0("time_to_", event)

    nelson_aalen <- .data %>%
      mutate(
        !!new_col_name := nelsonaalen(
          .,
          time = !!sym(time_col),
          status = !!sym(event)
        )
      ) |>
      select(!!new_col_name)

    return(nelson_aalen)
  }

  # Create cubic spline transformations
  predictors_with_cubic_splines <- predictors |>
    map(~ compute_cubic_splines(.data, .x)) |>
    list_cbind()

  # Create Nelson-Aalen estimates
  na_estimates <- unname(events) |>
    map(~ create_nelsonaalen_column(.data, .x)) |>
    list_cbind()

  # Combine all components
  data_prepared <- cbind(
    .data |>
      select(
        lopnr,
        all_of(covariates),
        "kidney_donor",
        "history_aki",
        all_of(unlist(unname(events))),
        -all_of(unlist(predictors)),
        starts_with("time_to"),
        "log_uacr_enriched"
      ),
    predictors_with_cubic_splines,
    na_estimates
  )

  return(data_prepared)
}

## Perform MICE imputation -----------------------------------------------------
perform_mice <- function(.data, .m = 50, .maxit = 10, .seed = NULL) {
  if (!is.null(.seed)) {
    set.seed(.seed)
  }

  # Blank mice to set up predictor matrix and methods
  imp0 <- mice(.data, m = 5, maxit = 0)

  # Set predictor matrix: only impute log_uacr_enriched and bmi, using everything but
  # time to event and lopnr (ID) as predictors
  imp0$predictorMatrix[] <- 0
  imp0$predictorMatrix["log_uacr_enriched", ] <- 1
  imp0$predictorMatrix["bmi", ] <- 1
  imp0$predictorMatrix[, "lopnr"] <- 0
  imp0$predictorMatrix[, grepl("^time_to", colnames(imp0$predictorMatrix))] <- 0

  # Perform imputation
  data_imp <- mice(
    .data,
    method = imp0$method,
    predictorMatrix = imp0$predictorMatrix,
    m = .m,
    maxit = .maxit,
    printFlag = FALSE
  )

  return(data_imp)
}

## Fit survival models on imputed datasets -------------------------------------
compute_fits_imp_data <- function(
  .outcome,
  .predictor,
  .covariates,
  .imp_data,
  .spline = "cubic",
  .linear_knots = "c(30, 60, 90)",
  .weights = NULL
) {
  .imp_predictor <- paste0(.predictor, "_s1")
  .time <- paste0("time_to_", .outcome)
  .covariates_str <- paste(.covariates, collapse = " + ")

  if (.spline == "linear") {
    .surv <- paste0(
      "Surv(",
      .time,
      ",",
      .outcome,
      ") ~ splines::bs(",
      .imp_predictor,
      ", knots =",
      .linear_knots,
      ", degree = 1) +"
    )
  } else if (.spline == "cubic") {
    .surv <- paste0(
      "Surv(",
      .time,
      ",",
      .outcome,
      ") ~ rcs(",
      .imp_predictor,
      ",4) +"
    )
  } else {
    rlang::abort("spline must be either 'linear' or 'cubic'")
  }

  .fits <- with(
    .imp_data,
    coxph(as.formula(paste0(.surv, .covariates_str)), weights = .weights)
  )

  return(.fits$analyses)
}

## Extract and pool termplots (BOOT-MI simplified, no Rubin pooling) -----------
extract_each_termplotof_imp <- function(fit, ref, truncate) {
  ptemp <- termplot(fit, se = FALSE, plot = FALSE)
  x1 <- ptemp[[1]][["x"]]
  y1 <- ptemp[[1]][["y"]]
  yref1 <- -y1[abs(x1 - ref) == min(abs(x1 - ref))]

  dat1 <- as.data.frame(cbind(x1, y1)) |>
    filter(x1 <= truncate) |>
    mutate(y1 = y1 + yref1) |>
    filter(x1 >= 15) |>
    mutate(predictor = names(ptemp)[[1]]) |>
    mutate(
      predictor = factor(
        predictor,
        levels = c(
          "mgfr_s1",
          "ckd_epi_2021_cr_s1",
          "ckd_epi_2012_cys_s1",
          "ckd_epi_2021_cr_cys_s1"
        ),
        labels = c("mGFR", "eGFR[cr]", "eGFR[cys]", "eGFR[cr-cys]")
      )
    )

  return(dat1)
}

extract_termplot_imp <- function(.imp_fits, ref, truncate) {
  .imp_fits |>
    map(~ extract_each_termplotof_imp(.x, ref, truncate))
}

pool_termplots <- function(.extracted_termplots) {
  .extracted_termplots |>
    bind_rows() |>
    # Pick the x1 value closest to the integer
    group_by(predictor, x1_round = round(x1)) |>
    slice_min(abs(x1 - x1_round), n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(x1 = round(x1), yhat = y1)
}


## Compute coefficient for all outcomes ----------------------------------------
compute_coeff <- function(
  .predictor,
  outcome_element,
  outcome_name,
  .covariates,
  .imp_data,
  .ref,
  .trunc,
  .spline = "cubic",
  .linear_knots = "c(30,60,90)",
  .weights = NULL
) {
  .predictor |>
    map(
      ~ compute_fits_imp_data(
        .outcome = outcome_element,
        .predictor = .x,
        .covariates = .covariates,
        .imp_data = .imp_data,
        .spline = .spline,
        .linear_knots = .linear_knots,
        .weights = .weights
      )
    ) |>
    map(~ extract_termplot_imp(.x, ref = .ref, truncate = .trunc)) |>
    map(~ pool_termplots(.x)) |>
    list_rbind() |>
    mutate(.outcome = outcome_name)
}


################################################################################
# BOOT-MI FUNCTION ############################################################
################################################################################

boot_mi_analysis <- function(data, indices) {
  # Draw bootstrap sample
  data_boot <- data[indices, ]

  # Prepare data for imputation (add splines and Nelson-Aalen)
  data_prep <- prepare_data_for_imputation(data_boot)

  # Perform MICE imputation
  data_imp <- perform_mice(data_prep, .m = 50, .maxit = 10, .seed = NULL)

  # Create subset mids objects for analyses excluding specific outcomes
  data_imp_wo_mace <- data_imp |>
    complete(action = "long", include = TRUE) |>
    filter(history_mi == 0, history_stroke == 0) |>
    mice::as.mids()

  data_imp_wo_hf <- data_imp |>
    complete(action = "long", include = TRUE) |>
    filter(hf == 0) |>
    mice::as.mids()

  data_imp_wo_aki <- data_imp |>
    complete(action = "long", include = TRUE) |>
    filter(history_aki == 0) |>
    mice::as.mids()

  # Run analysis for each outcome
  death <- compute_coeff(
    .predictor = predictors,
    outcome_element = "death",
    outcome_name = "death",
    .covariates = covariates,
    .imp_data = data_imp,
    .ref = 90,
    .trunc = 120
  )

  KFRT <- compute_coeff(
    .predictor = predictors,
    outcome_element = "rrt",
    outcome_name = "KFRT",
    .covariates = covariates,
    .imp_data = data_imp,
    .ref = 90,
    .trunc = 120
  )

  inc_hf <- compute_coeff(
    .predictor = predictors,
    outcome_element = "heart_failure",
    outcome_name = "HF",
    .covariates = covariates[!covariates %in% c("hf")],
    .imp_data = data_imp_wo_hf,
    .ref = 90,
    .trunc = 120
  )

  inc_MACE <- compute_coeff(
    .predictor = predictors,
    outcome_element = "MACE_wo_HF",
    outcome_name = "MACE",
    .covariates = covariates[
      !covariates %in% c("history_stroke", "history_mi")
    ],
    .imp_data = data_imp_wo_mace,
    .ref = 90,
    .trunc = 120
  )

  inc_AKI <- compute_coeff(
    .predictor = predictors,
    outcome_element = "aki",
    outcome_name = "AKI",
    .covariates = covariates,
    .imp_data = data_imp_wo_aki,
    .ref = 90,
    .trunc = 120
  )

  # Format output as named vector
  boot_output <- bind_rows(death, KFRT, inc_hf, inc_MACE, inc_AKI) |>
    mutate(
      predictor = recode(
        predictor,
        "mGFR" = "mGFR",
        "eGFR[cr]" = "creat",
        "eGFR[cys]" = "cys",
        "eGFR[cr-cys]" = "cr.cys"
      )
    ) |>
    pivot_wider(
      names_from = c(.outcome, predictor, x1),
      values_from = yhat
    ) |>
    unlist()

  return(boot_output)
}


################################################################################
# RUN BOOTSTRAP ################################################################
################################################################################

# Run boot with parallel processing
# This will perform B=500 bootstrap replicates on the HPC cluster
boot_results <- boot(
  data = data_raw,
  statistic = boot_mi_analysis,
  R = 500,
  parallel = "multicore",
  ncpus = parallelly::availableCores(logical = TRUE) - 1,
  sim = "ordinary"
)

# Save to HPC
boot_results_path_HPC <- file.path(results_path_HPC, "boot_results.rds")
write_rds(boot_results, file = boot_results_path_HPC)


################################################################################
# COPY RESULTS TO P: DRIVE #####################################################
################################################################################

# List result fits on the HPC drive
files <- list.files(
  path = results_path_HPC,
  pattern = paste0("^boot_"),
  full.names = TRUE
)

# Copy
file.copy(files, results_path_P, overwrite = TRUE)

message("Bootstrap-MICE analysis completed!")
message(paste("Results saved to:", results_path_HPC))
message(paste("Results copied to:", results_path_P))


################################################################################
# ANALYSIS ON THE LOCAL DRIVE ##################################################
################################################################################

boot <- read_rds(here::here(
  "output",
  "inc-events_MICE_logHR",
  "boot_results.rds"
))

## Extract bootstrap results ---------------------------------------------------

# Create column names for the bootstrapped results
boot_colnames <- crossing(
  boot_id = 1:500,
  outcome = c("death", "KFRT", "inc_hf", "inc_MACE", "inc_AKI"),
  predictor = c("mGFR", "creat", "cys", "cr.cys"),
  gfr = 15:120
) |>
  mutate(
    predictor = factor(predictor, levels = c("mGFR", "creat", "cys", "cr.cys")),
    outcome = factor(
      outcome,
      levels = c("death", "KFRT", "inc_hf", "inc_MACE", "inc_AKI")
    )
  ) |>
  arrange(boot_id, outcome, predictor, gfr)

# Add the actual results and compute the HR = exp(logHR)
res <- boot$t

res_tbl <- as_tibble(res) |>
  pivot_longer(
    cols = everything(),
    values_to = "logHR"
  ) |>
  select(-name) %>%
  cbind(boot_colnames, .) |>
  mutate(HR = exp(logHR))


## HR and logHR ----------------------------------------------------------------

indiv_estimates <- res_tbl |>
  pivot_longer(
    cols = c(logHR, HR),
    names_to = "parameter",
    values_to = "value"
  ) |>
  group_by(outcome, predictor, gfr, parameter) |>
  summarize(
    estimate = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  )

## HR ratios (= logHR differences) ---------------------------------------------

HR_ratios <- res_tbl |>
  select(-logHR) |>
  pivot_wider(
    names_from = predictor,
    values_from = HR
  ) |>
  mutate(
    ratio_creat = creat / mGFR,
    ratio_cys = cys / mGFR,
    ratio_cr.cys = cr.cys / mGFR
  ) |>
  pivot_longer(
    cols = starts_with("ratio_"),
    names_to = "contrast",
    values_to = "HR_ratio"
  ) |>
  group_by(outcome, gfr, contrast) |>
  summarize(
    estimate = median(HR_ratio),
    lower = quantile(HR_ratio, 0.025),
    upper = quantile(HR_ratio, 0.975)
  ) |>
  mutate(
    contrast = recode(
      contrast,
      "ratio_creat" = "creat / mGFR",
      "ratio_cys" = "cys / mGFR",
      "ratio_cr.cys" = "cr.cys / mGFR"
    )
  )
