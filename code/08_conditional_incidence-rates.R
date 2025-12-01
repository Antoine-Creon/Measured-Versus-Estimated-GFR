# ##############################################################################
#  Study: mGFR versus outcomes
#
#  Conditional incidence rates by eGFR category
#  Update: 2025-11-28
#  Author: Antoine Creon
#  Change logs: none
#
# ##############################################################################

#  ------------------- Load imputed dataset and functions ----------------------

source(here::here("code", "05_helper-functions-survival.R"))
source(here::here(
  "code",
  "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"
))

data_imp <- read_rds(
  file = here::here("data", "cleaned", "not_cens_kfrt_imp_cs.rds")
)

# Add GFR categories to the main dataset
gfr_cuts <- c(0, 15, 30, 45, 60, 75, 90, 105, 999) # mGFR alvays <150 but not eGFR
gfr_labels <- c(
  "<15",
  "15-30",
  "30-45",
  "45-60",
  "60-75",
  "75-90",
  "90-105",
  ">105"
)

data_imp <- data_imp |>
  complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
  mutate(
    across(
      c(
        mgfr_s1,
        ckd_epi_2021_cr_s1,
        ckd_epi_2012_cys_s1,
        ckd_epi_2021_cr_cys_s1
      ),
      \(x) cut(x, breaks = gfr_cuts, labels = gfr_labels, right = FALSE),
      .names = "{.col}_cat"
    )
  ) |>
  rename_with(
    \(x) str_remove(x, "_s1_cat$"),
    ends_with("_s1_cat")
  ) |>
  rename_with(
    \(x) paste0(x, "_cat"),
    matches("^(mgfr|ckd_epi)")
  ) |>
  mice::as.mids()

## Filter incident patients for each outcome -----------------------------------
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
# HELPER FUNCTION  #############################################################
# ##############################################################################

cond_IR <- function(
  imputed_dataset,
  .outcome,
  .predictor,
  .covariates
) {
  ## Step 0: Prepare variables -------------------------------------------------

  # Time to oucome
  .time_to_outcome <- paste0("time_to_", .outcome)

  # name of the predictor as in imputed data
  .imp_predictor <- paste0(.predictor, "_cat")

  # string of covariates
  .covariates_string <- paste(.covariates, collapse = " + ")

  # GFR categories
  gfr_labels <- imputed_dataset[[1]] |> pull(.imp_predictor) |> levels()

  ## Step 1: Model fitting in MI datasets and pool coefficients -----------------
  .poisson_formula <- paste0(
    .outcome,
    " == 1 ~ ",
    .imp_predictor,
    " + ",
    .covariates_string
  )

  fit_mi <- with(
    imputed_dataset,
    glm(
      formula = as.formula(.poisson_formula),
      family = poisson(link = "log"),
      offset = log(get(.time_to_outcome))
    )
  )

  # Stack all imputed datasets to compute overall medians
  all_imputed <- map_dfr(
    seq_len(imputed_dataset$m),
    ~ complete(imputed_dataset, .x)
  )

  # Compute median for each covariate (mode for factors)
  median_values <- list()

  for (cov in .covariates) {
    col <- all_imputed[[cov]]
    if (is.numeric(col)) {
      median_values[[cov]] <- median(col, na.rm = TRUE)
    } else if (is.factor(col)) {
      # Use mode (most frequent level) for factors
      median_values[[cov]] <- names(which.max(table(col)))
    } else {
      median_values[[cov]] <- col[1]
    }
  }

  # Step 2: Create prediction data at median values for target GFR categories
  new_data <- as_tibble_row(median_values) |>
    select(all_of(.covariates)) |>
    uncount(length(gfr_labels)) |>
    mutate(
      !!.imp_predictor := factor(gfr_labels, levels = gfr_labels),
      !!.time_to_outcome := 365.25,
      !!.outcome := 0
    )

  # Step 3: For each imputation, predict log-rate and its variance
  m <- length(fit_mi$analyses)
  n_pred <- nrow(new_data)

  # Store predictions and variances
  log_rates <- matrix(NA, nrow = n_pred, ncol = m)
  var_log_rates <- matrix(NA, nrow = n_pred, ncol = m)

  for (i in seq_len(m)) {
    fit_i <- fit_mi$analyses[[i]]

    # Predict on link scale (log)
    pred <- predict(fit_i, newdata = new_data, type = "link", se.fit = TRUE)

    log_rates[, i] <- pred$fit
    var_log_rates[, i] <- pred$se.fit^2
  }

  # Step 4: Pool using Rubin's rules (for each prediction row)
  # Q_bar = mean of point estimates
  Q_bar <- rowMeans(log_rates)

  # W = mean within-imputation variance
  W <- rowMeans(var_log_rates)

  # B = between-imputation variance
  B <- apply(log_rates, 1, var)

  # Total variance
  T_var <- W + (1 + 1 / m) * B
  T_se <- sqrt(T_var)

  # Step 5: Compute CIs and exponentiate
  log_rate_pooled <- Q_bar

  conditional_IR <- exp(log_rate_pooled)
  CI_lower <- exp(log_rate_pooled - qnorm(0.975) * T_se)
  CI_upper <- exp(log_rate_pooled + qnorm(0.975) * T_se)

  # Results
  tibble(
    outcome = .outcome,
    gfr_type = .predictor,
    gfr_cat = gfr_labels,
    conditional_IR = conditional_IR,
    CI_lower = CI_lower,
    CI_upper = CI_upper
  )
}


# ##############################################################################
# COMPUTE ESTIMATES  ###########################################################
# ##############################################################################

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
      outcome == "aki" ~ list(covariates[!covariates %in% c("history_aki")]),
      outcome == "MACE_wo_HF" ~ list(covariates[
        !covariates %in% c("history_stroke", "history_mi")
      ])
    )
  )

# Apply std_IR to each row using pmap
IR_res <- conditions |>
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
  IR_res,
  file = here::here("output", "r_objects", "conditional_incidence_rates.rda")
)
