################################################################################
# Project: mGFR and outcomes in SCREAM
# Purpose: eGFRcr-Mortality association in the source pop versus study pop
# Written by: Antoine Creon
# Date: 2025-10-19
################################################################################

################################################################################
# LOAD DATA AND PACKAGES #######################################################
################################################################################

# Load packages, helper functions and covariate/predictor lists
source(here::here("code", "02_analysis-preparation.R"))

# Datasets
## All individuals with one creatinine available
creat <- qread(file = here::here("data", "cleaned", "sample_pop_wo_mgfr.qs"))

## Study population
mgfr <- qread(file = here::here("data", "cleaned", "data_not_cens_kfrt.qs"))

################################################################################
# CREATININE VERSUS MORTALITY IN THE CREAT/CYS/mGFR POP  #######################
################################################################################
# Plot with eGFRcr versus mortality in:
# - Full population with creatinine available
# - Study population, with creatinine, cystatin C and mGFR available

## Hazard ratios for mortality -------------------------------------------------

# Pop with creat
fullpop_creat_HRs <- compute_fits(
   "death",
   "ckd_epi_2021_cr",
   c(
      "uacr_enriched_missing_indicator",
      covariates[!covariates %in% c("bmi", "log_uacr_enriched")]
   ),
   creat
) %>%
   extract_termplot(ref = 90, 120) |>
   mutate(predictor = "Creatinine")

# Keep one logHR per GFR value
fullpop_creat_logHRs <- fullpop_creat_HRs |>
   mutate(x1 = round(x1)) |>
   group_by(predictor, x1) |>
   dplyr::summarize(yhat = mean(log(y1))) |>
   ungroup()

# pop with creat + cys + mGFR
mgfrpop_creat_HRs <- compute_fits(
   "death",
   "ckd_epi_2021_cr",
   c(
      "uacr_enriched_missing_indicator",
      covariates[!covariates %in% c("bmi", "log_uacr_enriched")]
   ),
   mgfr
) %>%
   extract_termplot(ref = 90, 120) |>
   mutate(predictor = "mGFR")

creatinine_mortality <- rbind(
   fullpop_creat_HRs,
   # fullpop_creatcys_HRs,
   mgfrpop_creat_HRs
)

# save the data for the table
save(
   creatinine_mortality,
   file = here::here("data", "cleaned", "creatinine-pop_mortality.RData")
)

## CONDITIONAL INCIDENCE RATES - full pop --------------------------------------

# Prepare function call
covariates_creat_pop <- c(
   "uacr_enriched_missing_indicator",
   covariates[!covariates %in% c("bmi", "log_uacr_enriched")]
)
.covariates_string <- paste(covariates_creat_pop, collapse = " + ")
.poisson_formula <- paste0(
   "death",
   " == 1 ~ rcs(",
   "ckd_epi_2021_cr",
   ",4) + ",
   .covariates_string,
   " + offset(log(",
   "time_to_death",
   "))"
)

# Fit Poisson model
pois_fit_creat_pop <- rms::Glm(
   formula = as.formula(.poisson_formula),
   family = poisson(link = "log"),
   data = creat
)

# Compute median for each covariate in creat (mode for factors)
get_typical_value <- function(x) {
   x_non_na <- x[!is.na(x)]

   if (is.numeric(x) || is.integer(x)) {
      median(x_non_na, na.rm = TRUE)
   } else if (is.factor(x)) {
      mode_level <- names(which.max(table(x_non_na)))
      factor(mode_level, levels = levels(x))
   } else if (is.logical(x)) {
      as.logical(round(median(as.integer(x_non_na), na.rm = TRUE)))
   } else {
      x_non_na[[which.max(tabulate(match(x_non_na, unique(x_non_na))))]]
   }
}

median_values <- creat |>
   dplyr::summarise(
      dplyr::across(all_of(covariates_creat_pop), get_typical_value),
      .groups = "drop"
   ) |>
   as.list()

# Step 2: Create prediction data at median values for target GFR categories
new_data <- as_tibble_row(median_values) |>
   select(all_of(covariates_creat_pop)) |>
   uncount(length(c(15, 30, 45, 60, 75, 90, 120))) |>
   mutate(
      ckd_epi_2021_cr = c(15, 30, 45, 60, 75, 90, 120),
      time_to_death = 1, # Is ignored when using rms::Glm() (but not stats::glm())
      death = 0
   )

# Step 3: For each conditions, predict log-rate and its variance
# Predict on link scale (log)
# /!\ Here, only X*beta is returned, not X*beta + offset. The rate will be per day
pred <- predict(
   pois_fit_creat_pop,
   newdata = new_data,
   type = "lp",
   se.fit = TRUE
)

log_rates <- pred$linear.predictors + log(365.25) # Add the offset manually when using rms::Glm() to convert to person-year
se_log_rates <- pred$se.fit

conditional_IR <- exp(log_rates)
CI_lower <- exp(log_rates - qnorm(0.975) * se_log_rates)
CI_upper <- exp(log_rates + qnorm(0.975) * se_log_rates)

# Results
IR_fullpop <- tibble(
   outcome = "death",
   gfr_type = "ckd_epi_2021_cr",
   gfr_thresholds = c(15, 30, 45, 60, 75, 90, 120),
   conditional_IR = conditional_IR,
   CI_lower = CI_lower,
   CI_upper = CI_upper
)

## CONDITIONAL INCIDENCE RATES - study pop -------------------------------------

# Prepare function call
covariates_creat_pop <- c(
   "uacr_enriched_missing_indicator",
   covariates[!covariates %in% c("bmi", "log_uacr_enriched")]
)
.covariates_string <- paste(covariates_creat_pop, collapse = " + ")
.poisson_formula <- paste0(
   "death",
   " == 1 ~ rcs(",
   "ckd_epi_2021_cr",
   ",4) + ",
   .covariates_string,
   " + offset(log(",
   "time_to_death",
   "))"
)

# Fit Poisson model
pois_fit_mgfr_pop <- rms::Glm(
   formula = as.formula(.poisson_formula),
   family = poisson(link = "log"),
   data = mgfr
)

# Compute median for each covariate in mgfr (mode for factors)
get_typical_value <- function(x) {
   x_non_na <- x[!is.na(x)]

   if (is.numeric(x) || is.integer(x)) {
      median(x_non_na, na.rm = TRUE)
   } else if (is.factor(x)) {
      mode_level <- names(which.max(table(x_non_na)))
      factor(mode_level, levels = levels(x))
   } else if (is.logical(x)) {
      as.logical(round(median(as.integer(x_non_na), na.rm = TRUE)))
   } else {
      x_non_na[[which.max(tabulate(match(x_non_na, unique(x_non_na))))]]
   }
}

median_values <- mgfr |>
   dplyr::summarise(
      dplyr::across(all_of(covariates_creat_pop), get_typical_value),
      .groups = "drop"
   ) |>
   as.list()

# Step 2: Create prediction data at median values for target GFR categories
new_data <- as_tibble_row(median_values) |>
   select(all_of(covariates_creat_pop)) |>
   uncount(length(c(15, 30, 45, 60, 75, 90, 120))) |>
   mutate(
      ckd_epi_2021_cr = c(15, 30, 45, 60, 75, 90, 120),
      time_to_death = 1, # Is ignored when using rms::Glm() (but not stats::glm())
      death = 0
   )

# Step 3: For each conditions, predict log-rate and its variance
# Predict on link scale (log)
# /!\ Here, only X*beta is returned, not X*beta + offset. The rate will be per day
pred <- predict(
   pois_fit_mgfr_pop,
   newdata = new_data,
   type = "lp",
   se.fit = TRUE
)

log_rates <- pred$linear.predictors + log(365.25) # Add the offset manually when using rms::Glm() to convert to person-year
se_log_rates <- pred$se.fit

conditional_IR <- exp(log_rates)
CI_lower <- exp(log_rates - qnorm(0.975) * se_log_rates)
CI_upper <- exp(log_rates + qnorm(0.975) * se_log_rates)

# Results
IR_creatpop <- tibble(
   outcome = "death",
   gfr_type = "ckd_epi_2021_cr",
   gfr_thresholds = c(15, 30, 45, 60, 75, 90, 120),
   conditional_IR = conditional_IR,
   CI_lower = CI_lower,
   CI_upper = CI_upper
)

## RATIO OF HAZARD RATIOS ------------------------------------------------------
# We estimated the HRs in two ~independent samples
# So HR_fullpop / HR_mgfr = log(HR_fullpop_creat) - log(HR_mgfr)
# SE[log(HR_fullpop_creat) - log(HR_mgfr)] = sqrt[SE(log(HR_fullpop_creat))^2 + SE(log(HR_mgfr))^2]
# CI: Delta +/- 1.96 * SE[Delta]

# Pick the x1 value closest to the integer
logHRs_creat_mgfr <- creatinine_mortality |>
   mutate(logy1 = log(y1)) |>
   group_by(predictor, x1_round = round(x1)) |>
   slice_min(abs(x1 - x1_round), n = 1, with_ties = FALSE) |>
   ungroup() |>
   mutate(x1 = round(x1)) |>
   select(x1, logy1, se1, predictor)

# Compute the contrast and their SE
# logHRs_contrast <- logHRs_creat_mgfr |>
RHR_creat_mgfr <- logHRs_creat_mgfr |>
   pivot_wider(names_from = predictor, values_from = c(logy1, se1)) |>
   rowwise() |>
   mutate(
      logHR_contrast = logy1_Creatinine - logy1_mGFR,
      se_contrast = sqrt(se1_Creatinine^2 + se1_mGFR^2),
      CI_lower = logHR_contrast - qnorm(0.975) * se_contrast,
      CI_upper = logHR_contrast + qnorm(0.975) * se_contrast
   ) |>
   mutate(
      RHR = exp(logHR_contrast),
      CI_lower_exp = exp(CI_lower),
      CI_upper_exp = exp(CI_upper)
   ) |>
   select(x1, RHR, CI_lower_exp, CI_upper_exp)

save(
   RHR_creat_mgfr,
   file = here::here("data", "cleaned", "RHR_creat_mgfr_mortality.RData")
)
