################################################################################
## mGFR and outcomes in SCREAM ##
# Code for missing data description and imputation
# Written by: Antoine Creon
# Date: 2024-11-13
# Update: 2025-03-13 (imputation without converting PCR and disptick into UACR first)
################################################################################

################################################################################
# LOAD DATA AND PACKAGES #######################################################
################################################################################

source(here::here("code", "02_analysis-preparation.R"))

# Load data with events NOT CENSORED AT KFRT
data <- read_rds(
   file = here::here("data", "cleaned", "data_not_cens_kfrt_named.rds")
)


################################################################################
################################################################################
###                                                                          ###
###                                 UACR ENRICHED                            ###
###                                                                          ###
################################################################################
################################################################################

################################################################################
# MISSING DATA IMPUTATION WITH MICE ############################################
################################################################################
# Chosen based on Austin PC, Computational Statistics 2024
# (flavors of FCS including SMC-FCS perform similarly in realistic settings)

## Helper functions ------------------------------------------------------------
#  compute cubic splines
compute_cubic_splines <- function(.data, .variable) {
   rcspline.eval(
      x = .data[[.variable]],
      knots = quantile(.data[[.variable]], probs = c(0.05, 0.35, 0.65, 0.95)),
      nk = 4,
      norm = 2, # /!\ by default norm = 2, would then yield a different result than above
      pc = FALSE,
      inclx = TRUE
   ) %>%
      as.data.frame |>
      rename_with(~ paste0(.variable, "_s", seq_along(.)), everything())
}

# Create a Nelson Aalen column from the event name (with time_to_event and H0_event)
create_nelsonaalen_column <- function(.data, event) {
   # Create the name of the new column
   new_col_name <- paste0("H0_", event)
   time_col <- paste0("time_to_", event)

   # Use `mutate()` and `rlang` functions to create the new column dynamically
   nelson_aalen <- .data %>%
      mutate(
         !!new_col_name := nelsonaalen(
            .,
            time = !!sym(time_col),
            status = !!sym(event)
         )
      ) %>%
      select(!!new_col_name)

   return(nelson_aalen)
}

## Prepare the dataset with cubic splines --------------------------------------

# Create the cubic spline transformations of the predictors
predictors_with_cubic_splines <- predictors |>
   map(~ compute_cubic_splines(data, .x)) |>
   list_cbind()

# B. Create all cumulative hazard functions
nelson_aalen_estimates <- unname(events)[
   !events %in% c("mi", "stroke", "cv_death", "af", "pad")
] |>
   map(~ create_nelsonaalen_column(data, .x)) |>
   list_cbind()

# C. cbind splines and cumulative hazards with the rest of the data
data_with_cubic_splines <- cbind(
   data |>
      select(
         lopnr,
         all_of(covariates),
         "kidney_donor",
         "history_aki",
         all_of(unlist(unname(events))),
         -all_of(unlist(predictors)), # avoid redundancy between pred and pred_s1
         starts_with("time_to"),
         "log_uacr_enriched"
      ),
   predictors_with_cubic_splines,
   nelson_aalen_estimates
)

## Perform the imputation with cubic splines -----------------------------------

# A. Blank MICE to get predictor matrix and method list
imp0_cs <- mice(data_with_cubic_splines, m = 5, maxit = 0)

# B. Imputation method
imp0_cs$method["log_uacr_enriched"] # pmm
imp0_cs$method["bmi"] # pmm

# C. Predictor Matrix: only impute log_UACR and BMI
imp0_cs$predictorMatrix[] <- 0

# use everything except lopnr (ID column) and time to events
imp0_cs$predictorMatrix["log_uacr_enriched", ] <- 1
imp0_cs$predictorMatrix["bmi", ] <- 1
imp0_cs$predictorMatrix[, "lopnr"] <- 0
imp0_cs$predictorMatrix[, grepl(
   "^time_to",
   colnames(imp0_cs$predictorMatrix)
)] <- 0

# D. COMPUTE MICE IN PARALLEL
data_imp_cs <- futuremice(
   data_with_cubic_splines,
   method = imp0_cs$method, # imputation method
   predictorMatrix = imp0_cs$predictorMatrix, # predictor matrix
   m = 50, # nb imputations
   maxit = 10, # nb iterations : min 10
   parallelseed = 3091992, # seed over the parallel backend
   n.core = parallelly::availableCores(logical = TRUE) - 1 # nb of cores used
)

write_rds(
   data_imp_cs,
   file = here::here("data", "cleaned", "not_cens_kfrt_imp_cs.rds")
)


################################################################################
################################################################################
###                                                                          ###
###                               UACR NOT ENRICHED                          ###
###                                                                          ###
################################################################################
################################################################################

# Load the dataset
data <- qread(file = here::here("data", "cleaned", "UACR-not-enriched.qs"))

# Create the cubic spline transformations of the predictors
predictors_with_cubic_splines <- predictors |>
   map(~ compute_cubic_splines(data, .x)) |>
   list_cbind()

# B. Create all cumulative hazard functions
nelson_aalen_estimates <- unname(events)[
   !events %in% c("mi", "stroke", "cv_death", "af", "pad")
] |>
   map(~ create_nelsonaalen_column(data, .x)) |>
   list_cbind()

# C. cbind splines and cumulative hazards with the rest of the data
data_with_cubic_splines <- cbind(
   data |>
      select(
         all_of(covariates),
         "kidney_donor",
         "history_aki",
         all_of(unlist(unname(events))),
         -all_of(unlist(predictors)), # avoid redundancy between pred and pred_s1
         starts_with("time_to"),
         "log_uacr_enriched"
      ),
   predictors_with_cubic_splines,
   nelson_aalen_estimates
)

# Perform imputation
imp0_cs <- mice(data_with_cubic_splines, m = 5, maxit = 0)

imp0_cs$method["log_uacr_enriched"] # pmm
imp0_cs$method["bmi"] # pmm

imp0_cs$predictorMatrix[] <- 0
imp0_cs$predictorMatrix["log_uacr_enriched", ] <- 1 # use everything to impute log(UACR)
imp0_cs$predictorMatrix["bmi", ] <- 1 # use everything to impute BMI
imp0_cs$predictorMatrix[, "lopnr"] <- 0
imp0_cs$predictorMatrix[, grepl(
   "^time_to",
   colnames(imp0_cs$predictorMatrix)
)] <- 0

data_imp_cs <- futuremice(
   data_with_cubic_splines,
   method = imp0_cs$method, # imputation method
   predictorMatrix = imp0_cs$predictorMatrix, # predictor matrix
   m = 50, # nb imputations
   maxit = 10, # nb iterations : min 10
   parallelseed = 3091992, # seed over the parallel backend
   n.core = parallelly::availableCores(logical = TRUE) - 1 # nb of cores used
)

write_rds(
   data_imp_cs,
   file = here::here("data", "cleaned", "MICE_UACR-not-enriched.rds")
)
