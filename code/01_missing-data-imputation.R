################################################################################
## mGFR and outcomes in SCREAM ##
# Code for missing data description and imputation
# Written by: Antoine Creon
# Date: 2024-11-13
# Update: 2025-03-13 (imputation without converting PCR and disptick into UACR first)
################################################################################

################################################################################
###                         LOAD DATA AND PACKAGES                           ###
################################################################################

pacman::p_load(
   tidyverse,
   glue,
   here,
   labelled,
   mice,
   micemd,
   mitools,
   ggmice,
   smcfcs,
   rms,
   splines,
   gt,
   gtsummary,
   qs
)

# Load data with events NOT CENSORED AT KFRT
data <- read_rds(
   file = here::here("data", "cleaned", "data_not_cens_kfrt_named.rds")
)

source(here::here("code", "05_helper-functions-survival.R"))
source(here::here(
   "code",
   "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"
))


################################################################################
################################################################################
###                                                                          ###
###                                 UACR ENRICHED                            ###
###                                                                          ###
################################################################################
################################################################################

################################################################################
#                        MISSING DATA IMPUTATION WITH MICE                     #
################################################################################
# Chosen based on Austin PC, Computational Statistics 2024
# (flavors of FCS including SMC-FCS perform similarly in realistic settings)

# ---------------------------- Helper functions --------------------------------

# add cubic spline transformation of variables, 4 knots

# # add linear spline transformation of variables
# compute_linear_splines <- function(.data, .variable){
#
#    splines::bs(.data[[.variable]], knots = c(30, 60, 90), degree = 1, Boundary.knots = range(.data[[.variable]])) %>%
#       as.data.frame |>
#       rename_with(~ paste0(.variable, "_ls", seq_along(.)), everything())
# }

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

# Create interation between Nelson-Aalen estimate and covariates
create_interaction <- function(.data, .nelson, .covariates) {
   .data |>
      select(all_of(.covariates), all_of(.nelson)) |>
      mutate(across(everything(), ~ .x * .data[[.nelson]])) |>
      rename_with(.fn = ~ paste0(.x, "_int_", .nelson))
}


# ------------------Prepare the dataset with cubic splines ---------------------

# Create the cubic spline transformations of the predictors
predictors_with_cubic_splines <- predictors |>
   map(~ compute_cubic_splines(data, .x)) |>
   list_cbind()

# B. Create all cumulative hazard functions
neslson_aalen_estimates <- unname(events)[
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
   neslson_aalen_estimates
)

# D. create interactions between H0(t) and covariates (including GFRs and their spline transformation)
## Compile the names of all covariates/predictors
all_predictors <- c(covariates, names(predictors_with_cubic_splines))

## create an interaction with all covariates for each of the H(t) and cbind with the rest
# data_to_imp_cs <- names(neslson_aalen_estimates) |>
#    map(~ create_interaction(data_with_cubic_splines, .x, all_predictors)) |>
#    list_cbind() |>
#    cbind(data_with_cubic_splines)
#
# names(data_to_imp_cs)
# TOO MANY VARIABLES (248): interactions between H0(t) and all covariates not included

# ------------------ Perform the imputation with cubic splines -----------------

# A. Blank MICE to get predictor matrix and method list
imp0_cs <- mice(data_with_cubic_splines, m = 5, maxit = 0)

# B. Imputation method
imp0_cs$method["log_uacr_enriched"] # pmm
imp0_cs$method["bmi"] # pmm

# C. Predictor Matrix: only impute log_UACR and BMI
imp0_cs$predictorMatrix[] <- 0

# use everything except lopnr (ID)
imp0_cs$predictorMatrix["log_uacr_enriched", ] <- 1
imp0_cs$predictorMatrix["bmi", ] <- 1
imp0_cs$predictorMatrix[, "lopnr"] <- 0


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
neslson_aalen_estimates <- unname(events)[
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
   neslson_aalen_estimates
)

# Perform imputation
imp0_cs <- mice(data_with_cubic_splines, m = 5, maxit = 0)

imp0_cs$method["log_uacr_enriched"] # pmm
imp0_cs$method["bmi"] # pmm

imp0_cs$predictorMatrix[] <- 0
imp0_cs$predictorMatrix["log_uacr_enriched", ] <- 1 # use everything to impute log(UACR)
imp0_cs$predictorMatrix["bmi", ] <- 1 # use everything to impute BMI

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
