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

pacman::p_load(tidyverse, glue, 
               here,
               labelled,
               mice, micemd, mitools, ggmice, smcfcs,
               rms, splines,
               gt, gtsummary,
               qs)

# Load data with events NOT CENSORED AT KFRT
data <- read_rds(file = here::here("data","cleaned","data_not_cens_kfrt_named.rds")) 

source(here::here("code", "05_helper-functions-survival.R"))
source(here::here("code", "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"))



################################################################################
################################################################################
###                                                                          ###
###                                 UACR ENRICHED                            ###
###                                                                          ###
################################################################################
################################################################################

################################################################################
#                            MISSING DATA DESCRIPTION                          #
################################################################################

# ------------- NUMBER AND PROP. OF MISSING DATA FOR UACR ----------------------

# Add complete cases variable and label it
data <- data %>%
   mutate(complete_case = ifelse(complete.cases(data %>% select(all_of(covariates))) == TRUE, "Yes", "No"))

data <- data %>% 
   set_variable_labels(complete_case = "Complete cases")

data %>% dplyr::summarize("Prop UACR Missing (%)" = round(mean(is.na(log_uacr_enriched)) *100, 2),
                          "Count UACR Missing" = sum(is.na(log_uacr_enriched)),
                          "Prop BMI Missing (%)" = round(mean(is.na(bmi)) *100, 2),
                          "Count BMI Missing" = sum(is.na(bmi)),
                          "Prop of missing data (%)" = round(mean(as.numeric(complete_case == "No"))*100, 2)) |>
   gt()


# ---------------- COMPARISON OF COMPLETE AND INCOMPLETE CASES -----------------
# Imbalance between complete and incomplete cases, not compatible with MACR mechanism

uacr <- c("log_uacr_enriched", "uacr_enriched", "uacr_enriched_missing_indicator")

complete_cases_table1 <- data |>
   select(-all_of(uacr), -calendar_year, -starts_with("time_to_"), - all_of(unlist(events)), - history_pad,
          - ihd, - cevd, -arrh, - antihypertensives) |>
   tbl_summary(by = complete_case,
               statistic = list(all_continuous() ~ "{mean} ({sd})",
                                all_categorical() ~ "{n} ({p}%)"),
               missing = "no") %>%
   add_overall() %>%
   add_difference(everything() ~ "smd") %>%
   modify_spanning_header(c("stat_1", "stat_2") ~ "**Complete case**") %>%
   modify_header(label ~ "**Baseline Characteristics**",
                 estimate ~ "**SMD**") |> 
   modify_column_hide(conf.low)

save(complete_cases_table1, file = here::here("output","r_objects","CCA_description_tbl.rda"))


#  -------------------- PATTERN IN MISSING OBSERVATIONS ?-----------------------
# No pattern to observe as only one variable is missing

plot_pattern(data |> select(-starts_with("time_to"), -all_of(unlist(events))),
              vrb = "all",
              square = TRUE,
              rotate = TRUE,
              npat = NULL,
              caption = TRUE)


################################################################################
#                        MISSING DATA IMPUTATION WITH MICE                     #
################################################################################
# Chosen based on Austin PC, Computational Statistics 2024 
# (flavors of FCS including SMC-FCS perform similarly in realistic settings)
# Content of the model: 
# all predictors
# + all cumulative hazard functions H(t) to approximate H0(t) (Nelson Aalen estimator)
# + all event indicator variables
# + interactions between each H(t) and the predictors

# ---------------------------- Helper functions --------------------------------

#  compute cubic splines
compute_cubic_splines <- function(.data, .variable){
   
   rcspline.eval(x = .data[[.variable]], 
                 knots = quantile(.data[[.variable]],probs=c(0.05,0.35,0.65,0.95)),
                 nk = 4,
                 norm = 2, # /!\ by default norm = 2, would then yield a different result than above
                 pc = FALSE,
                 inclx = TRUE) %>%
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
      mutate(!!new_col_name := nelsonaalen(., time = !!sym(time_col), status = !!sym(event))) %>%
      select(!!new_col_name)
   
   return(nelson_aalen)
}

# Create interation between Nelson-Aalen estimate and covariates
create_interaction <- function(.data, .nelson, .covariates) {
   .data |>
      select(all_of(.covariates), all_of(.nelson)) |>
      mutate(across(everything(), ~ .x * .data[[.nelson]])) |>
      rename_with(.fn = ~paste0(.x, "_int_", .nelson))
}


# ------------------Prepare the dataset with cubic splines ---------------------

# Create the cubic spline transformations of the predictors
predictors_with_cubic_splines <- predictors |>
   map(~ compute_cubic_splines(data, .x)) |>
   list_cbind()

# B. Create all cumulative hazard functions
neslson_aalen_estimates <- unname(events)[!events %in% c("mi","stroke", "cv_death", "af", "pad")] |>
   map(~ create_nelsonaalen_column(data, .x)) |>
   list_cbind()

# C. cbind splines and cumulative hazards with the rest of the data
data_with_cubic_splines <- cbind(data |> select(all_of(covariates),
                                                "kidney_donor",
                                                "history_aki",
                                                all_of(unlist(unname(events))), 
                                                - all_of(unlist(predictors)), # avoid redundancy between pred and pred_s1
                                                starts_with("time_to"),
                                                "log_uacr_enriched"), 
                                 predictors_with_cubic_splines,
                                 neslson_aalen_estimates) 

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

# C. Predictor Matrix: only impute log_UACR
imp0_cs$predictorMatrix[] <- 0
imp0_cs$predictorMatrix["log_uacr_enriched", ] <- 1 # use everything to impute log(UACR)
imp0_cs$predictorMatrix["bmi", ] <- 1 # use everything to impute BMI

# D. COMPUTE MICE IN PARALLEL
data_imp_cs <- futuremice(data_with_cubic_splines,
                          method = imp0_cs$method, # imputation method
                          predictorMatrix = imp0_cs$predictorMatrix, # predictor matrix
                          m = 50, # nb imputations
                          maxit = 10, # nb iterations : min 10
                          parallelseed = 3091992, # seed over the parallel backend
                          n.core = parallelly::availableCores(logical = TRUE) - 1 # nb of cores used
)


# ----------------------------- Save the result --------------------------------

write_rds(data_imp_cs, file = here::here("data","cleaned","not_cens_kfrt_imp_cs.rds"))


################################################################################
#                      MISSING DATA IMPUTATION WITH SMC-FCS                    #
################################################################################
## SUBSTANTIVE MODEL-COMPATIBLE FULLY CONDITIONAL SPECIFICATION
## considered the gold standard for imputation of survival analysis with competing risks
## /!\ in smcfcs :
## - variables must be in the same order in data.frame, imp_method and formula
## - it MUST be a data.frame and NOT a tibble !!
## - the event must be numeric
## ONE SET OF IMPUTATIONS FOR EACH OUTCOME COMPETING WITH DEATH

# ---------------------------- Helper functions --------------------------------

# A. to write all the cox models
write_models <- function(outcome, predictors, covariates){
   
   # time to event variable
   outcome_compet <- paste0(outcome, "_compet")
   futime_compet <- paste0("time_to_", outcome_compet)
   
   # string of covariates
   .covariates <- paste0(covariates, collapse = " + ")
   
   # string of predictor and its spline transformations
   .predictors_spline <- predictors |>
      map(~ glue("{.x}_s1 + {.x}_s2 + {.x}_s3")) |>
      paste0(collapse = " + ")
   
   # Surv function and ~
   .formula1 <- glue("Surv({futime_compet}, {outcome_compet} == 1) ~ {.predictors_spline} + {.covariates}")
   .formula2 <- glue("Surv({futime_compet}, {outcome_compet} == 2) ~ {.predictors_spline} + {.covariates}")
   
   .formula <- list(.formula1, .formula2)
   return(.formula)
}


# B. to select the data for each set of imputations
create_smcfcs_dataset <- function(.data, .outcome, .covariates){
   
   # define the new columns
   time_to_outcome <- paste0("time_to_", .outcome)
   outcome_compet <- paste0(.outcome, "_compet")
   futime_compet <- paste0("time_to_", outcome_compet)
   outcome_and_futime <- c(futime_compet, outcome_compet)
   
   .data_organized <- .data |>
      mutate({{futime_compet}} := pmin(time_to_death, get({{time_to_outcome}})),
             {{outcome_compet}} := case_when(get({{.outcome}}) == 1 & get({{time_to_outcome}}) == get({{futime_compet}}) ~ 1,
                                             death == 1 & time_to_death == get({{futime_compet}}) ~ 2,
                                             TRUE ~ 0)) |>
      
      
      select(-all_of(unlist(unname(events))), 
             -starts_with("time_to"),
             all_of(c(outcome_and_futime, .covariates)))  |>
      relocate(all_of(.covariates), .after = last_col()) |>
      relocate(outcome_and_futime) |>
      as.data.frame() # must be a dataframe for smc fcs to run
   
   return(.data_organized)
}


# ------------------------ Write the substantive models ------------------------

# each model include all GFRs and their spline transformation
models_rrt <- write_models("rrt", predictors, covariates)
models_aki <- write_models("aki", predictors, covariates)
models_hf <- write_models("heart_failure", predictors, covariates)
models_MACEwo <- write_models("MACE_wo_HF", predictors, covariates)
models_MACEwi <- write_models("MACE_with_HF", predictors, covariates)
model_death <- "Surv(time_to_death, death) ~ mgfr_s1 + mgfr_s2 + mgfr_s3 + ckd_epi_2021_cr_s1 + ckd_epi_2021_cr_s2 + ckd_epi_2021_cr_s3 + ckd_epi_2012_cys_s1 + ckd_epi_2012_cys_s2 + ckd_epi_2012_cys_s3 + ckd_epi_2021_cr_cys_s1 + ckd_epi_2021_cr_cys_s2 + ckd_epi_2021_cr_cys_s3 + age + female + bmi + history_mi + hyperten + hf + history_stroke + history_AF + pvd + dm + copd + cancer + liver + transplant + antihypertensives + lipid + nsaid + log_uacr_enriched + glucocorticosteroids"

# ---------------------------- Prepare datasets --------------------------------

# Remove H0(t) from dataset
smcfcs_data <- data_with_cubic_splines |>
   select(-starts_with("H0_"))

# One dataset per outcome competing with death (no loop to modify for each variable the set of covariates)
smcfcs_data_rrt <- create_smcfcs_dataset(smcfcs_data, "rrt", covariates)
smcfcs_data_aki <- create_smcfcs_dataset(smcfcs_data, "aki", covariates)
smcfcs_data_hf <- create_smcfcs_dataset(smcfcs_data, "heart_failure", covariates)
smcfcs_data_MACEwo <- create_smcfcs_dataset(smcfcs_data, "MACE_wo_HF", covariates)
smcfcs_data_MACEwi <- create_smcfcs_dataset(smcfcs_data, "MACE_with_HF", covariates)
smcfcs_data_death <- smcfcs_data |>
   select(time_to_death, death, all_of(covariates), ends_with("_s1"), ends_with("_s2"), ends_with("_s3")) |>
   data.frame()

# ----------------- SET UP IMPUTATION  METHOD FOR UACR AND BMI------------------

# Blank mice to get predictor matrix and method list
imp_method_smcfcs <- list(smcfcs_data_rrt, smcfcs_data_aki, smcfcs_data_hf, smcfcs_data_MACEwo, smcfcs_data_MACEwi, smcfcs_data_death) |>
   map(~ mice(., m = 5, maxit = 0)) |>
   map(~ .x$method) |> 
   map(~ if_else(.x == "pmm", "norm", .x)) # linear regression for BMI and UACR

#  PREDICTOR MATRIX
# default option : all variables except outcomes (which is  is implicitly 
# conditioned on by the rejection sampling scheme used by smcfcs, and should not
# be included.)

# --------------------------------- RUN SMC-FCS --------------------------------

seed <- 3091992

# A. KFRT
set.seed(seed)
smcfcs_rrt <- smcfcs(smcfcs_data_rrt, 
                     smtype = "compet",
                     smformula = unlist(models_rrt),
                     method = imp_method_smcfcs[[1]],
                     m = 50,
                     numit = 10,
                     rjlimit = 10000)

smcfcs_rrt <- mitools::imputationList(smcfcs_rrt$impDatasets)

write_rds(smcfcs_rrt, file = here::here("data","cleaned","smcfcs_rrt.rds"))

# B. AKI
set.seed(seed)
smcfcs_aki <- smcfcs(smcfcs_data_aki, 
                     smtype = "compet",
                     smformula = unlist(models_aki),
                     method = imp_method_smcfcs[[2]],
                     m = 50,
                     numit = 10,
                     rjlimit = 10000)

smcfcs_aki <- mitools::imputationList(smcfcs_aki$impDatasets)

write_rds(smcfcs_aki, file = here::here("data","cleaned","smcfcs_aki.rds"))

# C. HF
set.seed(seed)
smcfcs_hf <- smcfcs(smcfcs_data_hf, 
                    smtype = "compet",
                    smformula = unlist(models_hf),
                    method = imp_method_smcfcs[[3]],
                    m = 50,
                    numit = 10,
                    rjlimit = 10000)

smcfcs_hf <- mitools::imputationList(smcfcs_hf$impDatasets)

write_rds(smcfcs_hf, file = here::here("data","cleaned","smcfcs_hf.rds"))

# D. 3-point MACE
set.seed(seed)
smcfcs_MACEwo <- smcfcs(smcfcs_data_MACEwo, 
                        smtype = "compet",
                        smformula = unlist(models_MACEwo),
                        method = imp_method_smcfcs[[4]],
                        m = 50,
                        numit = 10,
                        rjlimit = 10000)

smcfcs_MACEwo <- mitools::imputationList(smcfcs_MACEwo$impDatasets)

write_rds(smcfcs_MACEwo, file = here::here("data","cleaned","smcfcs_MACEwo.rds"))

# E. 4-point MACE
set.seed(seed)
smcfcs_MACEwi <- smcfcs(smcfcs_data_MACEwi, 
                        smtype = "compet",
                        smformula = unlist(models_MACEwi),
                        method = imp_method_smcfcs[[5]],
                        m = 50,
                        numit = 10,
                        rjlimit = 10000)

smcfcs_MACEwi <- mitools::imputationList(smcfcs_MACEwi$impDatasets)

write_rds(smcfcs_MACEwi, file = here::here("data","cleaned","smcfcs_MACEwi.rds"))

# F. Death
set.seed(seed)
smcfcs_death <- smcfcs(smcfcs_data_death, 
                       smtype = "coxph",
                       smformula = model_death,
                       method = imp_method_smcfcs[[6]],
                       m = 50,
                       numit = 10,
                       rjlimit = 10000)

smcfcs_death <- mitools::imputationList(smcfcs_death$impDatasets)

write_rds(smcfcs_death, file = here::here("data","cleaned","smcfcs_death.rds"))



################################################################################
################################################################################
###                                                                          ###
###                               UACR NOT ENRICHED                          ###
###                                                                          ###
################################################################################
################################################################################

# Load the dataset
data <- qread(file = here::here("data","cleaned","UACR-not-enriched.qs"))

# Describe missing data
data <- data %>%
   mutate(complete_case = ifelse(complete.cases(data %>% select(all_of(covariates))) == TRUE, "Yes", "No"))

data <- data %>% 
   set_variable_labels(complete_case = "Complete cases")

data %>% dplyr::summarize("Prop UACR Missing (%)" = round(mean(is.na(log_uacr_enriched)) *100, 2),
                          "Count UACR Missing" = sum(is.na(log_uacr_enriched)),
                          "Prop BMI Missing (%)" = round(mean(is.na(bmi)) *100, 2),
                          "Count BMI Missing" = sum(is.na(bmi)),
                          "Prop of missing data (%)" = round(mean(as.numeric(complete_case == "No"))*100, 2)) |>
   gt()

uacr <- c("log_uacr_enriched", "uacr_enriched", "uacr_enriched_missing_indicator")

data |>
   select(-all_of(uacr)) |>
   tbl_summary(by = complete_case,
               statistic = list(all_continuous() ~ "{mean} ({sd})",
                                all_categorical() ~ "{n} ({p}%)"),
               missing = "no") %>%
   add_overall() %>%
   add_difference(everything() ~ "smd") %>%
   modify_spanning_header(c("stat_1", "stat_2") ~ "**Complete case**") %>%
   modify_header(label ~ "**Baseline Characteristics**",
                 estimate ~ "**SMD**")


# Create the cubic spline transformations of the predictors
predictors_with_cubic_splines <- predictors |>
   map(~ compute_cubic_splines(data, .x)) |>
   list_cbind()

# B. Create all cumulative hazard functions
neslson_aalen_estimates <- unname(events)[!events %in% c("mi","stroke", "cv_death", "af", "pad")] |>
   map(~ create_nelsonaalen_column(data, .x)) |>
   list_cbind()

# C. cbind splines and cumulative hazards with the rest of the data
data_with_cubic_splines <- cbind(data |> select(all_of(covariates),
                                                "kidney_donor",
                                                "history_aki",
                                                all_of(unlist(unname(events))), 
                                                - all_of(unlist(predictors)), # avoid redundancy between pred and pred_s1
                                                starts_with("time_to"),
                                                "log_uacr_enriched"), 
                                 predictors_with_cubic_splines,
                                 neslson_aalen_estimates) 


# Perform imputation
imp0_cs <- mice(data_with_cubic_splines, m = 5, maxit = 0)

imp0_cs$method["log_uacr_enriched"] # pmm 
imp0_cs$method["bmi"] # pmm 

imp0_cs$predictorMatrix[] <- 0
imp0_cs$predictorMatrix["log_uacr_enriched", ] <- 1 # use everything to impute log(UACR)
imp0_cs$predictorMatrix["bmi", ] <- 1 # use everything to impute BMI

data_imp_cs <- futuremice(data_with_cubic_splines,
                          method = imp0_cs$method, # imputation method
                          predictorMatrix = imp0_cs$predictorMatrix, # predictor matrix
                          m = 50, # nb imputations
                          maxit = 10, # nb iterations : min 10
                          parallelseed = 3091992, # seed over the parallel backend
                          n.core = parallelly::availableCores(logical = TRUE) - 1 # nb of cores used
)

write_rds(data_imp_cs, file = here::here("data","cleaned","MICE_UACR-not-enriched.rds"))


