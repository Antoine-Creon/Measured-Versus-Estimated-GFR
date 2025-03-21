################################################################################
# Project: mGFR and outcomes in SCREAM 
# Purpose: Code for analysis GFR versus outcomes, after imputation by  SMC-FCS
# Written by: Antoine Creon
# Date: 2024-11-13
################################################################################

#  ------------------- Load imputed dataset and functions ----------------------

source(here::here("code", "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"))

smcfcs_rrt <- read_rds(file = here::here("data","cleaned","smcfcs_rrt.rds"))
smcfcs_aki <- read_rds(file = here::here("data","cleaned","smcfcs_aki.rds"))
smcfcs_hf <- read_rds(file = here::here("data","cleaned","smcfcs_hf.rds"))
smcfcs_MACEwo <- read_rds(file = here::here("data","cleaned","smcfcs_MACEwo.rds"))
smcfcs_death <- read_rds(file = here::here("data","cleaned","smcfcs_death.rds"))


smcfcs_datasets <- list("rrt" = smcfcs_rrt, 
                        "aki" = smcfcs_aki, 
                        "heart_failure" = smcfcs_hf, 
                        "MACE_wo_HF" = smcfcs_MACEwo,
                        "death" = smcfcs_death)

# filter patients without HF at baseline
## helper function
smcfcs_hf <- smcfcs_hf[["imputations"]] |>
      map(~ filter(.x, hf == 0)) |> 
      imputationList()

#  ---------------------------- CUBIC SPLINE PLOTS -----------------------------


# KFRT
kfrt_smcfcs_cs_plot <- plot_imputed_outcomes(.predictor = predictors, 
                                             outcome_element = "rrt", outcome_name = "KFRT", 
                                             .covariates =  covariates, 
                                             .imp_data = smcfcs_rrt, 
                                             .ref = 90, .trunc = 120) +  
   scale_y_log10()  

# AKI
aki_smcfcs_cs_plot <- plot_imputed_outcomes(.predictor = predictors, 
                                            outcome_element = "aki", outcome_name = "AKI", 
                                            .covariates =  covariates, 
                                            .imp_data = smcfcs_aki, 
                                            .ref = 90, .trunc = 120)

# MACE without HF
MACE_wo_HF_smcfcs_cs_plot <- plot_imputed_outcomes(.predictor = predictors, 
                                                   outcome_element = "MACE_wo_HF", outcome_name = "MACE", 
                                                   .covariates =  covariates, 
                                                   .imp_data = smcfcs_MACEwo, 
                                                   .ref = 90, .trunc = 120)

# Heart failure
HF_smcfcs_cs_plot <- plot_imputed_outcomes(.predictor = predictors, 
                                           outcome_element = "heart_failure", outcome_name = "Heart Failure", 
                                           .covariates =  covariates[!covariates %in% c("hf")], 
                                           .imp_data = smcfcs_hf, 
                                           .ref = 90, .trunc = 120)


# Death 
## Modify the function to fit the specificity of the data
compute_fits_imp_data_death <- function(.outcome, .predictor, .covariates, .imp_data, .spline = "cubic", .linear_knots = "c(30, 60, 90)", .weights = NULL) {
   
   # name of the predictor as in imputed data
   .imp_predictor <- paste0(.predictor, "_s1")
   
   # Name the variables of the model depending on the imputation method
   .time <- paste0("time_to_", .outcome)   # time to event variable
   .range_predictor <- paste0(range(.imp_data[["imputations"]][[1]][[.imp_predictor]]), collapse = ",")
   
   # string of covariates
   .covariates <- paste(.covariates, collapse = " + ")
   
   # Surv function and ~
   if (.spline == "linear"){
      .surv <- paste0("Surv(", .time, ",", .outcome, ") ~ splines::bs(", .imp_predictor,", knots =", .linear_knots, ", degree = 1) +")
      
   } else if (.spline == "cubic"){
      .surv <- paste0("Surv(", .time, ",", .outcome, ") ~ rcs(", .imp_predictor,",4) +")
      
   } else {
      rlang::abort("spline must be either 'linear' or 'cubic'")
   }
   
   # Compute the fit in each imputed dataset
   .fits <- with(.imp_data, coxph(as.formula(paste0(.surv, .covariates)), weights = .weights))
   
   # keep only each analysis
   if (class(.imp_data) == "mids") {
      .imp_fits <- .fits$analyses
      return(.imp_fits)
   } else if (class(.imp_data) == "imputationList"){
      return(.fits)
   } else {
      rlang::abort(".imp_data must be either of class 'mids' (for MICE) or 'imputationList' (for SMC-FCS)")
   }
   
}

## Plot death
death_smcfcs_cs_plot <- predictors %>%
   map( ~ compute_fits_imp_data_death(.outcome = "death", 
                                .predictor = .x, 
                                .covariates = covariates, 
                                .imp_data = smcfcs_death)) |>
   map(~ extract_termplot_imp(.x, 90, 120)) |>
   map(~ pool_termplots(.x)) |>
   list_rbind() |>
   plot_HRs(.outcome = "All-cause death")

smcfcs_withKTR_cs_plots <- list(kfrt_smcfcs_cs_plot, aki_smcfcs_cs_plot, death_smcfcs_cs_plot, 
                                MACE_wo_HF_smcfcs_cs_plot, HF_smcfcs_cs_plot)

save(smcfcs_withKTR_cs_plots, file = here::here("output","r_objects","gfr-vs-outcomes_smcfcs_withKTR_cs_plots.rda"))


#  ----------------------------- CUBIC SPLINE TABLES ---------------------------

# KFRT
kfrt_smcfcs_cs_tbl <- summarize_HR_MICE(.predictor = predictors, 
                                        outcome_element = "rrt", outcome_name = "KFRT", 
                                        .covariates =  covariates, 
                                        .imp_data = smcfcs_rrt, 
                                        .ref = 90, .trunc = 120)

# AKI
aki_smcfcs_cs_tbl <- summarize_HR_MICE(.predictor = predictors, 
                                       outcome_element = "aki", outcome_name = "AKI", 
                                       .covariates =  covariates, 
                                       .imp_data = smcfcs_aki, 
                                       .ref = 90, .trunc = 120)

# Death
death_smcfcs_cs_tbl  <- predictors %>%
   map( ~ compute_fits_imp_data_death(.outcome = "death", 
                                .predictor = .x, 
                                .covariates = covariates, 
                                .imp_data = smcfcs_death)) |>
   map(~ extract_termplot_imp(.x, 90, 120)) |>
   map(~ pool_termplots(.x)) |>
   list_rbind() |>
   format_table(.outcome = "All-cause death")

# MACE without HF
MACE_wo_HF_smcfcs_cs_tbl <- summarize_HR_MICE(.predictor = predictors, 
                                              outcome_element = "MACE_wo_HF", outcome_name = "MACE", 
                                              .covariates =  covariates, 
                                              .imp_data = smcfcs_MACEwo, 
                                              .ref = 90, .trunc = 120)

# Heart failure
HF_smcfcs_cs_tbl <- summarize_HR_MICE(.predictor = predictors, 
                                      outcome_element = "heart_failure", outcome_name = "Heart Failure", 
                                      .covariates =  covariates[!covariates %in% c("hf")], 
                                      .imp_data = smcfcs_hf, 
                                      .ref = 90, .trunc = 120)

# Merge all summarized tables
smcfcs_withKTR_cs_tbl <- purrr::reduce(.x = list(kfrt_smcfcs_cs_tbl, aki_smcfcs_cs_tbl, death_smcfcs_cs_tbl, 
                                                  MACE_wo_HF_smcfcs_cs_tbl, HF_smcfcs_cs_tbl), 
                                      .f = function(x, y) {left_join(x, y, by = c("x1", "predictor"))}) %>%
   relocate(predictor, GFR = x1, `All-cause death`,  KFRT, AKI, MACE, `Heart Failure`)

# save them 
save(smcfcs_withKTR_cs_tbl, file = here::here("output","r_objects","gfr-vs-outcomes_smcfcs_withKTR_cs_tbl.rda"))
