################################################################################
# Project: mGFR and outcomes in SCREAM 
# Purpose: sensitivity analysis including individuals with a history of HF 
#          or excluding patients with a history of MACE/AKI
# Written by: Antoine Creon
# Date: 2025-03-21
################################################################################

################################################################################
###                    Prepare datasets and helper functions                 ###
################################################################################

source(here::here("code", "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"))

data_imp <- read_rds(file = here::here("data","cleaned","not_cens_kfrt_imp_cs.rds"))

# DON'T FILTER patients with history of heart failure

# FILTER patients without history of MACE (MI and stroke)
data_imp_wo_mace <- data_imp |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(history_mi == 0, history_stroke == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed

# FILTER patients without history of AKI
data_imp_wo_aki <- data_imp |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(history_aki == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed



################################################################################
###           INCLUDING INDIVIDUALS WITH HISTORY OF HEART FAILURE            ###
################################################################################

#  --------------------------------- Plot --------------------------------------

with_hf <- plot_imputed_outcomes(.predictor = predictors, 
                                outcome_element = "heart_failure", 
                                outcome_name = "Heart failure", 
                                .covariates = covariates, 
                                .imp_data = data_imp, 
                                .ref = 90, 
                                .trunc = 120) 


#  -------------------------------- Table --------------------------------------

with_hf_tbl <- summarize_HR_MICE(.predictor = predictors, 
                                 outcome_element = "heart_failure", 
                                 outcome_name = "Heart failure", 
                                 .covariates = covariates, 
                                 .imp_data = data_imp, 
                                 .ref = 90, 
                                 .trunc = 120) 


################################################################################
###           EXCLUDING INDIVIDUALS WITH HISTORY OF MACE / AKI               ###
################################################################################

#  --------------------------------- Plots -------------------------------------


# MACE excluding patients with history of stroke or MI
wo_MACE <- plot_imputed_outcomes(.predictor = predictors, 
                                 outcome_element = "MACE_wo_HF", 
                                 outcome_name = "MACE", 
                                 .covariates = covariates[!covariates %in% c("history_stroke", "history_mi")], 
                                 .imp_data = data_imp_wo_mace, 
                                 .ref = 90, 
                                 .trunc = 120) 

# AKI excluding patients with history of AKI
wo_AKI <- plot_imputed_outcomes(.predictor = predictors, 
                                 outcome_element = "aki", 
                                 outcome_name = "AKI", 
                                 .covariates = covariates[!covariates %in% c("history_aki")], 
                                 .imp_data = data_imp_wo_aki, 
                                 .ref = 90, 
                                 .trunc = 120) 


#  --------------------------------- Tables -------------------------------------

# MACE excluding patients with history of stroke or MI
wo_MACE_tbl <- summarize_HR_MICE(.predictor = predictors, 
                                 outcome_element = "MACE_wo_HF", 
                                 outcome_name = "MACE", 
                                 .covariates = covariates[!covariates %in% c("history_stroke", "history_mi")], 
                                 .imp_data = data_imp_wo_mace, 
                                 .ref = 90, 
                                 .trunc = 120) 

# AKI excluding patients with history of AKI
wo_AKI_tbl <- summarize_HR_MICE(.predictor = predictors, 
                                outcome_element = "aki", 
                                outcome_name = "AKI", 
                                .covariates = covariates[!covariates %in% c("history_aki")], 
                                .imp_data = data_imp_wo_aki, 
                                .ref = 90, 
                                .trunc = 120) 


################################################################################
###                                SAVE OBJECTS                              ###
################################################################################

#  --------------------------------- Plots -------------------------------------

# Append all lists
incident_or_prevalent_indiv_plot <- list("HF" = with_hf, "MACE" = wo_MACE, "AKI" = wo_AKI)

# Save them as R objects to be called and modified in analysis reports
save(incident_or_prevalent_indiv_plot, file = here::here("output","r_objects","incident_or_prevalent_indiv_plot.rda"))


#  --------------------------------- Tables ------------------------------------

save(with_hf_tbl, file = here::here("output","r_objects","with_hf_tbl.rda"))


incident_indiv_tbl <- left_join(wo_MACE_tbl, wo_AKI_tbl, by = c("x1", "predictor"))

save(incident_indiv_tbl, file = here::here("output","r_objects","incident_indiv_tbl.rda"))
