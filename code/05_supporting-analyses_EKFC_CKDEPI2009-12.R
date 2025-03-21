################################################################################
# Project: mGFR and outcomes in SCREAM 
# Purpose: code for supporting analysis using EKFC and CKD-EPI 2009 estimating equations
# Written by: Antoine Creon
# Date: 2025-02-06
################################################################################

################################################################################
###                              LOAD DATA                                   ###
################################################################################

# Load packages, helper functions and covariate/predictor lists
source(here::here("code", "06_gfr-versus-outcomes_uncensored-KFRT_analysis-preparation.R"))

# Load data with events NOT CENSORED AT KFRT
data <- read_rds(file = here::here("data","cleaned","data_not_cens_kfrt_named.rds")) 

# MICE-imputed data
data_imp <- read_rds(file = here::here("data","cleaned","not_cens_kfrt_imp_cs.rds"))


################################################################################
###                            PREPARE THE ANALYSIS                          ###
################################################################################


#  ------------------------ UNCHANGED HELPER FUNCTIONS -------------------------

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

plot_imputed_outcomes <- function(.predictor, outcome_element, outcome_name, .covariates, .imp_data, .ref, .trunc, .spline = "cubic", .linear_knots = "c(30,60,90)"){
   
   .predictor %>%
      map(~ compute_fits_imp_data(.outcome = outcome_element, 
                                  .predictor = .x, 
                                  .covariates, 
                                  .imp_data = .imp_data,
                                  .spline = .spline,
                                  .linear_knots = .linear_knots)) %>%
      map(~ extract_termplot_imp(.x, ref = .ref, truncate = .trunc)) %>%
      map(~ pool_termplots(.x)) %>%
      list_rbind() %>%
      plot_HRs(., .outcome = outcome_name)
}


summarize_HR_MICE <- function(.predictor, outcome_element, outcome_name, .covariates, .imp_data, .ref, .trunc, .spline = "cubic", .linear_knots = "c(30,60,90)"){
   
   .predictor %>%
      map(~ compute_fits_imp_data(.outcome = outcome_element, 
                                  .predictor = .x, 
                                  .covariates, 
                                  .imp_data = .imp_data,
                                  .spline = .spline,
                                  .linear_knots = .linear_knots)) %>%
      map(~ extract_termplot_imp(.x, ref = .ref, truncate = .trunc)) %>%
      map(~ pool_termplots(.x)) %>%
      list_rbind() %>%
      format_table(., .outcome = {{outcome_name}})
}


#  --------------------- ADD EKFC and CKDEPI 2009 TO DATASET -------------------


# Add EKFC and CKDEPI2009 eGFR to original dataset
data <- data |>
   mutate(ekfc_cr = ekfc_cr(creat, age, female),
          ekfc_cys = ekfc_cys(cystatin, age),
          ckd_epi_2009_cr = ckd_epi_2009_cr(creat, age, female),
          ckd_epi_2012_cr_cys = ckd_epi_2012_cr_cys(creat, cystatin, age, female)) |> 
   rowwise() |> 
   mutate(ekfc_combined = mean(c(ekfc_cr, ekfc_cys))) |>
   ungroup() |> 
   relocate(ekfc_cr, ekfc_combined, ekfc_cys)

# Add EKFC eGFR with cubic splines to the imputed dataset
data_imp_sensitivity <- data_imp |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   mutate(ekfc_cr = rep(data$ekfc_cr, 51),
          ekfc_cys = rep(data$ekfc_cys, 51),
          ekfc_combined = rep(data$ekfc_combined, 51),
          ckd_epi_2009_cr = rep(data$ckd_epi_2009_cr, 51),
          ckd_epi_2012_cr_cys = rep(data$ckd_epi_2012_cr_cys, 51)) |> 
   cbind(compute_cubic_splines(data, "ekfc_cr")) |> 
   cbind(compute_cubic_splines(data, "ekfc_cys")) |> 
   cbind(compute_cubic_splines(data, "ekfc_combined")) |> 
   cbind(compute_cubic_splines(data, "ckd_epi_2009_cr")) |> 
   cbind(compute_cubic_splines(data, "ckd_epi_2012_cr_cys")) |> 
   mice::as.mids() # reconstruct into a mids object if needed

# Imputed data without individuals with HF at baseline
data_imp_se_wo_hf <- data_imp_sensitivity |>
   complete(action = 'long', include = TRUE) |> # Extract all completed datasets into a long df
   filter(hf == 0) |> # Apply filter to each dataset
   mice::as.mids() # reconstruct into a mids object if needed



################################################################################
###                           EKFC eGFR distribution                         ###
################################################################################

data <- data |>
   mutate(ekfc_cr = ekfc_cr(creat, age, female),
          ekfc_cys = ekfc_cys(cystatin, age)) |> 
   rowwise() |> 
   mutate(ekfc_combined = mean(c(ekfc_cr, ekfc_cys))) |>
   ungroup() |> 
   relocate(ekfc_cr, ekfc_combined, ekfc_cys)


gfr_distrib_EKFC <- data |> 
   select(mgfr, ekfc_cr, ekfc_cys, ekfc_combined) |> 
   pivot_longer(cols = everything(), names_to = "Equation", values_to = "GFR") |> 
   mutate(Equation = factor(Equation, levels = c("mgfr", "ekfc_cr", "ekfc_cys", "ekfc_combined")),
          linesize = if_else(Equation == "mgfr", 1, 0))

.palette <- palette_okabe_ito(c(1,3,5,7))

.custom_labels <- c(expression("mGFR"),
                    expression(paste("eGFR"[EKFCcr])),
                    expression(paste("eGFR"[EKFCcys])),
                    expression(paste("eGFR"[EKFCcr-cys])))

gfr_distrib_EKFC_plot <- gfr_distrib_EKFC |> 
   ggplot(aes(x = GFR, y =  after_stat(density), color = Equation, fill = Equation, linewidth = factor(linesize))) +
   geom_density(alpha = 0.2, position="identity") +
   labs(x = expression("GFR, mL/min/1.73m"^2),
        y = "Density",
        tag = "A") +
   scale_color_manual(name = NULL, values = .palette, labels = .custom_labels) +
   scale_fill_manual(name = NULL, values = .palette, labels = .custom_labels) +
   scale_x_continuous(breaks = c(15, 30, 45, 60, 90, 120), limits = c(5,120)) +
   scale_linewidth_manual(name = NULL, values = c(0.5, 2)) +
   guides(linewidth = "none", color = "none", fill = "none") +
   theme_bw() +
   theme(aspect.ratio = 0.75,
         legend.position = "inside",
         legend.justification.inside = c(0, 1),# Position in upper-right corner
         legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),  # Black border, white background
         legend.box.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "mm"),
         legend.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "mm"),
         legend.text = element_text(size = 8),
         legend.key.size = unit(0.2, "cm"),            # Adjust the size of the legend keys
         legend.spacing = unit(0.1, "cm"),
         axis.title = element_text(size = 10),
         plot.title = element_text(size = 11),
         plot.margin = unit(c(0, 0, 0, 0), "cm"))

save(gfr_distrib_EKFC_plot, file = here::here("output","r_objects","eGFR-EKFC_distribution.rda"))



################################################################################
###                           EKFC EGFR VERSUS OUTCOMES                      ###
################################################################################

# ---------- MODIFY HELPER FUNCTIONS TO HANDLE THE NEW PREDICTORS --------------

# Extract termplots from fits on imputed datasets
## For each imputation
extract_each_termplotof_imp <- function(fit,ref,truncate) {
   ptemp <- termplot(fit, se=TRUE, plot = FALSE)
   x1 <- ptemp[[1]][["x"]]
   y1 <- ptemp[[1]][["y"]]
   yref1 <- -y1[abs(x1-ref)==min(abs(x1-ref))]
   se1 <- ptemp[[1]][["se"]]
   
   dat1 <- as.data.frame(cbind(x1,y1,se1)) %>% 
      filter(x1 <= truncate) 
   
   dat1 <- dat1 %>% 
      mutate(y1 = y1 + yref1, # dont't compute HR and CI right away (pooling needed first)
             var1 = se1^2) %>%  # compute var instead of SE
      select(-se1)
   
   # Only include data where x1 is greater than or equal to 15
   dat1 <- dat1[dat1$x1 >= 15, ]
   
   # Specify the predictor and name it to suit the graph
   dat1 <- dat1 %>%
      mutate(predictor = names(ptemp)[[1]],
             predictor = factor(predictor, 
                                levels = c("mgfr_s1", "ekfc_cr_s1", "ekfc_cys_s1", "ekfc_combined_s1"),
                                labels = c("mGFR", 
                                           "eGFR[EKFCcr]",
                                           "eGFR[EKFCcys]",
                                           "eGFR[EKFCcr-cys]")))
   return(dat1)
}

## Second, map the previous function to all elements of the list
extract_termplot_imp <- function(.imp_fits, ref, truncate) {
   
   .imp_fits %>%
      map(~ extract_each_termplotof_imp(.x, ref, truncate))
}

# Plot the results
plot_HRs <- function(.extracted_termplots, .outcome) {
   
   .palette <- palette_okabe_ito(c(1,3,5,7))
   # palette <- c("#D55E00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2")
   
   .custom_labels <- c(expression("mGFR"),
                       expression(paste("eGFR"[EKFCcr])),
                       expression(paste("eGFR"[EKFCcys])),
                       expression(paste("eGFR"[EKFCcr-cys])))
   
   .plot <- .extracted_termplots %>%
      mutate(bold_line = case_when(predictor == "mGFR" ~ 1,
                                   TRUE ~ 0)) %>%
      ggplot(aes(x = x1, y = y1, group = predictor, color = predictor, fill = predictor, size = bold_line)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper, group = predictor), alpha = 0.1, color = NA) +     # area inside 95% CI
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
      labs(x = expression("GFR, mL/min/1.73m"^2),
           y = paste0("Hazard ratio for ", .outcome),
           fill = "", color = "", group = "",
           title = paste0(.outcome)) +
      guides(size = "none") + # remove the legend for the thickest line
      theme_bw() + 
      theme(aspect.ratio = 0.75,
            legend.position = "inside",
            legend.justification.inside = c(1, 1),# Position in upper-right corner
            legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),  # Black border, white background
            legend.box.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "mm"),
            legend.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "mm"),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.2, "cm"),            # Adjust the size of the legend keys
            legend.spacing = unit(0.1, "cm"),
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 11),
            plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      scale_color_manual(name = NULL, values = .palette, labels = .custom_labels) +  
      scale_fill_manual(name = NULL, values = .palette, labels = .custom_labels) +
      scale_size_continuous(range = c(0.5, 2)) +
      scale_x_continuous(breaks = c(15, 30, 45, 60, 90, 120))
   
   return(.plot)
}


#  --------------------------- PERFORM THE ANALYSIS ----------------------------

predictors_ekfc = list("mgfr", "ekfc_cr", "ekfc_cys", "ekfc_combined")


# All events except heart failure
ref90_EKFC_plots <- events[!events %in% c("heart_failure")] %>%
   imap(~ plot_imputed_outcomes(.predictor = predictors_ekfc, 
                                outcome_element = .x, 
                                outcome_name = .y, 
                                .covariates = covariates, 
                                .imp_data = data_imp_sensitivity, 
                                .ref = 90, 
                                .trunc = 120)
   )

hf_EKFC_plot <- plot_imputed_outcomes(.predictor = predictors_ekfc, 
                                      outcome_element = "heart_failure", 
                                      outcome_name = "Heart failure", 
                                      .covariates = covariates[!covariates %in% c("hf")], 
                                      .imp_data = data_imp_se_wo_hf, 
                                      .ref = 90, 
                                      .trunc = 120) 

# Append all lists
EKFC_plots <- append(x = ref90_EKFC_plots, list("HF" = hf_EKFC_plot))

# Save them as R objects to be called and modified in analysis reports
save(EKFC_plots, file = here::here("output","r_objects","EKFC_plots.rda"))

# Summarize HRs in tables
ref90_EKFC_tbl <- events[!events %in% c("heart_failure")] %>%
   map2(.x = .,
        .y = names(.),
        ~ summarize_HR_MICE(.predictor = predictors_ekfc, 
                            outcome_element = .x, 
                            outcome_name = {{.y}}, 
                            .covariates = covariates, 
                            .imp_data = data_imp_sensitivity, 
                            .ref = 90, 
                            .trunc = 120)
   )

hf_EKFC_tbl <- summarize_HR_MICE(.predictor = predictors_ekfc, 
                                 outcome_element = "heart_failure", 
                                 outcome_name = "Heart failure",
                                 .covariates = covariates[!covariates %in% c("hf")], 
                                 .imp_data = data_imp_se_wo_hf, 
                                 .ref = 90, 
                                 .trunc = 120)


# merge all summarized tables
EKFC_tbl <- purrr::reduce(.x = append(x = ref90_EKFC_tbl, list("HF" = hf_EKFC_tbl)),
                          .f = function(x, y) {left_join(x, y, by = c("x1", "predictor"))}) %>%
   relocate(predictor, GFR = x1, `All-cause Death`,  KFRT, AKI, MACE, `Heart failure`)


# Save  as R object to be called and modified in analysis reports
save(EKFC_tbl, file = here::here("output","r_objects","EKFC_tbl.rda"))



################################################################################
###                       CKD-EPI 2009-12 eGFR distribution                  ###
################################################################################


gfr_distrib_CKDEPI2009 <- data |> 
   select(mgfr, ckd_epi_2009_cr, ckd_epi_2012_cys, ckd_epi_2012_cr_cys) |> 
   pivot_longer(cols = everything(), names_to = "Equation", values_to = "GFR") |> 
   mutate(Equation = factor(Equation, levels = c("mgfr", "ckd_epi_2009_cr", "ckd_epi_2012_cys", "ckd_epi_2012_cr_cys")),
          linesize = if_else(Equation == "mgfr", 1, 0))

.palette <- palette_okabe_ito(c(1,3,5,7))

.custom_labels <- c(expression("mGFR"),
                    expression(paste("eGFR"[CKDEPI2009-cr])),
                    expression(paste("eGFR"[CKDEPI2009-cys])),
                    expression(paste("eGFR"[CKDEPI2009-cr-cys])))

gfr_distrib_CKDEPI2009_plot <- gfr_distrib_CKDEPI2009 |> 
   ggplot(aes(x = GFR, y =  after_stat(density), color = Equation, fill = Equation, linewidth = factor(linesize))) +
   geom_density(alpha = 0.2, position="identity") +
   labs(x = expression("GFR, mL/min/1.73m"^2),
        y = "Density",
        tag = "A") +
   scale_color_manual(name = NULL, values = .palette, labels = .custom_labels) +
   scale_fill_manual(name = NULL, values = .palette, labels = .custom_labels) +
   scale_x_continuous(breaks = c(15, 30, 45, 60, 90, 120), limits = c(5,120)) +
   scale_linewidth_manual(name = NULL, values = c(0.5, 2)) +
   guides(linewidth = "none", color = "none", fill = "none") +
   theme_bw() +
   theme(aspect.ratio = 0.75,
         legend.position = "inside",
         legend.justification.inside = c(0, 1),# Position in upper-right corner
         legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),  # Black border, white background
         legend.box.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "mm"),
         legend.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "mm"),
         legend.text = element_text(size = 8),
         legend.key.size = unit(0.2, "cm"),            # Adjust the size of the legend keys
         legend.spacing = unit(0.1, "cm"),
         axis.title = element_text(size = 10),
         plot.title = element_text(size = 11),
         plot.margin = unit(c(0, 0, 0, 0), "cm"))

save(gfr_distrib_CKDEPI2009_plot, file = here::here("output","r_objects","eGFR-CKDEPI2009_distribution.rda"))



################################################################################
###                    CKDEPI 2009-12 EGFR VERSUS OUTCOMES                   ###
################################################################################

# ---------- MODIFY HELPER FUNCTIONS TO HANDLE THE NEW PREDICTORS --------------

# Extract termplots from fits on imputed datasets
## For each imputation
extract_each_termplotof_imp <- function(fit,ref,truncate) {
   ptemp <- termplot(fit, se=TRUE, plot = FALSE)
   x1 <- ptemp[[1]][["x"]]
   y1 <- ptemp[[1]][["y"]]
   yref1 <- -y1[abs(x1-ref)==min(abs(x1-ref))]
   se1 <- ptemp[[1]][["se"]]
   
   dat1 <- as.data.frame(cbind(x1,y1,se1)) %>% 
      filter(x1 <= truncate) 
   
   dat1 <- dat1 %>% 
      mutate(y1 = y1 + yref1, # dont't compute HR and CI right away (pooling needed first)
             var1 = se1^2) %>%  # compute var instead of SE
      select(-se1)
   
   # Only include data where x1 is greater than or equal to 15
   dat1 <- dat1[dat1$x1 >= 15, ]
   
   # Specify the predictor and name it to suit the graph
   dat1 <- dat1 %>%
      mutate(predictor = names(ptemp)[[1]],
             predictor = factor(predictor, 
                                levels = c("mgfr_s1", "ckd_epi_2009_cr_s1", "ckd_epi_2012_cys_s1", "ckd_epi_2012_cr_cys_s1"),
                                labels = c("mGFR", 
                                           "eGFR[CKDEPI2009-cr]",
                                           "eGFR[CKDEPI2012-cys]",
                                           "eGFR[CKDEPI2012-cr-cys]")))
   return(dat1)
}

## Second, map the previous function to all elements of the list
extract_termplot_imp <- function(.imp_fits, ref, truncate) {
   
   .imp_fits %>%
      map(~ extract_each_termplotof_imp(.x, ref, truncate))
}

# Plot the results
plot_HRs <- function(.extracted_termplots, .outcome) {
   
   .palette <- palette_okabe_ito(c(1,3,5,7))
   # palette <- c("#D55E00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2")
   
   .custom_labels <- c(expression("mGFR"),
                       expression(paste("eGFR"[CKDEPI2009-cr])),
                       expression(paste("eGFR"[CKDEPI2012-cys])),
                       expression(paste("eGFR"[CKDEPI2012-cr-cys])))
   
   .plot <- .extracted_termplots %>%
      mutate(bold_line = case_when(predictor == "mGFR" ~ 1,
                                   TRUE ~ 0)) %>%
      ggplot(aes(x = x1, y = y1, group = predictor, color = predictor, fill = predictor, size = bold_line)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper, group = predictor), alpha = 0.1, color = NA) +     # area inside 95% CI
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
      labs(x = expression("GFR, mL/min/1.73m"^2),
           y = paste0("Hazard ratio for ", .outcome),
           fill = "", color = "", group = "",
           title = paste0(.outcome)) +
      guides(size = "none") + # remove the legend for the thickest line
      theme_bw() + 
      theme(aspect.ratio = 0.75,
            legend.position = "inside",
            legend.justification.inside = c(1, 1),# Position in upper-right corner
            legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),  # Black border, white background
            legend.box.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "mm"),
            legend.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "mm"),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.2, "cm"),            # Adjust the size of the legend keys
            legend.spacing = unit(0.1, "cm"),
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 11),
            plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      scale_color_manual(name = NULL, values = .palette, labels = .custom_labels) +  
      scale_fill_manual(name = NULL, values = .palette, labels = .custom_labels) +
      scale_size_continuous(range = c(0.5, 2)) +
      scale_x_continuous(breaks = c(15, 30, 45, 60, 90, 120))
   
   return(.plot)
}


#  --------------------------- PERFORM THE ANALYSIS ----------------------------

predictors_ckdepi2009 <-  list("mgfr", "ckd_epi_2009_cr", "ckd_epi_2012_cys", "ckd_epi_2012_cr_cys")

# All events except heart failure
ref90_CKDEPI2009_plots <- events[!events %in% c("heart_failure")] %>%
   imap(~ plot_imputed_outcomes(.predictor = predictors_ckdepi2009, 
                                outcome_element = .x, 
                                outcome_name = .y, 
                                .covariates = covariates, 
                                .imp_data = data_imp_sensitivity, 
                                .ref = 90, 
                                .trunc = 120)
   )

hf_CKDEPI2009_plot <- plot_imputed_outcomes(.predictor = predictors_ckdepi2009, 
                                            outcome_element = "heart_failure", 
                                            outcome_name = "Heart failure", 
                                            .covariates = covariates[!covariates %in% c("hf")], 
                                            .imp_data = data_imp_se_wo_hf, 
                                            .ref = 90, 
                                            .trunc = 120) 

# Append all lists
CKDEPI2009_plots <- append(x = ref90_CKDEPI2009_plots, list("HF" = hf_CKDEPI2009_plot))

# Save them as R objects to be called and modified in analysis reports
save(CKDEPI2009_plots, file = here::here("output","r_objects","CKDEPI2009_plots.rda"))

# Summarize HRs in tables
ref90_CKDEPI2009_tbl <- events[!events %in% c("heart_failure")] %>%
   map2(.x = .,
        .y = names(.),
        ~ summarize_HR_MICE(.predictor = predictors_ckdepi2009, 
                            outcome_element = .x, 
                            outcome_name = {{.y}}, 
                            .covariates = covariates, 
                            .imp_data = data_imp_sensitivity, 
                            .ref = 90, 
                            .trunc = 120)
   )

hf_CKDEPI2009_tbl <- summarize_HR_MICE(.predictor = predictors_ckdepi2009, 
                                       outcome_element = "heart_failure", 
                                       outcome_name = "Heart failure",
                                       .covariates = covariates[!covariates %in% c("hf")], 
                                       .imp_data = data_imp_se_wo_hf, 
                                       .ref = 90, 
                                       .trunc = 120)


# merge all summarized tables
CKDEPI2009_tbl <- purrr::reduce(.x = append(x = ref90_CKDEPI2009_tbl, list("HF" = hf_CKDEPI2009_tbl)),
                                .f = function(x, y) {left_join(x, y, by = c("x1", "predictor"))}) %>%
   relocate(predictor, GFR = x1, `All-cause Death`,  KFRT, AKI, MACE, `Heart failure`)


# Save  as R object to be called and modified in analysis reports
save(CKDEPI2009_tbl, file = here::here("output","r_objects","CKDEPI2009_tbl.rda"))
