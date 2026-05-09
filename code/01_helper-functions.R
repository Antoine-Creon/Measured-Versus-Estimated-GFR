################################################################################
# Project: mGFR and outcomes in SCREAM
# Purpose: Code for helper functions to perform the survival analyses
# Written by: Antoine Creon
# Date: 2024-11-06
################################################################################

################################################################################
# eGFR EQUATIONS ###############################################################
################################################################################

ckd_epi_2009_cr <- function(creatinine, age, female) {
   k <- ifelse(female == 1, 62, 80) # this differs from the original formula because we have in umol/L, not mg/dL
   alpha <- ifelse(female == 1, -0.329, -0.411)
   return(ifelse(
      female == 1,
      141 *
         (pmin(creatinine / k, 1)^alpha) *
         (pmax(creatinine / k, 1)^(-1.209)) *
         (0.9929^age) *
         1.018,
      141 *
         (pmin(creatinine / k, 1)^alpha) *
         (pmax(creatinine / k, 1)^(-1.209)) *
         (0.9929^age)
   ))
}

ckd_epi_2012_cys <- function(cystatin, age, female) {
   return(ifelse(
      female == 1,
      133 *
         (pmin(cystatin / 0.8, 1)^(-0.499)) *
         (pmax(cystatin / 0.8, 1)^(-1.328)) *
         (0.9962^age) *
         0.932,
      133 *
         (pmin(cystatin / 0.8, 1)^(-0.499)) *
         (pmax(cystatin / 0.8, 1)^(-1.328)) *
         (0.9962^age)
   ))
}

ckd_epi_2021_cr <- function(creatinine, age, female) {
   k <- ifelse(female == 1, 62, 80)
   alpha <- ifelse(female == 1, -0.241, -0.302)
   return(ifelse(
      female == 1,
      142 *
         (pmin(creatinine / k, 1)^alpha) *
         (pmax(creatinine / k, 1)^(-1.200)) *
         (0.9938^age) *
         1.012,
      142 *
         (pmin(creatinine / k, 1)^alpha) *
         (pmax(creatinine / k, 1)^(-1.200)) *
         (0.9938^age)
   ))
}

ckd_epi_2012_cr_cys <- function(creatinine, cystatin, age, female) {
   k <- ifelse(female == 1, 62, 80)
   alpha <- ifelse(female == 1, -0.248, -0.207)
   return(ifelse(
      female == 1,
      135 *
         (pmin(creatinine / k, 1)^alpha) *
         (pmax(creatinine / k, 1)^(-0.601)) *
         (pmin(cystatin / 0.8, 1)^(-0.375)) *
         (pmax(cystatin / 0.8, 1)^(-0.711)) *
         (0.9952^age) *
         0.969,
      135 *
         (pmin(creatinine / k, 1)^alpha) *
         (pmax(creatinine / k, 1)^(-0.601)) *
         (pmin(cystatin / 0.8, 1)^(-0.375)) *
         (pmax(cystatin / 0.8, 1)^(-0.711)) *
         (0.9952^age)
   ))
}

ckd_epi_2021_cr_cys <- function(creatinine, cystatin, age, female) {
   k <- ifelse(female == 1, 62, 80)
   alpha <- ifelse(female == 1, -0.219, -0.144)
   return(ifelse(
      female == 1,
      135 *
         (pmin(creatinine / k, 1)^alpha) *
         (pmax(creatinine / k, 1)^(-0.544)) *
         (pmin(cystatin / 0.8, 1)^(-0.323)) *
         (pmax(cystatin / 0.8, 1)^(-0.778)) *
         (0.9961^age) *
         0.963,
      135 *
         (pmin(creatinine / k, 1)^alpha) *
         (pmax(creatinine / k, 1)^(-0.544)) *
         (pmin(cystatin / 0.8, 1)^(-0.323)) *
         (pmax(cystatin / 0.8, 1)^(-0.778)) *
         (0.9961^age)
   ))
}

ekfc_cr <- function(creatinine, age, female) {
   Q <- ifelse(
      age <= 25 & female == 0,
      exp(
         3.200 +
            0.259 * age -
            0.543 * log(age) -
            0.00763 * age^2 +
            0.0000790 * age^3
      ),
      ifelse(
         age <= 25 & female == 1,
         exp(
            3.080 +
               0.177 * age -
               0.223 * log(age) -
               0.00596 * age^2 +
               0.0000686 * age^3
         ),
         ifelse(
            age > 25 & female == 0,
            80,
            ifelse(age > 25 & female == 1, 62, NA)
         )
      )
   )

   return(ifelse(
      creatinine / Q < 1 & age <= 40,
      107.3 * (creatinine / Q)^-0.322,
      ifelse(
         creatinine / Q < 1 & age > 40,
         107.3 * (creatinine / Q)^-0.322 * 0.990^(age - 40),
         ifelse(
            creatinine / Q >= 1 & age <= 40,
            107.3 * (creatinine / Q)^-1.132,
            ifelse(
               creatinine / Q >= 1 & age > 40,
               107.3 * (creatinine / Q)^-1.132 * 0.990^(age - 40),
               NA
            )
         )
      )
   ))
}

ekfc_cys <- function(cystatin, age) {
   Q <- if_else(age > 50, 0.83 + 0.005 * (age - 50), 0.83)

   return(if_else(
      age >= 18 & age <= 40 & cystatin / Q < 1,
      107.3 * (cystatin / Q)^-0.322,
      if_else(
         age >= 18 & age <= 40 & cystatin / Q >= 1,
         107.3 * (cystatin / Q)^-1.132,
         if_else(
            age > 40 & age <= 50 & cystatin / Q < 1,
            107.3 * (cystatin / Q)^-0.322 * 0.990^(age - 40),
            if_else(
               age > 40 & age <= 50 & cystatin / Q >= 1,
               107.3 * (cystatin / Q)^-1.132 * 0.990^(age - 40),
               if_else(
                  age > 50 & cystatin / Q < 1,
                  107.3 * (cystatin / Q)^-0.322 * 0.990^(age - 40),
                  if_else(
                     age > 50 & cystatin / Q >= 1,
                     107.3 * (cystatin / Q)^-1.132 * 0.990^(age - 40),
                     NA
                  )
               )
            )
         )
      )
   ))
}

################################################################################
# TO PLOT OR MAKE TABLES  ######################################################
################################################################################

## to plot each model  ---------------------------------------------------------
# For one outcome, we want the 4 predictors to be on the same graph.
# the input is an extracted termplot as given by the functions below

plot_HRs <- function(.extracted_termplots, .outcome) {
   .palette <- palette_okabe_ito(c(1, 3, 5, 7))
   # palette <- c("#D55E00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2")

   .custom_labels <- c(
      expression("mGFR"),
      expression(paste("eGFR"[cr])),
      expression(paste("eGFR"[cys])),
      expression(paste("eGFR"[cr - cys]))
   )

   .plot <- .extracted_termplots %>%
      mutate(bold_line = case_when(predictor == "mGFR" ~ 1, TRUE ~ 0)) %>%
      ggplot(aes(
         x = x1,
         y = y1,
         group = predictor,
         color = predictor,
         fill = predictor,
         linewidth = factor(bold_line)
      )) +
      geom_line() +
      geom_ribbon(
         aes(ymin = lower, ymax = upper, group = predictor),
         alpha = 0.1,
         color = NA
      ) + # area inside 95% CI
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
      labs(
         x = expression("Glomerular filtration rate, mL/min/1.73m"^2),
         y = "Adjusted hazard ratio",
         fill = "",
         color = "",
         group = "",
         title = paste0(.outcome)
      ) +
      guides(size = "none") + # remove the legend for the thickest line
      theme_bw() +
      theme(
         aspect.ratio = 0.75,
         legend.position = "inside",
         legend.justification.inside = c(1, 1), # Position in upper-right corner
         legend.background = element_rect(
            color = "black",
            fill = "white",
            linewidth = 0.1
         ), # Black border, white background
         legend.box.margin = margin(
            t = 0.5,
            r = 0.5,
            b = 0.5,
            l = 0.5,
            unit = "mm"
         ),
         legend.margin = margin(
            t = 0.5,
            r = 0.5,
            b = 0.5,
            l = 0.5,
            unit = "mm"
         ),
         legend.text = element_text(size = 8),
         legend.key.size = unit(0.2, "cm"), # Adjust the size of the legend keys
         legend.spacing = unit(0.1, "cm"),
         axis.title = element_text(size = 10),
         plot.title = element_text(size = 11),
         plot.margin = unit(c(0, 0, 0, 0), "cm")
      ) +
      scale_color_manual(
         name = NULL,
         values = .palette,
         labels = .custom_labels
      ) +
      scale_fill_manual(
         name = NULL,
         values = .palette,
         labels = .custom_labels
      ) +
      scale_linewidth_manual(name = NULL, values = c(0.5, 1)) +
      scale_x_continuous(breaks = c(15, 30, 45, 60, 75, 90, 120))

   return(.plot)
}


## to summarize each model  ----------------------------------------------------
# the input is an extracted termplot as given by the functions below
# a column for GFR level and each outcome
# a row for 15-30-45-45-60 ml/min
# 4 groups of rows cr, cys, cys-cr and mGFR

# A. intermediate function to make a tibble from an extracted termplot
format_table <- function(.bound_termplots, .outcome) {
   .bound_termplots %>%
      mutate(
         x1 = scales::number(x1, accuracy = 1),
         across(
            all_of(c("y1", "upper", "lower")),
            ~ scales::number(.x, accuracy = 0.01)
         )
      ) %>%
      group_by(predictor, x1) %>%
      slice(1) %>%
      filter(
         x1 == 15 |
            x1 == 30 |
            x1 == 45 |
            x1 == 60 |
            x1 == 75 |
            x1 == 90 |
            x1 == 120
      ) %>%
      ungroup() %>%
      mutate({{ .outcome }} := glue::glue("{y1} ({lower} - {upper})")) %>%
      select(x1, predictor, {{ .outcome }})
}

################################################################################
# TO DEAL WITH STANDARD (not multiply imputed) DATASETS ########################
################################################################################

## to compute the fits of specific predictors and covariates -------------------
# The function takes the outcome, the predictor, the covariates and the data
# It returns the coxph fit
# Spline option: linear or cubic (default)

compute_fits <- function(
   outcome,
   predictor,
   covariates,
   data,
   .spline = "cubic",
   .linear_knots = "c(30, 60, 90)"
) {
   # time to event variable
   .time <- paste0("time_to_", outcome)

   # string of covariates
   .covariates <- paste(covariates, collapse = " + ")

   # Range of the predictor, for boundary knots of the linear splines
   .range_predictor <- paste0(range(data[[predictor]]), collapse = ",")

   # Surv function and ~, depending on the type of spline chosen
   if (.spline == "linear") {
      .surv <- paste0(
         "Surv(",
         .time,
         ",",
         outcome,
         ") ~ splines::bs(",
         predictor,
         ", knots =",
         .linear_knots,
         ", degree = 1) +"
      )
   } else if (.spline == "cubic") {
      .surv <- paste0(
         "Surv(",
         .time,
         ",",
         outcome,
         ") ~ rcs(",
         predictor,
         ",4) +"
      )
   } else {
      rlang::abort("spline must be either 'linear' or 'cubic'")
   }

   # final fit
   .fit <- coxph(as.formula(paste0(.surv, .covariates)), data = data)
   return(.fit)
}


## to extract the termplot of the predictor from the coxph fit -----------------

extract_termplot <- function(fit, ref, truncate) {
   ptemp <- termplot(fit, se = TRUE, plot = FALSE)
   x1 <- ptemp[[1]][["x"]]
   y1 <- ptemp[[1]][["y"]]
   yref1 <- -y1[abs(x1 - ref) == min(abs(x1 - ref))]
   se1 <- ptemp[[1]][["se"]]

   dat1 <- as.data.frame(cbind(x1, y1, se1)) %>%
      filter(x1 <= truncate)

   dat1 <- dat1 %>%
      mutate(
         y1 = y1 + yref1,
         upper = exp(y1 + 1.96 * se1),
         lower = exp(y1 - 1.96 * se1),
         y1 = exp(y1)
      )

   # Only include data where x1 is greater than or equal to 15
   dat1 <- dat1[dat1$x1 >= 15, ]

   # Specify the predictor and name it to suit the graph
   dat1 <- dat1 %>%
      mutate(
         predictor = names(ptemp)[[1]],
         predictor = factor(
            predictor,
            levels = c(
               "mgfr",
               "ckd_epi_2021_cr",
               "ckd_epi_2012_cys",
               "ckd_epi_2021_cr_cys"
            ),
            labels = c("mGFR", "eGFR[cr]", "eGFR[cys]", "eGFR[cr-cys]")
         )
      )
   return(dat1)
}


## to extract fit + termplot + plot them  --------------------------------------
# To be used with a named list of outcomes + imap() or map2(outcome, names(outcome))

plot_outcomes <- function(
   .predictor,
   outcome_element,
   outcome_clean,
   .covariates,
   .data,
   .ref,
   .trunc,
   .spline = "cubic"
) {
   .predictor %>%
      map(~ compute_fits(outcome_element, .x, .covariates, .data, .spline)) %>%
      map(~ extract_termplot(.x, ref = .ref, .trunc)) %>%
      list_rbind() %>%
      plot_HRs(., .outcome = outcome_clean)
}


## to extract fit + termplot + summarize HRs  ----------------------------------

summarize_HR_CCA <- function(
   .predictor,
   outcome_element,
   outcome_name,
   .covariates,
   .data,
   .ref,
   .trunc,
   .spline = "cubic"
) {
   .predictor %>%
      map(~ compute_fits(outcome_element, .x, .covariates, .data, .spline)) %>%
      map(~ extract_termplot(.x, ref = .ref, .trunc)) %>%
      list_rbind() %>%
      format_table(., .outcome = {{ outcome_name }})
}


################################################################################
# FOR IMPUTED DATASETS #########################################################
################################################################################

## to compute the fits of predictors and covariates in each imputed dataset ----

compute_fits_imp_data <- function(
   .outcome,
   .predictor,
   .covariates,
   .imp_data,
   .spline = "cubic",
   .linear_knots = "c(30, 60, 90)",
   .weights = NULL
) {
   # name of the predictor as in imputed data
   .imp_predictor <- paste0(.predictor, "_s1")

   # Name the variables of the model depending on the imputation method
   if (class(.imp_data) == "mids") {
      .time <- paste0("time_to_", .outcome) # time to event variable
      .range_predictor <- paste0(
         range(.imp_data[["data"]][[.imp_predictor]]),
         collapse = ","
      )
   } else if (class(.imp_data) == "imputationList") {
      if (.outcome == "death") {
         .time <- paste0("time_to_", "MACE_wo_HF", "_compet")
         .outcome <- paste0("MACE_wo_HF", "_compet == 2") # for death outcome, use the MACE dataset
      } else {
         .time <- paste0("time_to_", .outcome, "_compet") # time to event variable
         .outcome <- paste0(.outcome, "_compet == 1") # outcome as written in SMC FCS
      }

      .range_predictor <- paste0(
         range(.imp_data[["imputations"]][[1]][[.imp_predictor]]),
         collapse = ","
      )
   } else {
      rlang::abort(
         ".imp_data must be either of class 'mids' (for MICE) or 'imputationList' (for SMC-FCS)"
      )
   }

   # string of covariates
   .covariates <- paste(.covariates, collapse = " + ")

   # Surv function and ~
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

   # Compute the fit in each imputed dataset
   .fits <- with(
      .imp_data,
      coxph(as.formula(paste0(.surv, .covariates)), weights = .weights)
   )

   # keep only each analysis
   if (class(.imp_data) == "mids") {
      .imp_fits <- .fits$analyses
      return(.imp_fits)
   } else if (class(.imp_data) == "imputationList") {
      return(.fits)
   } else {
      rlang::abort(
         ".imp_data must be either of class 'mids' (for MICE) or 'imputationList' (for SMC-FCS)"
      )
   }
}


## to extract the termplot of the predictor from the coxph fit -----------------
# The difference with the first function is that
# The function extract the termplot of each imputed dataset
# HR and CI are not computed (must be pooled first)
# and the names are changed to hamdle the inputs

# First extract the termplot of one dataset
extract_each_termplotof_imp <- function(fit, ref, truncate) {
   ptemp <- termplot(fit, se = TRUE, plot = FALSE)
   x1 <- ptemp[[1]][["x"]]
   y1 <- ptemp[[1]][["y"]]
   yref1 <- -y1[abs(x1 - ref) == min(abs(x1 - ref))]
   se1 <- ptemp[[1]][["se"]]

   dat1 <- as.data.frame(cbind(x1, y1, se1)) %>%
      filter(x1 <= truncate)

   dat1 <- dat1 %>%
      mutate(
         y1 = y1 + yref1, # dont't compute HR and CI right away (pooling needed first)
         var1 = se1^2
      ) %>% # compute var instead of SE
      select(-se1)

   # Only include data where x1 is greater than or equal to 15
   dat1 <- dat1[dat1$x1 >= 15, ]

   # Specify the predictor and name it to suit the graph
   dat1 <- dat1 %>%
      mutate(
         predictor = names(ptemp)[[1]],
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

# Second, map the previous function to all elements of the list
extract_termplot_imp <- function(.imp_fits, ref, truncate) {
   .imp_fits %>%
      map(~ extract_each_termplotof_imp(.x, ref, truncate))
}

## to pool the termplots from each imputed dataset -----------------------------
# First all data set must be col bound
# so that for each value of x1, the estimate and its SE can be pooled
# It needs to work rowwise
# Rubin's rules :
# pooled estimate = sum(Q)/m ; Q the estimate in each imputation, m nb of imputations
# Within-imputations variance W = sum(se^2)/m ; se the standard error of the estimate in each imputation
# Between-imputations variance B = Var(Q)
# Total variance of the pooled estimate: V = W + B + B/m
# Then HR = exp(coeff) and CI = exp(coeff +/- 1.96*se)

pool_termplots <- function(.extracted_termplots) {
   .pooled_termplots <- .extracted_termplots %>%
      map(~ as_tibble(.x)) %>%
      reduce(left_join, by = c("x1", "predictor")) %>% # merge imputed datasets
      rowwise(predictor, x1) %>%
      mutate(
         yhat = mean(c_across(starts_with("y1"))), # Pooled estimate: Q/m
         within_imp_var = mean(c_across(starts_with("var1"))), # Within-imputations variance W = se^2/m
         between_imp_var = var(c_across(starts_with("y1"))), # Between-imputations variance B = Var(Q)
         total_var = within_imp_var +
            between_imp_var +
            between_imp_var / length(.extracted_termplots), # Total variance V = W + B + B/m
         total_se = sqrt(total_var), # SE = sqrt(Var)
         HR_centered = exp(yhat), # HR = exp(coeff)
         HR_CI_high = exp(yhat + 1.96 * total_se),
         HR_CI_low = exp(yhat - 1.96 * total_se)
      ) %>%
      select(
         x1,
         y1 = HR_centered,
         lower = HR_CI_low,
         upper = HR_CI_high,
         predictor
      )

   return(.pooled_termplots)
}


## to do all steps above + plot ------------------------------------------------

plot_imputed_outcomes <- function(
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
   .predictor %>%
      map(
         ~ compute_fits_imp_data(
            .outcome = outcome_element,
            .predictor = .x,
            .covariates,
            .imp_data = .imp_data,
            .spline = .spline,
            .linear_knots = .linear_knots,
            .weights = .weights
         )
      ) %>%
      map(~ extract_termplot_imp(.x, ref = .ref, truncate = .trunc)) %>%
      map(~ pool_termplots(.x)) %>%
      list_rbind() %>%
      plot_HRs(., .outcome = outcome_name)
}


## to extract fit + termplot + pool summarize HRs  -----------------------------

summarize_HR_MICE <- function(
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
   .predictor %>%
      map(
         ~ compute_fits_imp_data(
            .outcome = outcome_element,
            .predictor = .x,
            .covariates,
            .imp_data = .imp_data,
            .spline = .spline,
            .linear_knots = .linear_knots,
            .weights = .weights
         )
      ) %>%
      map(~ extract_termplot_imp(.x, ref = .ref, truncate = .trunc)) %>%
      map(~ pool_termplots(.x)) %>%
      list_rbind() %>%
      format_table(., .outcome = {{ outcome_name }})
}


################################################################################
# CONDITIONAL INCIDENCE RATES ##################################################
################################################################################

cond_IR <- function(
   imputed_dataset,
   .outcome,
   .predictor,
   .covariates,
   .thresholds_to_estimate = c(15, 30, 45, 60, 75, 90, 120)
) {
   ## Step 0: Prepare variables -------------------------------------------------

   # Time to oucome
   .time_to_outcome <- paste0("time_to_", .outcome)

   # name of the predictor as in imputed data
   .imp_predictor <- paste0(.predictor, "_s1")

   # string of covariates
   .covariates_string <- paste(.covariates, collapse = " + ")

   ## Step 1: Model fitting in MI datasets and pool coefficients -----------------
   # /!\ rms::Glm requires the offset to be passed in the formula
   .poisson_formula <- paste0(
      .outcome,
      " == 1 ~ rcs(",
      .imp_predictor,
      ",4) + ",
      .covariates_string,
      " + offset(log(",
      .time_to_outcome,
      "))"
   )

   fit_mi <- with(
      imputed_dataset,
      rms::Glm(
         formula = as.formula(.poisson_formula),
         family = poisson(link = "log")
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
      uncount(length(.thresholds_to_estimate)) |>
      mutate(
         !!.imp_predictor := .thresholds_to_estimate,
         !!.time_to_outcome := 365.25, # Is ignored when using rms::Glm() (but not stats::glm())
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
      # /!\ Here, only X*beta is returned, not X*beta + offset. The rate will be per day
      pred <- predict(fit_i, newdata = new_data, type = "lp", se.fit = TRUE)

      log_rates[, i] <- pred$linear.predictors + log(365.25) # Add the offset manually when using rms::Glm() to convert to person-year
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
      gfr_thresholds = .thresholds_to_estimate,
      conditional_IR = conditional_IR,
      CI_lower = CI_lower,
      CI_upper = CI_upper
   )
}
