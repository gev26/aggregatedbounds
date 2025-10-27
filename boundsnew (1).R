# Clear the workspace
rm(list = ls())

# ---- Setup ------------------------------------------------------
my_path <- "~/Fall 2024/Empirical/"

req_pkgs <- c(
  "hdm","feather","tidyverse","miceadds","AER","ivmodel",
  "expm","questionr","nnet", "ivmodel"
)
to_get <- setdiff(req_pkgs, rownames(installed.packages()))
if (length(to_get)) install.packages(to_get, repos = "https://cloud.r-project.org")
invisible(lapply(req_pkgs, library, character.only = TRUE))

source(file.path(my_path, "utilsnew.R"))  # <- updated utils with orthogonal/plug-in functions

# ---- Load data --------------------------------------------------
prepared_data <- read.csv(file.path(my_path, "prepared_data.csv"))

# Filter to 12m mail survey sample (adjust the label if needed)
prepared_data_sample_12m <- dplyr::filter(prepared_data, sample_12m == "In 12m mail survey sample")

# ---- REQUIRED: Lee "selection" indicator -----------------------
# 1 if outcome observed, 0 otherwise (e.g., completed survey)
selection_var <- "ohp_all_ever_survey"   # <-- set if different
stopifnot(selection_var %in% names(prepared_data_sample_12m))

# Weights
weight_name <- "weight_12m"
stopifnot(weight_name %in% names(prepared_data_sample_12m))

# ---- Outcomes (binary) -----------------------------------------
binary_outcome_names <- c(
  "rx_any_12m","doc_any_12m","er_any_12m","hosp_any_12m",
  "health_genflip_bin_12m","health_notpoor_12m","health_chgflip_bin_12m","nodep_screen_12m"
)

binary_outcome_names <- c(
  "rx_any_12m","doc_any_12m"
)


# ---- Covariates -------------------------------------------------
pre_vars <- grep("^pre", colnames(prepared_data_sample_12m), value = TRUE)
additional_vars <- c("female_list","english_list","zip_msa")
stratification_controls <- grep("ddd", colnames(prepared_data_sample_12m), value = TRUE)

X_names <- unique(c(pre_vars, additional_vars, stratification_controls))
covariates <- prepared_data_sample_12m[, X_names, drop = FALSE]
covariates[is.na(covariates)] <- 0 #----------------

# Fill numeric NAs with 0 (leave factors as-is)
num_cols <- sapply(covariates, is.numeric)
if (any(num_cols)) {
  covariates[ , num_cols] <- lapply(covariates[ , num_cols, drop=FALSE],
                                    function(v) { v[is.na(v)] <- 0; v })
}
# Optional: keep factor NAs explicit to avoid row drops
# fac_cols <- sapply(covariates, is.factor)
# for (cn in names(covariates)[fac_cols]) covariates[[cn]] <- addNA(covariates[[cn]])

# ---- Build base lee-data frame ---------------------------------
base_df <- dplyr::tibble(
  treat     = prepared_data_sample_12m[["treatment"]],
  selection = prepared_data_sample_12m[[selection_var]],
  weights   = prepared_data_sample_12m[[weight_name]]
) %>% cbind(covariates)

# Selection formula used repeatedly
form_selection <- as.formula(
  paste0("selection ~ treat * (", paste0(colnames(covariates), collapse = "+"), ")")
)

# Results containers
initialize_matrices <- function(n_outcomes) {
  list(
    baseline          = matrix(NA_real_, n_outcomes, 1),
    share             = matrix(NA_real_, n_outcomes, 1),
    basic_bounds      = matrix(NA_real_, n_outcomes, 2),
    basic_uncond      = matrix(NA_real_, n_outcomes, 2),
    basic_CR          = matrix(NA_real_, n_outcomes, 2),
    basic_uncond_CR   = matrix(NA_real_, n_outcomes, 2),
    baseline_ci       = matrix(NA_real_, n_outcomes, 2),
    sharp_bounds      = matrix(NA_real_, n_outcomes, 2),
    sharp_CR          = matrix(NA_real_, n_outcomes, 2),
    p_values          = matrix(NA_real_, n_outcomes, 2)
  )
}
results_matrices <- initialize_matrices(length(binary_outcome_names))


prepared_data$rx_any_12m

# ---- Main loop (BINARY) ----------------------------------------
for (j in seq_along(binary_outcome_names)) {
  outcome_name <- binary_outcome_names[j]
  cat("\n>>> Working on outcome:", outcome_name, "\n")
  
  # Build leedata for this outcome
  if (!outcome_name %in% names(prepared_data_sample_12m)) {
    warning("Outcome not found: ", outcome_name); next
  }
  this_outcome <- prepared_data_sample_12m[[outcome_name]]
  # Coerce to 0/1 and guard NA
  if (is.factor(this_outcome)) this_outcome <- as.numeric(as.character(this_outcome))
  if (is.logical(this_outcome)) this_outcome <- as.numeric(this_outcome)
  this_outcome[is.na(this_outcome)] <- 0
  leedata <- cbind(base_df, outcome = this_outcome)
  

  # ---- Baseline estimation (run on leedata with cleaned outcome) ----
  # make sure leedata has 'treat', 'weights', and the stratification controls
  rhs_controls <- if (length(stratification_controls) > 0) {
    paste0("+", paste(stratification_controls, collapse = "+"))
  } else {
    ""
  }
  
  # response is the cleaned 'outcome' column inside leedata
  form_baseline <- as.formula(paste0("outcome ~ treat", rhs_controls))
  
  res <- lm(
    formula = form_baseline,
    data    = leedata,
    weights = as.numeric(leedata$weights)  # or leedata[[weight_name]] if that's different
  )
  
  # extract by name to avoid position mistakes
  coefs <- summary(res)$coefficients
  if (!("treat" %in% rownames(coefs))) {
    warning("Treatment column 'treat' not found in model matrix for ", outcome_name)
  } else {
    beta_hat <- coefs["treat", "Estimate"]
    se_hat   <- coefs["treat", "Std. Error"]
    z        <- qnorm(c(0.025, 0.975))
    
    results_matrices$baseline[j]        <- beta_hat
    results_matrices$baseline_ci[j, ]   <- beta_hat + z * se_hat
  }
  
  
  
  
  
  # Selection model & s.hat (WEIGHTED)
  glm_fit <- estimate_selection(
    leedata = leedata,
    form = form_selection,
    selection_function_name = "glm",
    myweights = leedata$weights
  )
  s.hat <- as.data.frame(predict_selection(glm_fit, leedata = leedata))
  results_matrices$share[j] <- mean(s.hat$s.0.hat > s.hat$s.1.hat)
  
  # -------- Basic Lee bounds (plug-in) -----------
  basic_res <- basic_lee_bound_binary(leedata, ortho_d = FALSE)
  results_matrices$basic_bounds[j, ] <- GetBounds(basic_res)
  
  # Bootstrap CR for plug-in bounds
  ATE_bb_basic <- main_bb_unified(
    leedata = leedata,
    N_rep   = 500,
    mode    = "binary",
    ortho   = FALSE,
    flag_helps = TRUE
  )
  results_matrices$basic_CR[j, ] <- compute_confidence_region(
    ATE_bb_basic,
    ATE_est = results_matrices$basic_bounds[j, ]
  )
  
  # -------- Orthogonal Lee bounds (DML) ---------
  # Outcome model for orthogonal path: fit on D=1 & S=1
  dsub <- (leedata$treat == 1 & leedata$selection == 1)
  if (sum(dsub) < 5) {
    warning("Too few treated&selected to fit outcome model; skipping orthogonal for ", outcome_name)
    next
  }
  form_outcome <- as.formula(paste0("outcome ~ ", paste0(colnames(covariates), collapse = "+")))
  glm_y_fit <- glm(form_outcome, data = leedata[dsub, , drop = FALSE], family = binomial(),
                   weights = leedata$weights[dsub])
  y.hat <- predict(glm_y_fit, newdata = leedata, type = "response")
  
  # choose anchor by lower average selection rate (control-anchored if s0 <= s1)
  auto_flag <- weighted.mean(s.hat$s.0.hat, w = leedata$weights) <=
    weighted.mean(s.hat$s.1.hat, w = leedata$weights)
  
  ortho_res <- ortho_leebounds(
    leedata     = leedata,
    s.hat       = s.hat,
    y.hat       = y.hat,
    flag_binary = TRUE,
    flag_helps  = auto_flag,   # <-- boolean
    ortho_d     = TRUE
  )
  
  results_matrices$sharp_bounds[j, ] <- c(ortho_res$lower_bound, ortho_res$upper_bound)
  
  # Bootstrap CR for orthogonal bounds
  ATE_bb_ortho <- main_bb_unified(
    leedata = leedata,
    N_rep   = 500,
    mode    = "binary",
    ortho   = TRUE,
    flag_helps = auto_flag
  )
  results_matrices$sharp_CR[j, ] <- compute_confidence_region(
    ATE_bb_ortho,
    ATE_est = results_matrices$sharp_bounds[j, ]
  )
  
  # Diagnostics
  cat("share(s0>s1)=", round(results_matrices$share[j],3),
      "  basic=[", paste(round(results_matrices$basic_bounds[j,],3), collapse=", "), "]",
      "  ortho=[", paste(round(results_matrices$sharp_bounds[j,],3), collapse=", "), "]\n")
}

# ---- Save binary results ---------------------------------------
results <- data.frame(
  name                 = binary_outcome_names,
  baseline             = results_matrices$baseline,
  baseline_lower       = results_matrices$baseline_ci[, 1],
  baseline_upper       = results_matrices$baseline_ci[, 2],
  basic_lower_bound    = results_matrices$basic_bounds[, 1],
  basic_upper_bound    = results_matrices$basic_bounds[, 2],
  basic_lower_bound_CR = results_matrices$basic_CR[, 1],
  basic_upper_bound_CR = results_matrices$basic_CR[, 2],
  share                = results_matrices$share,
  sharp_lower_bound    = results_matrices$sharp_bounds[, 1],
  sharp_upper_bound    = results_matrices$sharp_bounds[, 2]
)
write.csv(results, file.path(my_path, "binary.csv"), row.names = FALSE)

# =================================================================
# ===============  DISCRETE OUTCOME(S): auto-levels  ===============
# =================================================================

discrete_outcomes <- c("doc_num_mod_12m", "rx_num_mod_12m")  # <- set as needed
disc_results <- list()

for (outcome_name in discrete_outcomes) {
  cat("\n>>> DISCRETE outcome:", outcome_name, "\n")
  if (!outcome_name %in% names(prepared_data_sample_12m)) {
    warning("Outcome not found: ", outcome_name); next
  }
  
  # Build leedata for this outcome
  y_raw <- prepared_data_sample_12m[[outcome_name]]
  
  # Infer support (Y_levels) from the column safely
  if (is.factor(y_raw)) {
    Y_levels <- levels(y_raw)
    Y_levels_num <- suppressWarnings(as.numeric(Y_levels))
    if (!any(is.na(Y_levels_num))) Y_levels <- Y_levels_num
  } else {
    vals_num <- suppressWarnings(as.numeric(y_raw))
    Y_levels <- sort(unique(vals_num[!is.na(vals_num)]))
  }
  
  # Coerce outcome consistent with Y_levels
  if (is.factor(y_raw)) {
    if (is.numeric(Y_levels)) {
      y_val <- as.numeric(as.character(y_raw))
    } else {
      labs <- levels(y_raw)
      map  <- setNames(seq_along(labs), labs)
      y_val <- as.numeric(map[as.character(y_raw)])
      Y_levels <- seq_along(labs)
    }
  } else {
    y_val <- suppressWarnings(as.numeric(y_raw))
  }
  if (anyNA(y_val)) y_val[is.na(y_val)] <- min(Y_levels, na.rm = TRUE)
  
  leedata_disc <- cbind(base_df, outcome = y_val)
  auto_flag <- weighted.mean(s.hat$s.0.hat, w = leedata_disc$weights) <=
    weighted.mean(s.hat$s.1.hat, w = leedata_disc$weights)
  
  # Selection model → s.hat (WEIGHTED)
  # Selection model → s.hat (WEIGHTED)
  sel_fit <- estimate_selection(
    leedata = leedata_disc,
    form = form_selection,
    selection_function_name = "glm",
    myweights = leedata_disc$weights          # <- FIXED
  )
  s.hat <- as.data.frame(predict_selection(sel_fit, leedata_disc))
  
  # Multinomial on D=1 & S=1 to get pi_hat (n x K) aligned to Y_levels
  dsub <- (leedata_disc$treat == 1 & leedata_disc$selection == 1)
  if (sum(dsub) < max(5, length(unique(leedata_disc$outcome[dsub])))) {
    warning("Too few treated&selected to fit multinomial for ", outcome_name); next
  }
  sub_df <- leedata_disc[dsub, , drop = FALSE]
  
  # Factor with levels in Y_levels order
  if (is.numeric(Y_levels)) {
    sub_df$y_fac <- factor(sub_df$outcome, levels = Y_levels, labels = Y_levels)
  } else {
    sub_df$y_fac <- factor(sub_df$outcome, levels = Y_levels)
  }
  
  form_outcome_mult <- as.formula(paste0("y_fac ~ ", paste0(colnames(covariates), collapse = "+")))
  mfit <- nnet::multinom(form_outcome_mult, data = sub_df,
                         weights = leedata_disc$weights[dsub], trace = FALSE)
  
  # Predict class probabilities; pad missing classes with zeros and reorder
  pi_hat <- predict(mfit, newdata = leedata_disc, type = "probs")
  if (is.vector(pi_hat)) {
    p1 <- as.numeric(pi_hat)
    pi_hat <- cbind(p1, 1 - p1)
    colnames(pi_hat) <- as.character(Y_levels)  # assumes binary support
  } else {
    want <- as.character(Y_levels)
    have <- colnames(pi_hat)
    missing <- setdiff(want, have)
    if (length(missing) > 0) {
      pi_hat <- cbind(pi_hat,
                      matrix(0, nrow = nrow(pi_hat), ncol = length(missing),
                             dimnames = list(NULL, missing)))
    }
    pi_hat <- pi_hat[, want, drop = FALSE]
  }
  pi_hat <- as.matrix(pi_hat)
  
  # choose anchor once (boolean) for orthogonal call
  auto_flag <- weighted.mean(s.hat$s.0.hat, w = leedata_disc$weights) <=
    weighted.mean(s.hat$s.1.hat, w = leedata_disc$weights)
  
  # -------- Plug-in discrete bounds --------
  plug_res <- lee_bounds_discrete(
    leedata    = leedata_disc,
    s.hat      = s.hat,
    pi_hat     = pi_hat,
    Y_levels   = as.numeric(Y_levels),
    flag_helps = TRUE,
    ortho_d    = FALSE
  )
  
  # -------- Orthogonal discrete bounds --------
  ortho_res <- lee_bounds_discrete(
    leedata    = leedata_disc,
    s.hat      = s.hat,
    pi_hat     = pi_hat,
    Y_levels   = as.numeric(Y_levels),
    flag_helps = auto_flag,
    ortho_d    = TRUE
  )
  
  # Bootstraps (unchanged)
  # -------- Bootstraps --------
  ATE_bb_plug  <- main_bb_unified(
    leedata = leedata_disc,
    N_rep   = 500,
    mode    = "discrete",
    ortho   = FALSE,
    flag_helps = TRUE,
    Y_levels  = as.numeric(Y_levels)   
  )
  
  ATE_bb_ortho <- main_bb_unified(
    leedata = leedata_disc,
    N_rep   = 500,
    mode    = "discrete",
    ortho   = TRUE,
    flag_helps = auto_flag,
    Y_levels  = as.numeric(Y_levels)   
  )
  
  
  
  disc_results[[outcome_name]] <- list(
    basic_bounds = c(plug_res$lower_bound,  plug_res$upper_bound),
    basic_CR     = compute_confidence_region(ATE_bb_plug,
                                             ATE_est = c(plug_res$lower_bound,  plug_res$upper_bound)),
    sharp_bounds = c(ortho_res$lower_bound, ortho_res$upper_bound),
    sharp_CR     = compute_confidence_region(ATE_bb_ortho,
                                             ATE_est = c(ortho_res$lower_bound, ortho_res$upper_bound))
  )
  
  cat("   plug-in   =[", paste(round(disc_results[[outcome_name]]$basic_bounds,3), collapse=", "), "]",
      "  orthogonal=[", paste(round(disc_results[[outcome_name]]$sharp_bounds,3), collapse=", "), "]\n")
}

# Save discrete results (if any)
if (length(disc_results)) {
  res_disc_df <- do.call(rbind, lapply(names(disc_results), function(nm) {
    x <- disc_results[[nm]]
    data.frame(
      name = nm,
      basic_lower_bound      = x$basic_bounds[1],
      basic_upper_bound      = x$basic_bounds[2],
      basic_lower_bound_CR   = x$basic_CR[1],
      basic_upper_bound_CR   = x$basic_CR[2],
      sharp_lower_bound      = x$sharp_bounds[1],
      sharp_upper_bound      = x$sharp_bounds[2],
      sharp_lower_bound_CR   = x$sharp_CR[1],
      sharp_upper_bound_CR   = x$sharp_CR[2]
    )
  }))
  write.csv(res_disc_df, file.path(my_path, "discrete.csv"), row.names = FALSE)
}
