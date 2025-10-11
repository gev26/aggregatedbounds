### ============================
### orthogonal_bounds.R  (final)
### ============================

rm(list = ls())
set.seed(42)

my_path <- "~/Fall 2024/Empirical/"

suppressPackageStartupMessages({
  library(hdm)
  library(feather)
  library(tidyverse)
  library(radiant.data)
  library(miceadds)
  library(AER)
  library(ivmodel)
  library(expm)
  library(questionr)
})

source(paste0(my_path, "new_utils.R"))

# Load data
prepared_data <- read.csv(paste0(my_path, "prepared_data.csv"))
prepared_data_sample_12m <- dplyr::filter(prepared_data, sample_12m == "In 12m mail survey sample")
prepared_data_sample_12m$weight_12m <- as.numeric(prepared_data_sample_12m$weight_12m)

# Outcomes (binary here)
outcome_names <- c(
  "rx_any_12m", "doc_any_12m", "er_any_12m", "hosp_any_12m",
  "health_genflip_bin_12m", "health_notpoor_12m", "health_chgflip_bin_12m", "nodep_screen_12m"
)

# Covariates
pre_vars                <- grep("pre", colnames(prepared_data), value = TRUE)
additional_vars         <- c("female_list", "english_list", "zip_msa")
stratification_controls <- grep("ddd", colnames(prepared_data), value = TRUE)
covariates <- prepared_data_sample_12m[, c(pre_vars, additional_vars, stratification_controls)]
covariates[is.na(covariates)] <- 0

# Formulas (reused)
form_selection <- as.formula(paste0("selection ~ treat * (", paste0(colnames(covariates), collapse = "+"), ")"))
form_outcome   <- as.formula(paste0("outcome ~ (",   paste0(colnames(covariates), collapse = "+"), ")"))

initialize_matrices <- function(n_outcomes) {
  list(
    baseline          = matrix(0, n_outcomes, 1),
    share             = matrix(0, n_outcomes, 1),
    basic_bounds      = matrix(0, n_outcomes, 2),
    basic_uncond      = matrix(0, n_outcomes, 2),
    basic_CR          = matrix(0, n_outcomes, 2),
    basic_uncond_CR   = matrix(0, n_outcomes, 2),
    baseline_ci       = matrix(0, n_outcomes, 2),
    sharp_bounds      = matrix(0, n_outcomes, 2),
    sharp_CR          = matrix(0, n_outcomes, 2),
    p_values          = matrix(0, n_outcomes, 2)
  )
}

results_matrices <- initialize_matrices(length(outcome_names))
N_rep_basic <- 200
N_rep_ortho <- 200

for (j in seq_along(outcome_names)) {
  outcome_name <- outcome_names[j]
  
  ## Baseline ITT with strat controls
  form_baseline <- as.formula(paste0(outcome_name, " ~ treatment + ", paste0(stratification_controls, collapse = "+")))
  res <- lm(
    formula = form_baseline,
    data    = prepared_data_sample_12m,
    weights = prepared_data_sample_12m$weight_12m
  )
  results_matrices$baseline[j]    <- summary(res)$coefficients[2, 1]
  ci_vec <- qnorm(c(0.025, 0.975)) * summary(res)$coefficients[2, 2]
  results_matrices$baseline_ci[j,] <- summary(res)$coefficients[2, 1] + ci_vec
  
  ## Build Lee data
  y_raw <- prepared_data_sample_12m[[outcome_name]]
  leedata <- data.frame(
    treat     = prepared_data_sample_12m$treatment,
    selection = as.numeric(!is.na(y_raw)),
    outcome   = ifelse(is.na(y_raw), 0, y_raw),
    weights   = prepared_data_sample_12m$weight_12m,
    covariates
  )
  leedata$weights <- as.numeric(leedata$weights)
  
  ## Detect if outcome is binary 0/1 (if you later add non-binary discrete, this will switch automatically)
  is_binary <- all(leedata$outcome %in% c(0,1))
  
  ## Unconditional basic Lee bounds + CR
  if (is_binary) {
    res_basic <- basic_lee_bound_binary(leedata)
    fun_basic <- basic_lee_bound_binary
  } else {
    res_basic <- basic_lee_bound_discrete(leedata)
    fun_basic <- basic_lee_bound_discrete
  }
  results_matrices$basic_uncond[j, ] <- GetBounds(res_basic)
  
  ATE_bb_uncond <- main_bb(
    function_name = fun_basic,
    sample_size   = nrow(leedata),
    N_rep         = N_rep_basic,
    leedata       = leedata
  )
  results_matrices$basic_uncond_CR[j, ] <- compute_confidence_region(ATE_bb_uncond, ATE_est = results_matrices$basic_uncond[j, ])
  
  ## Intersection test (for binary outcomes this is meaningful; we still report for completeness)
  results_matrices$p_values[j, ] <- intersection_test(leedata, leedata$weights)
  
  ## Selection model and share (weighted)
  sel_fit <- estimate_selection(leedata = leedata, form = form_selection, selection_function_name = "glm")
  sel_hat <- predict_selection(sel_fit, leedata = leedata)
  results_matrices$share[j] <- weighted.mean(sel_hat$s.0.hat > sel_hat$s.1.hat, w = leedata$weights)
  
  ## Sharp subset: s0 > s1
  keep_idx <- (leedata$weights != 0) & (sel_hat$s.0.hat > sel_hat$s.1.hat)
  myleedata <- leedata[keep_idx, , drop = FALSE]
  
  ## Basic sharp bounds + CR (non-orthogonal)
  if (is_binary) {
    sharp_basic <- basic_lee_bound_binary(myleedata)
    fun_sharp   <- basic_lee_bound_binary
  } else {
    sharp_basic <- basic_lee_bound_discrete(myleedata)
    fun_sharp   <- basic_lee_bound_discrete
  }
  results_matrices$basic_bounds[j, ] <- GetBounds(sharp_basic)
  
  ATE_bb_sharp <- main_bb(
    function_name = fun_sharp,
    sample_size   = nrow(myleedata),
    N_rep         = N_rep_basic,
    leedata       = myleedata
  )
  results_matrices$basic_CR[j, ] <- compute_confidence_region(ATE_bb_sharp, ATE_est = results_matrices$basic_bounds[j, ])
  
  ## ORTHOGONAL SHARP BOUNDS with CROSS-FITTING
  # Cross-fit selection & outcome within sharp subset
  sel_vars_cf <- setdiff(colnames(myleedata), c("treat","selection","outcome","weights"))
  sel_form_cf <- as.formula(paste0("selection ~ treat * (", paste0(sel_vars_cf, collapse = "+"), ")"))
  y_form_cf   <- as.formula(paste0("outcome ~ (", paste0(sel_vars_cf, collapse = "+"), ")"))
  
  cf_preds <- crossfit_kfold(df = myleedata, K = 5,
                             sel_form = sel_form_cf, y_form = y_form_cf, weight_col = "weights")
  s.hat.sub <- cf_preds[, c("s.0.hat", "s.1.hat")]
  y.hat.sub <- cf_preds$y.hat
  
  # Orthogonal bounds
  ortho_res <- ortho_leebounds(
    leedata = myleedata,
    s.hat   = s.hat.sub,
    y.hat   = y.hat.sub,
    flag_binary   = is_binary,
    flag_discrete = !is_binary,
    flag_helps    = FALSE,
    ortho_d       = TRUE
  )
  results_matrices$sharp_bounds[j, ] <- c(ortho_res$lower_bound, ortho_res$upper_bound)
  
  # Bootstrap CR for orthogonal sharp bounds (re-estimate nuisances inside)
  ATE_bb_ortho <- main_bb_ortho(
    sample_size = nrow(myleedata),
    N_rep       = N_rep_ortho,
    leedata     = myleedata,
    s.hat       = s.hat.sub,
    y.hat       = y.hat.sub,
    flag_binary = is_binary,
    flag_discrete = !is_binary,
    flag_helps  = FALSE
  )
  results_matrices$sharp_CR[j, ] <- compute_confidence_region(ATE_bb_ortho,
                                                              ATE_est = results_matrices$sharp_bounds[j, ])
}

# Collect results
results <- data.frame(
  name                    = outcome_names,
  baseline                = results_matrices$baseline,
  baseline_lower          = results_matrices$baseline_ci[, 1],
  baseline_upper          = results_matrices$baseline_ci[, 2],
  basic_lower_bound       = results_matrices$basic_uncond[, 1],
  basic_upper_bound       = results_matrices$basic_uncond[, 2],
  basic_lower_bound_CR    = results_matrices$basic_uncond_CR[, 1],
  basic_upper_bound_CR    = results_matrices$basic_uncond_CR[, 2],
  share                   = results_matrices$share,
  sharp_lower_bound       = results_matrices$sharp_bounds[, 1],
  sharp_upper_bound       = results_matrices$sharp_bounds[, 2],
  sharp_lower_bound_CR    = results_matrices$sharp_CR[, 1],
  sharp_upper_bound_CR    = results_matrices$sharp_CR[, 2],
  p_value_h1              = results_matrices$p_values[, 1],
  p_value_h2              = results_matrices$p_values[, 2]
)

write.csv(results, paste0(my_path, "binary.csv"), row.names = FALSE)
print("Analysis completed successfully!")
print(results)




# --- Test Proposition 8 (HML, discrete Y) ---
source("utils.R")

set.seed(2025)
N  <- 4000
X1 <- rnorm(N)
X2 <- rnorm(N)
mu_true <- 0.5
D  <- rbinom(N, 1, mu_true)

# Monotone selection S(1) >= S(0)
s0_0 <- plogis(-0.5 + 0.7*X1)     # Pr(S=1|D=0,X)
s0_1 <- plogis( 0.2 + 0.7*X1)     # Pr(S=1|D=1,X) >= s0_0
S    <- rbinom(N, 1, ifelse(D==1, s0_1, s0_0))

# Discrete outcomes on T = {0,1,2,3}; define Y|D=1,S=1,X via multinomial logits
eta1 <-  0.3 + 0.6*X2
eta2 <- -0.4 + 0.3*X2
eta3 <- -1.0 + 0.2*X2
den  <- 1 + exp(eta1) + exp(eta2) + exp(eta3)
p1   <- exp(eta1)/den; p2 <- exp(eta2)/den; p3 <- exp(eta3)/den; p0 <- 1 - p1 - p2 - p3
Y1   <- apply(cbind(p0, p1, p2, p3), 1, function(pr) sample(0:3, size = 1, prob = pr))
Y0   <- rbinom(N, 3, plogis(-0.3 + 0.5*X2))  # for D=0,S=1 (not used in PMF step)
Y    <- ifelse(D==1, Y1, Y0)
Y[S==0] <- NA_integer_

df   <- data.frame(y = Y, d = D, s = S, X1 = X1, X2 = X2)

res  <- estimate_hml_discrete(
  df, y = "y", d = "d", s = "s",
  x_vars = c("X1","X2"),
  K = 5, seed = 99,
  support_Y = 0:3,
  known_propensity = mu_true
)

cat("\n=== Proposition 8 (HML discrete) ===\n")
cat("Bounds (L,U): ", sprintf("(%.4f, %.4f)\n", res$bounds["L"], res$bounds["U"]))
cat("Sorted bounds: ", sprintf("(%.4f, %.4f)\n", res$sorted_bounds["L"], res$sorted_bounds["U"]))
cat("95%% CI for L: ", sprintf("(%.4f, %.4f)\n", res$ci_marginal$L[1], res$ci_marginal$L[2]))
cat("95%% CI for U: ", sprintf("(%.4f, %.4f)\n", res$ci_marginal$U[1], res$ci_marginal$U[2]))

