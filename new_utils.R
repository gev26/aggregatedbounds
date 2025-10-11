### ============================
### utils.R  (final)
### ============================

library(stats)

## -------------------------
## Small helpers (safe, weighted)
## -------------------------
.wmean <- function(x, weights_vec = NULL) {
  if (is.null(weights_vec)) weights_vec <- rep(1, length(x))
  sum(weights_vec * x) / sum(weights_vec)
}

## Safe binary outcome fit on controls & selected; falls back to intercept-only
## ---- helper: safe binary outcome fit on controls & selected (patched)
## ---- helper: safe binary outcome fit on controls & selected (final patch)
.safe_binom_yfit <- function(ctrl_sel, y_form, weight_col = "weights") {
  has_rows <- nrow(ctrl_sel) > 0
  has_var  <- has_rows && length(unique(ctrl_sel$outcome)) > 1
  
  if (!has_rows || !has_var) {
    # Fallback: intercept-only with weighted mean (if no rows, use 0.5)
    wv <- if (has_rows && weight_col %in% names(ctrl_sel)) as.numeric(ctrl_sel[[weight_col]]) else rep(1, max(nrow(ctrl_sel), 1))
    p_hat <- if (has_rows) sum(wv * ctrl_sel$outcome) / sum(wv) else 0.5
    return(list(kind = "intercept", p = p_hat))
  }
  
  # Ensure a concrete weights column exists INSIDE the data frame used by glm()
  if (weight_col %in% names(ctrl_sel)) {
    ctrl_sel$.__w <- as.numeric(ctrl_sel[[weight_col]])
  } else {
    ctrl_sel$.__w <- rep(1, nrow(ctrl_sel))
  }
  
  # IMPORTANT: pass weights as a bare column name so glm() finds it in `data`
  fit <- glm(
    formula = y_form,
    data    = ctrl_sel,
    weights = .__w,
    family  = binomial()
  )
  list(kind = "glm", fit = fit)
}


## Weighted one-sided trimming mean for discrete outcomes
.wtrim_mean <- function(y, weights_vec, alpha, side = c("upper","lower")) {
  side <- match.arg(side)
  stopifnot(length(y) == length(weights_vec), alpha >= 0, alpha < 1)
  if (alpha == 0 || length(y) == 0) return(sum(weights_vec * y) / sum(weights_vec))
  
  o <- order(y)  # ascending
  y <- y[o]; weights_vec <- weights_vec[o]
  cumw <- cumsum(weights_vec); totw <- tail(cumw, 1)
  
  if (side == "upper") {
    cut_w <- (1 - alpha) * totw
    keep <- cumw <= cut_w + 1e-12
    if (any(!keep)) {
      k <- which.max(!keep)
      if (k > 1 && cumw[k-1] < cut_w && cumw[k] > cut_w) {
        frac_keep <- (cut_w - cumw[k-1]) / weights_vec[k]
        weights_vec[k] <- weights_vec[k] * max(min(frac_keep, 1), 0)
        keep[k] <- TRUE
      }
    }
    yk <- y[keep]; wk <- weights_vec[keep]
  } else {
    cut_w <- alpha * totw
    keep <- cumw > cut_w - 1e-12
    if (any(keep) && any(!keep)) {
      k <- which.max(keep)
      if (k > 1 && cumw[k-1] < cut_w && cumw[k] > cut_w) {
        frac_keep <- (cumw[k] - cut_w) / weights_vec[k]
        weights_vec[k] <- weights_vec[k] * max(min(frac_keep, 1), 0)
      }
    }
    yk <- y[keep]; wk <- weights_vec[keep]
  }
  sum(wk * yk) / sum(wk)
}

## -------------------------
## Selection model (logit); optional rlassologit screen
## -------------------------
estimate_selection <- function(leedata, form = NULL, selection_function_name,
                               variables_for_selection = NULL, names_to_include = c(),
                               treat_name = "treat+", yname = "selection",
                               myweights = NULL, ...) {
  
  if (is.null(variables_for_selection)) {
    variables_for_selection <- setdiff(
      colnames(leedata),
      c("treat","selection","outcome","X.Intercept.","(Intercept)","weights","prop1","prop0")
    )
  } else {
    variables_for_selection <- unique(setdiff(
      variables_for_selection,
      c("outcome","X.Intercept.","(Intercept)","weights","prop1","prop0")
    ))
  }
  
  if (is.null(form)) {
    form <- as.formula(paste0("selection ~ (treat) * (", paste0(variables_for_selection, collapse = "+"), ")"))
  }
  
  if (is.null(myweights)) {
    myweights <- rep(1, nrow(leedata))
  }
  
  if (selection_function_name == "rlassologit") {
    glm.fit <- rlassologit(form, leedata[, c("treat", "selection", variables_for_selection)], family = "binomial", ...)
    non_zero_coefs <- glm.fit$coefficients[glm.fit$coefficients != 0]
    selected_names <- setdiff(names(non_zero_coefs), c("(Intercept)", "treat"))
    selected_names <- unique(c(selected_names, names_to_include))
    
    if (length(selected_names) > 0) {
      form <- as.formula(paste0(yname, "~", treat_name, paste0(selected_names, collapse = "+")))
    } else {
      form <- as.formula(paste0(yname, "~", treat_name))
    }
  }
  
  leedata$myweights <- myweights
  glm.postlasso <- glm(
    form,
    data   = leedata[, c("treat", "selection", variables_for_selection)],
    family = binomial()
  )
  glm.postlasso
}

predict_selection <- function(fit, leedata, ...) {
  leedata_0treat <- leedata; leedata_0treat$treat <- 0
  leedata_1treat <- leedata; leedata_1treat$treat <- 1
  s.0.hat <- predict(fit, leedata_0treat, type = "response")
  s.1.hat <- predict(fit, leedata_1treat, type = "response")
  data.frame(s.0.hat = s.0.hat, s.1.hat = s.1.hat)
}

## -------------------------
## Basic Lee bounds: binary outcome
## -------------------------
basic_lee_bound_binary <- function(leedata, inds = NULL, ...) {
  if (!is.null(inds)) leedata <- leedata[inds, ]
  
  d  <- leedata$treat
  s  <- leedata$selection
  sy <- leedata$outcome
  weights_vec <- if ("weights" %in% colnames(leedata)) as.numeric(leedata$weights) else rep(1, length(d))
  
  prop_control_nonmissing <- weighted.mean(s[d == 0] == 1, weights_vec[d == 0])
  prop_treat_nonmissing   <- weighted.mean(s[d == 1] == 1, weights_vec[d == 1])
  p0 <- prop_control_nonmissing / prop_treat_nonmissing
  
  if (p0 <= 1) {
    p0_star <- p0
    pY <- weighted.mean(sy[d == 1 & s == 1] == 1, weights_vec[d == 1 & s == 1])
    lower_bound <- max(1 - (1 - pY) / p0, 0) - weighted.mean(sy[d == 0 & s == 1], weights_vec[d == 0 & s == 1])
    upper_bound <- min(pY / p0 - 1, 0) + 1   - weighted.mean(sy[d == 0 & s == 1], weights_vec[d == 0 & s == 1])
  } else {
    p0_star <- 1 / p0
    pY <- weighted.mean(sy[d == 0 & s == 1] == 1, weights_vec[d == 0 & s == 1])
    lower_bound <- weighted.mean(sy[d == 1 & s == 1], weights_vec[d == 1 & s == 1]) - (min(pY / p0_star - 1, 0) + 1)
    upper_bound <- weighted.mean(sy[d == 1 & s == 1], weights_vec[d == 1 & s == 1]) - (max(1 - (1 - pY) / p0_star, 0))
  }
  
  list(lower_bound = lower_bound,
       upper_bound = upper_bound,
       s0 = prop_control_nonmissing,
       s1 = prop_treat_nonmissing,
       p0_star = p0_star, pY = pY)
}

## -------------------------
## Denominator (with orthogonal correction)
## -------------------------
denominator <- function(s.hat, leedata, weights = NULL, flag_helps = FALSE, inds = NULL, ortho_d = FALSE, ...) {
  if (!is.null(inds)) { s.hat <- s.hat[inds, ]; leedata <- leedata[inds, ] }
  if (is.null(weights)) {
    weights <- if ("weights" %in% colnames(leedata)) as.numeric(leedata$weights) else rep(1, nrow(leedata))
  }
  
  d <- leedata$treat
  s <- leedata$selection
  
  if (!ortho_d) {
    return(if (flag_helps) weighted.mean(s[d == 0] == 1, weights[d == 0])
           else            weighted.mean(s[d == 1] == 1, weights[d == 1]))
  }
  
  if (flag_helps) {
    weighted.mean(s.hat$s.0.hat + (1 - d) * (s - s.hat$s.0.hat), weights)
  } else {
    weighted.mean(s.hat$s.1.hat + d       * (s - s.hat$s.1.hat), weights)
  }
}

## -------------------------
## Numerator (binary): orthogonal envelope score
## -------------------------
numerator_binary <- function(s.hat, y.hat, weights = NULL, leedata = NULL,
                             flag_helps = FALSE, inds = NULL, ortho_d = FALSE, ...) {
  if (!is.null(inds)) { s.hat <- s.hat[inds, ]; y.hat <- y.hat[inds]; leedata <- leedata[inds, ] }
  
  d  <- leedata$treat
  s  <- leedata$selection
  sy <- leedata$outcome
  if (is.null(weights)) {
    weights <- if ("weights" %in% colnames(leedata)) as.numeric(leedata$weights) else rep(1, nrow(leedata))
  }
  
  if (!ortho_d) {
    if (flag_helps) {
      nl <- weighted.mean(pmax(s.hat$s.0.hat - s.hat$s.1.hat * (1 - y.hat), 0), weights)
      nu <- weighted.mean(pmin(s.hat$s.1.hat * y.hat - s.hat$s.0.hat, 0) + s.hat$s.0.hat, weights)
    } else {
      nl <- weighted.mean(pmax(s.hat$s.1.hat - s.hat$s.0.hat * (1 - y.hat), 0), weights)
      nu <- weighted.mean(pmin(s.hat$s.0.hat * y.hat - s.hat$s.1.hat, 0) + s.hat$s.1.hat, weights)
    }
    return(list(nl = nl, nu = nu))
  }
  
  if (flag_helps) {
    g1 <- s.hat$s.0.hat - s.hat$s.1.hat * (1 - y.hat)
    g2 <- s.hat$s.1.hat * y.hat - s.hat$s.0.hat
    I1_l <- as.numeric(g1 > 0)
    I2_u <- as.numeric(g2 < 0)
    base_nl <- pmax(g1, 0)
    base_nu <- pmin(g2, 0) + s.hat$s.0.hat
    
    ortho_s0_nl <- (1 - d) * (s - s.hat$s.0.hat) * I1_l
    ortho_s0_nu <- (1 - d) * (s - s.hat$s.0.hat) * (1 - I2_u)
    
    ortho_s1_nl <- d * (s - s.hat$s.1.hat) * I1_l * (-(1 - y.hat))
    ortho_s1_nu <- d * (s - s.hat$s.1.hat) * I2_u * y.hat
    
    ortho_y_nl  <- (1 - d) * s * (sy - y.hat) * I1_l * s.hat$s.1.hat
    ortho_y_nu  <- (1 - d) * s * (sy - y.hat) * I2_u * (-s.hat$s.1.hat)
    
    nl <- weighted.mean(base_nl + ortho_s0_nl + ortho_s1_nl + ortho_y_nl, weights)
    nu <- weighted.mean(base_nu + ortho_s0_nu + ortho_s1_nu + ortho_y_nu, weights)
  } else {
    g1 <- s.hat$s.1.hat - s.hat$s.0.hat * (1 - y.hat)
    g2 <- s.hat$s.0.hat * y.hat - s.hat$s.1.hat
    I1_l <- as.numeric(g1 > 0)
    I2_u <- as.numeric(g2 < 0)
    base_nl <- pmax(g1, 0)
    base_nu <- pmin(g2, 0) + s.hat$s.1.hat
    
    ortho_s1_nl <- d * (s - s.hat$s.1.hat) * I1_l
    ortho_s1_nu <- d * (s - s.hat$s.1.hat) * (1 - I2_u)
    
    ortho_s0_nl <- (1 - d) * (s - s.hat$s.0.hat) * I1_l * (-(1 - y.hat))
    ortho_s0_nu <- (1 - d) * (s - s.hat$s.0.hat) * I2_u * y.hat
    
    ortho_y_nl  <- (1 - d) * s * (sy - y.hat) * I1_l * s.hat$s.0.hat
    ortho_y_nu  <- (1 - d) * s * (sy - y.hat) * I2_u * (-s.hat$s.0.hat)
    
    nl <- weighted.mean(base_nl + ortho_s1_nl + ortho_s0_nl + ortho_y_nl, weights)
    nu <- weighted.mean(base_nu + ortho_s1_nu + ortho_s0_nu + ortho_y_nu, weights)
  }
  list(nl = nl, nu = nu)
}

## -------------------------
## Basic Lee bounds: discrete outcomes (finite support)
## -------------------------
basic_lee_bound_discrete <- function(leedata, inds = NULL, ...) {
  if (!is.null(inds)) leedata <- leedata[inds, , drop = FALSE]
  
  d  <- leedata$treat
  s  <- leedata$selection
  y  <- leedata$outcome
  weights_vec <- if ("weights" %in% colnames(leedata)) as.numeric(leedata$weights) else rep(1, length(d))
  
  s1 <- weighted.mean(s[d == 1] == 1, weights_vec[d == 1])
  s0 <- weighted.mean(s[d == 0] == 1, weights_vec[d == 0])
  if (s1 <= 0 || s0 <= 0) {
    return(list(lower_bound = NA_real_, upper_bound = NA_real_,
                s0 = s0, s1 = s1, p0_star = NA_real_))
  }
  
  p0 <- s0 / s1
  E0 <- { idx <- (d == 0 & s == 1); sum(weights_vec[idx] * y[idx]) / sum(weights_vec[idx]) }
  E1 <- { idx <- (d == 1 & s == 1); sum(weights_vec[idx] * y[idx]) / sum(weights_vec[idx]) }
  
  if (p0 <= 1) {
    alpha <- 1 - p0
    idx_t <- (d == 1 & s == 1)
    mu1_L <- .wtrim_mean(y[idx_t], weights_vec[idx_t], alpha = alpha, side = "upper")
    mu1_U <- .wtrim_mean(y[idx_t], weights_vec[idx_t], alpha = alpha, side = "lower")
    lower_bound <- mu1_L - E0
    upper_bound <- mu1_U - E0
    p0_star <- p0
  } else {
    alpha_c <- 1 - 1 / p0
    idx_c <- (d == 0 & s == 1)
    mu0_max <- .wtrim_mean(y[idx_c], weights_vec[idx_c], alpha = alpha_c, side = "lower")
    mu0_min <- .wtrim_mean(y[idx_c], weights_vec[idx_c], alpha = alpha_c, side = "upper")
    lower_bound <- E1 - mu0_max
    upper_bound <- E1 - mu0_min
    p0_star <- 1 / p0
  }
  
  list(lower_bound = lower_bound, upper_bound = upper_bound, s0 = s0, s1 = s1, p0_star = p0_star)
}

## -------------------------
## Numerator (discrete): selection-orthogonal mapping
## -------------------------
numerator_discrete <- function(s.hat, y.hat, leedata = NULL,
                               flag_helps = FALSE, inds = NULL,
                               ortho_d = FALSE, ...) {
  if (!is.null(inds)) leedata <- leedata[inds, , drop = FALSE]
  
  d  <- leedata$treat
  s  <- leedata$selection
  y  <- leedata$outcome
  weights_vec <- if ("weights" %in% colnames(leedata)) as.numeric(leedata$weights) else rep(1, nrow(leedata))
  
  bds <- basic_lee_bound_discrete(leedata)
  E0 <- { idx <- (d == 0 & s == 1); if (sum(idx)) sum(weights_vec[idx] * y[idx]) / sum(weights_vec[idx]) else NA_real_ }
  E1 <- { idx <- (d == 1 & s == 1); if (sum(idx)) sum(weights_vec[idx] * y[idx]) / sum(weights_vec[idx]) else NA_real_ }
  
  if (flag_helps) {
    denom <- weighted.mean(s[d == 0] == 1, weights_vec[d == 0])
    nl <- denom * (bds$lower_bound + E0)
    nu <- denom * (bds$upper_bound + E0)
  } else {
    denom <- weighted.mean(s[d == 1] == 1, weights_vec[d == 1])
    nu <- denom * (E1 - bds$lower_bound)
    nl <- denom * (E1 - bds$upper_bound)
  }
  list(nl = nl, nu = nu)
}

## -------------------------
## Unified orthogonal bounds wrapper
## -------------------------
ortho_leebounds <- function(leedata, s.hat, y.hat,
                            flag_binary = TRUE, flag_helps = FALSE,
                            flag_discrete = FALSE, ortho_d = TRUE, weights = NULL, ...) {
  d  <- leedata$treat
  s  <- leedata$selection
  sy <- leedata$outcome
  
  if (is.null(weights)) {
    weights <- if ("weights" %in% colnames(leedata)) as.numeric(leedata$weights) else rep(1, length(d))
  }
  
  if (flag_binary) {
    num <- numerator_binary(
      s.hat = s.hat, y.hat = y.hat, leedata = leedata,
      flag_helps = flag_helps, weights = weights, ortho_d = ortho_d, ...
    )
  } else if (flag_discrete) {
    num <- numerator_discrete(
      s.hat = s.hat, y.hat = y.hat, leedata = leedata,
      flag_helps = flag_helps, weights = weights, ortho_d = ortho_d, ...
    )
  } else {
    stop("Specify flag_binary or flag_discrete.")
  }
  
  denom <- denominator(
    s.hat = s.hat, leedata = leedata, flag_helps = flag_helps,
    weights = weights, ortho_d = ortho_d, ...
  )
  
  if (flag_helps) {
    lower_bound <- num$nl / denom - weighted.mean(sy[d == 0 & s == 1], weights[d == 0 & s == 1])
    upper_bound <- num$nu / denom - weighted.mean(sy[d == 0 & s == 1], weights[d == 0 & s == 1])
  } else {
    lower_bound <- weighted.mean(sy[d == 1 & s == 1], weights[d == 1 & s == 1]) - num$nu / denom
    upper_bound <- weighted.mean(sy[d == 1 & s == 1], weights[d == 1 & s == 1]) - num$nl / denom
  }
  
  list(lower_bound = lower_bound, upper_bound = upper_bound)
}

## -------------------------
## Bootstrap helpers
## -------------------------
main_bb <- function(function_name, sample_size, N_rep = 200, ...) {
  ATE_bb <- matrix(NA_real_, N_rep, 2)
  for (b in 1:N_rep) {
    set.seed(b)
    inds <- sample(1:sample_size, sample_size, replace = TRUE)
    resultb <- try(function_name(inds, ...), silent = TRUE)
    if (!inherits(resultb, "try-error")) {
      ATE_bb[b, ] <- c(resultb$lower_bound, resultb$upper_bound)
    }
  }
  ATE_bb
}

compute_confidence_region <- function(ATE_boot, ATE_est, ci_alpha = 0.05, tol = 1e-5) {
  Omega.hat <- matrix(0, 2, 2)
  if (sum(is.na(ATE_boot)) + sum(is.na(ATE_est)) > 0) {
    return(c(lower_bound = NA, upper_bound = NA))
  }
  ATE_boot_centered <- cbind(ATE_boot[,1] - ATE_est[1], ATE_boot[,2] - ATE_est[2])
  
  Omega.hat[1,1] <- var(ATE_boot_centered[,1])
  Omega.hat[2,2] <- var(ATE_boot_centered[,2])
  Omega.hat[1,2] <- cov(ATE_boot_centered[,1], ATE_boot_centered[,2])
  Omega.hat[2,1] <- Omega.hat[1,2]
  
  sm <- sqrtm(Omega.hat)
  if (max(abs(Im(sm))) > tol) stop("Non-trivial imaginary part!")
  sv <- Re(sm %*% c(-qnorm(sqrt(1 - ci_alpha)), qnorm(sqrt(1 - ci_alpha))))
  
  c(lower_bound = ATE_est[1] + sv[1],
    upper_bound = ATE_est[2] + sv[2])
}

GetBounds <- function(x) c(x$lower_bound, x$upper_bound)

intersection_test <- function(data, weights) {
  weights <- as.numeric(weights)
  weights <- weights / sum(weights)
  
  prop_treat   <- sum(weights[data$treat == 1] * data$selection[data$treat == 1]) / sum(weights[data$treat == 1])
  prop_control <- sum(weights[data$treat == 0] * data$selection[data$treat == 0]) / sum(weights[data$treat == 0])
  prop_outcome <- sum(weights[data$treat == 0 & data$selection == 1] * data$outcome[data$treat == 0 & data$selection == 1]) /
    sum(weights[data$treat == 0 & data$selection == 1])
  
  n_treat   <- sum(data$treat == 1)
  n_control <- sum(data$treat == 0)
  n_outcome <- sum(data$treat == 0 & data$selection == 1)
  
  var_prop_treat   <- (prop_treat * (1 - prop_treat)) / max(n_treat, 1)
  var_prop_control <- (prop_control * (1 - prop_control)) / max(n_control, 1)
  var_prop_outcome <- (prop_outcome * (1 - prop_outcome)) / max(n_outcome, 1)
  
  hypothesis_1 <- 1 - prop_treat / prop_control
  hypothesis_2 <- prop_treat / prop_control
  
  g1_treat   <- -1 / prop_control
  g1_control <-  prop_treat / (prop_control^2)
  g1_outcome <- -1
  var_g1 <- g1_treat^2 * var_prop_treat + g1_control^2 * var_prop_control + g1_outcome^2 * var_prop_outcome
  se_g1  <- sqrt(var_g1)
  test_stat_h1 <- (prop_outcome - hypothesis_1) / se_g1
  
  g2_treat   <-  1 / prop_control
  g2_control <- -prop_treat / (prop_control^2)
  g2_outcome <- -1
  var_g2 <- g2_treat^2 * var_prop_treat + g2_control^2 * var_prop_control + g2_outcome^2 * var_prop_outcome
  se_g2  <- sqrt(var_g2)
  test_stat_h2 <- (prop_outcome - hypothesis_2) / se_g2
  
  p_value_h1 <- 2 * (1 - pnorm(abs(test_stat_h1)))
  p_value_h2 <- 2 * (1 - pnorm(abs(test_stat_h2)))
  c(p_value_h1 = p_value_h1, p_value_h2 = p_value_h2)
}

## -------------------------
## Cross-fitting (selection + outcome)
## -------------------------
crossfit_kfold <- function(df, K = 5, sel_form, y_form, weight_col = "weights") {
  set.seed(123)
  n <- nrow(df)
  fold_id <- sample(rep(1:K, length.out = n))
  s0 <- s1 <- yhat <- numeric(n)
  sel_vars <- intersect(colnames(df), all.vars(sel_form))
  
  for (k in 1:K) {
    train <- df[fold_id != k, , drop = FALSE]
    test  <- df[fold_id == k, , drop = FALSE]
    
    sel_fit <- glm(sel_form,
                   data   = train[, sel_vars, drop = FALSE],
                   family = binomial())
    
    t0 <- test; t0$treat <- 0
    t1 <- test; t1$treat <- 1
    s0[fold_id == k] <- predict(sel_fit, newdata = t0, type = "response")
    s1[fold_id == k] <- predict(sel_fit, newdata = t1, type = "response")
    
    ctrl_sel <- subset(train, treat == 0 & selection == 1)
    yfit <- .safe_binom_yfit(ctrl_sel, y_form, weight_col)
    
    if (yfit$kind == "glm") {
      yhat[fold_id == k] <- predict(yfit$fit, newdata = test, type = "response")
    } else {
      yhat[fold_id == k] <- yfit$p
    }
  }
  
  data.frame(s.0.hat = s0, s.1.hat = s1, y.hat = yhat)
}

## -------------------------
## Bootstrap for orthogonal estimator (re-estimates nuisances inside)
## -------------------------
main_bb_ortho <- function(sample_size, N_rep = 200, leedata, s.hat, y.hat,
                          flag_binary = TRUE, flag_helps = FALSE, flag_discrete = FALSE, ...) {
  ATE_bb <- matrix(NA_real_, N_rep, 2)
  for (b in 1:N_rep) {
    set.seed(b)
    inds <- sample(1:sample_size, sample_size, replace = TRUE)
    boot_leedata <- leedata[inds, ]
    
    sel_vars <- setdiff(colnames(boot_leedata), c("treat","selection","outcome","weights"))
    sel_form <- as.formula(paste0("selection ~ treat * (", paste0(sel_vars, collapse = "+"), ")"))
    y_form   <- as.formula(paste0("outcome ~ (", paste0(sel_vars, collapse = "+"), ")"))
    
    cf_preds <- crossfit_kfold(boot_leedata, K = 5, sel_form = sel_form, y_form = y_form, weight_col = "weights")
    s.hat.boot <- cf_preds[, c("s.0.hat","s.1.hat")]
    y.hat.boot <- cf_preds$y.hat
    
    resultb <- try(ortho_leebounds(
      leedata = boot_leedata, s.hat = s.hat.boot, y.hat = y.hat.boot,
      flag_binary = flag_binary, flag_discrete = flag_discrete,
      flag_helps = flag_helps, ortho_d = TRUE
    ), silent = TRUE)
    
    if (!inherits(resultb, "try-error")) {
      ATE_bb[b, ] <- c(resultb$lower_bound, resultb$upper_bound)
    }
  }
  ATE_bb
}
