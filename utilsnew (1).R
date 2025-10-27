##### ------------------------------------------------------------
##### Selection model: post-lasso (optional) + low-dim GLM (logit)
##### ------------------------------------------------------------
estimate_selection <- function(
    leedata,
    form = NULL,
    selection_function_name,
    variables_for_selection = NULL,
    names_to_include = c(),
    treat_name = "treat",         # fixed: no trailing '+'
    yname = "selection",
    myweights = NULL,
    ...
) {
  # read in data (not used directly below, but checked)
  d  <- leedata$treat
  s  <- leedata$selection
  sy <- leedata$outcome
  
  # covariates set
  if (is.null(variables_for_selection)) {
    variables_for_selection <- setdiff(
      colnames(leedata),
      c("treat","selection","outcome","X.Intercept.","(Intercept)","weights","prop1","prop0","myweights")
    )
  } else {
    variables_for_selection <- unique(setdiff(
      variables_for_selection,
      c("outcome","X.Intercept.","(Intercept)","weights","prop1","prop0","myweights")
    ))
  }
  
  # default formula if not provided
  if (is.null(form)) {
    if (length(variables_for_selection) > 0) {
      form <- as.formula(paste0("selection ~ ", treat_name, " * (",
                                paste0(variables_for_selection, collapse = "+"), ")"))
    } else {
      form <- selection ~ treat
    }
  }
  
  # weights
  if (is.null(myweights)) myweights <- rep(1, nrow(leedata))
  leedata$myweights <- myweights
  
  # optional: rlassologit selection to reduce X, then low-dim logit
  if (identical(selection_function_name, "rlassologit")) {
    glm.fit <- rlassologit(
      form,
      leedata[, c("treat","selection", variables_for_selection)],
      family = "binomial",
      ...
    )
    nz <- glm.fit$coefficients[glm.fit$coefficients != 0]
    selected_names <- setdiff(names(nz), c("(Intercept)", treat_name))
    selected_names <- unique(c(selected_names, names_to_include))
    
    # final low-dim formula
    if (length(selected_names) > 0) {
      form <- as.formula(paste0(yname, " ~ ", treat_name, " + ",
                                paste0(selected_names, collapse = "+")))
    } else {
      form <- as.formula(paste0(yname, " ~ ", treat_name))
    }
  }
  
  # final post-lasso low-dim GLM (logit) — now honoring weights
  glm.postlasso <- glm(
    formula = form,
    data    = leedata,        # <— DON'T subset here
    family  = binomial(),
    weights = leedata$myweights  # <— explicit vector that exists
  )
  
  glm.postlasso
}


##### --------------------------------------------
##### Predict s0(x)=P(S=1|D=0,X), s1(x)=P(S=1|D=1,X)
##### --------------------------------------------
predict_selection <- function(fit, leedata, ...) {
  leedata_0treat <- leedata; leedata_0treat$treat <- 0
  leedata_1treat <- leedata; leedata_1treat$treat <- 1
  s.0.hat <- predict(fit, leedata_0treat, type = "response")
  s.1.hat <- predict(fit, leedata_1treat, type = "response")
  list(s.0.hat = s.0.hat, s.1.hat = s.1.hat)
}


##### --------------------------------------------
##### Basic Lee bounds — Binary Y (plug-in + orthogonal)
##### --------------------------------------------
basic_lee_bound_binary <- function(
    leedata,                    # treat (0/1), selection (0/1), outcome (0/1), optional weights, prop1
    inds       = NULL,
    flag_helps = TRUE,          # TRUE=control-anchored, FALSE=treated-anchored, "auto"=choose by mean s0<=s1 (orthogonal only)
    ortho_d    = FALSE,
    s.hat      = NULL,          # required if ortho_d=TRUE (data.frame with s.0.hat,s.1.hat)
    y.hat      = NULL,          # required if ortho_d=TRUE: P(Y=1 | D=1,S=1,X)
    mu1        = NULL,          # optional P(D=1|X)
    eps        = 1e-6
){
  if (!is.null(inds)) leedata <- leedata[inds, , drop = FALSE]
  stopifnot(all(c("treat","selection","outcome") %in% names(leedata)))
  
  d <- as.numeric(leedata$treat)
  s <- as.numeric(leedata$selection)
  y <- leedata$outcome
  if (is.logical(y)) y <- as.numeric(y)
  if (is.factor(y))  y <- as.numeric(as.character(y))
  y <- pmin(pmax(y, 0), 1)
  w <- if ("weights" %in% names(leedata)) leedata$weights else rep(1, length(d))
  
  ok <- is.finite(d) & is.finite(s) & is.finite(y) & is.finite(w)
  d <- d[ok]; s <- s[ok]; y <- y[ok]; w <- w[ok]
  
  clip01 <- function(x) pmin(pmax(x, eps), 1 - eps)
  
  ## ---------- NON-ORTHOGONAL (plug-in) ----------
  if (!ortho_d) {
    s0 <- if (any(d == 0)) weighted.mean(s[d == 0] == 1, w[d == 0]) else 0
    s1 <- if (any(d == 1)) weighted.mean(s[d == 1] == 1, w[d == 1]) else 0
    
    s0c <- max(s0, eps); s1c <- max(s1, eps)
    p0  <- s0c / s1c
    
    if (p0 <= 1) {
      # Trim treated to control rate
      pY_t1 <- if (any(d==1 & s==1)) weighted.mean(y[d==1 & s==1] == 1, w[d==1 & s==1]) else 0
      EY0   <- if (any(d==0 & s==1)) weighted.mean(y[d==0 & s==1], w[d==0 & s==1]) else 0
      
      LB_treat <- max(1 - (1 - pY_t1) / p0, 0)
      UB_treat <- min(pY_t1 / p0, 1)
      
      lower <- LB_treat - EY0
      upper <- UB_treat - EY0
    } else {
      # Trim control to treated rate
      p0_inv <- 1 / p0
      pY_c1  <- if (any(d==0 & s==1)) weighted.mean(y[d==0 & s==1] == 1, w[d==0 & s==1]) else 0
      EY1    <- if (any(d==1 & s==1)) weighted.mean(y[d==1 & s==1], w[d==1 & s==1]) else 0
      
      LB_ctrl <- max(1 - (1 - pY_c1) / p0_inv, 0)
      UB_ctrl <- min(pY_c1 / p0_inv, 1)
      
      lower <- EY1 - UB_ctrl
      upper <- EY1 - LB_ctrl
    }
    
    if (lower > upper) { tmp <- lower; lower <- upper; upper <- tmp }
    return(list(lower_bound = lower, upper_bound = upper, mode = "non-orthogonal"))
  }
  
  ## ---------- ORTHOGONAL / DML ----------
  stopifnot(!is.null(s.hat), !is.null(y.hat))
  stopifnot(all(c("s.0.hat","s.1.hat") %in% names(s.hat)))
  s0hat <- clip01(as.numeric(s.hat$s.0.hat[ok]))
  s1hat <- clip01(as.numeric(s.hat$s.1.hat[ok]))
  yhat  <- clip01(as.numeric(y.hat[ok]))
  
  # auto-anchor if requested
  if (identical(flag_helps, "auto")) {
    flag_helps <- (weighted.mean(s0hat, w) <= weighted.mean(s1hat, w))
  } else {
    flag_helps <- isTRUE(flag_helps)
  }
  
  # propensities
  if (is.null(mu1)) {
    if ("prop1" %in% names(leedata)) mu1 <- as.numeric(leedata$prop1[ok]) else mu1 <- rep(mean(d), length(d))
  } else if (length(mu1) == 1L) {
    mu1 <- rep(mu1, length(d))
  } else {
    mu1 <- mu1[ok]
  }
  mu1  <- clip01(as.numeric(mu1))
  mu0  <- clip01(1 - mu1)
  mu11 <- clip01(mu1 * s1hat)
  
  if (flag_helps) {
    # control-anchored envelopes & corrections
    gL <- s0hat - s1hat * (1 - yhat)
    gU <- s1hat * yhat - s0hat
    I_L <- as.numeric(gL > 0)
    I_U <- as.numeric(gU < 0)
    
    base_nl <- pmax(gL, 0)
    base_nu <- pmin(gU, 0) + s0hat
    
    corr_s0_L <- ((1 - d) / mu0) * (s - s0hat) * I_L
    corr_s0_U <- ((1 - d) / mu0) * (s - s0hat) * (1 - I_U)
    
    corr_s1_L <- (d / mu1) * (s - s1hat) * (-(1 - yhat)) * I_L
    corr_s1_U <- (d / mu1) * (s - s1hat) * ( yhat ) * I_U
    
    corr_y_L  <- (d * s / mu11) * (y - yhat) * ( s1hat ) * I_L
    corr_y_U  <- (d * s / mu11) * (y - yhat) * (-s1hat) * I_U
    
    nl  <- weighted.mean(base_nl + corr_s0_L + corr_s1_L + corr_y_L, w)
    nu  <- weighted.mean(base_nu + corr_s0_U + corr_s1_U + corr_y_U, w)
    den <- weighted.mean(s0hat + ((1 - d) / mu0) * (s - s0hat), w)
  } else {
    # treated-anchored
    gL <- s1hat - s0hat * (1 - yhat)
    gU <- s0hat * yhat - s1hat
    I_L <- as.numeric(gL > 0)
    I_U <- as.numeric(gU < 0)
    
    base_nl <- pmax(gL, 0)
    base_nu <- pmin(gU, 0) + s1hat
    
    corr_s1_L <- (d / mu1) * (s - s1hat) * I_L
    corr_s1_U <- (d / mu1) * (s - s1hat) * (1 - I_U)
    
    corr_s0_L <- ((1 - d) / mu0) * (s - s0hat) * (-(1 - yhat)) * I_L
    corr_s0_U <- ((1 - d) / mu0) * (s - s0hat) * ( yhat ) * I_U
    
    corr_y_L  <- (d * s / mu11) * (y - yhat) * ( s0hat ) * I_L
    corr_y_U  <- (d * s / mu11) * (y - yhat) * (-s0hat) * I_U
    
    nl  <- weighted.mean(base_nl + corr_s1_L + corr_s0_L + corr_y_L, w)
    nu  <- weighted.mean(base_nu + corr_s1_U + corr_s0_U + corr_y_U, w)
    den <- weighted.mean(s1hat + (d / mu1) * (s - s1hat), w)
  }
  
  lower <- nl / max(den, eps)
  upper <- nu / max(den, eps)
  if (lower > upper) { tmp <- lower; lower <- upper; upper <- tmp }
  
  list(lower_bound = lower, upper_bound = upper, mode = "orthogonal",
       diagnostics = list(denominator = den, nl = nl, nu = nu))
}


##### --------------------------------------------
##### Denominator helper (plug-in & orthogonal)
##### --------------------------------------------
denominator <- function(
    s.hat, leedata,
    weights   = NULL,
    flag_helps = TRUE,   # TRUE control-anchored
    inds      = NULL,
    ortho_d   = FALSE,
    mu1       = NULL
) {
  if (!is.null(inds)) {
    s.hat   <- s.hat[inds, , drop = FALSE]
    leedata <- leedata[inds, , drop = FALSE]
  }
  
  d <- as.numeric(leedata$treat)
  s <- as.numeric(leedata$selection)
  
  if (is.null(weights)) {
    weights <- if ("weights" %in% names(leedata)) leedata$weights else rep(1, nrow(leedata))
  }
  
  clip01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)
  s0hat <- clip01(as.numeric(s.hat$s.0.hat))
  s1hat <- clip01(as.numeric(s.hat$s.1.hat))
  
  # Plug-in: use empirical E[S|D=a] (CONSISTENT across codebase)
  if (!ortho_d) {
    if (isTRUE(flag_helps)) {
      return(weighted.mean(s[d == 0] == 1, w = weights[d == 0]))
    } else {
      return(weighted.mean(s[d == 1] == 1, w = weights[d == 1]))
    }
  }
  
  # Orthogonal
  if (is.null(mu1)) {
    if ("prop1" %in% names(leedata)) mu1 <- as.numeric(leedata$prop1) else mu1 <- rep(mean(d), length(d))
  } else if (length(mu1) == 1L) {
    mu1 <- rep(mu1, length(d))
  }
  mu1 <- clip01(mu1); mu0 <- clip01(1 - mu1)
  
  if (isTRUE(flag_helps)) {
    base <- s0hat + ((1 - d) / mu0) * (s - s0hat)
  } else {
    base <- s1hat + (d / mu1) * (s - s1hat)
  }
  weighted.mean(base, w = weights)
}


##### --------------------------------------------
##### Numerator — Binary Y (plug-in & orthogonal)
##### --------------------------------------------
numerator_binary <- function(
    s.hat, y.hat, leedata,
    weights   = NULL,
    flag_helps = TRUE,   # TRUE control-anchored; "auto" allowed in orthogonal paths
    inds      = NULL,
    ortho_d   = FALSE,
    mu1       = NULL
) {
  if (!is.null(inds)) {
    s.hat   <- s.hat[inds, , drop = FALSE]
    y.hat   <- y.hat[inds]
    leedata <- leedata[inds, , drop = FALSE]
  }
  
  d  <- as.numeric(leedata$treat)
  s  <- as.numeric(leedata$selection)
  y  <- as.numeric(leedata$outcome)
  
  if (is.null(weights)) {
    weights <- if ("weights" %in% names(leedata)) leedata$weights else rep(1, nrow(leedata))
  }
  
  clip01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)
  s0hat <- clip01(as.numeric(s.hat$s.0.hat))
  s1hat <- clip01(as.numeric(s.hat$s.1.hat))
  yhat  <- clip01(as.numeric(y.hat))
  
  # Plug-in envelopes
  if (!ortho_d) {
    if (isTRUE(flag_helps)) {
      g1 <- s0hat - s1hat * (1 - yhat)        # lower
      g2 <- s1hat * yhat - s0hat              # upper
      nl <- weighted.mean(pmax(g1, 0), weights)
      nu <- weighted.mean(pmin(g2, 0) + s0hat, weights)
    } else {
      g1 <- s1hat - s0hat * (1 - yhat)
      g2 <- s0hat * yhat - s1hat
      nl <- weighted.mean(pmax(g1, 0), weights)
      nu <- weighted.mean(pmin(g2, 0) + s1hat, weights)
    }
    return(list(nl = nl, nu = nu))
  }
  
  # Orthogonal
  if (identical(flag_helps, "auto")) {
    flag_helps <- (weighted.mean(s0hat, weights) <= weighted.mean(s1hat, weights))
  } else {
    flag_helps <- isTRUE(flag_helps)
  }
  
  if (is.null(mu1)) {
    if ("prop1" %in% names(leedata)) {
      mu1 <- clip01(as.numeric(leedata$prop1))
    } else {
      mu1 <- rep(clip01(mean(d)), length(d))
    }
  } else if (length(mu1) == 1L) {
    mu1 <- rep(clip01(mu1), length(d))
  } else {
    mu1 <- clip01(as.numeric(mu1))
  }
  mu0  <- clip01(1 - mu1)
  mu11 <- clip01(mu1 * s1hat)
  
  if (flag_helps) {
    gL <- s0hat - s1hat * (1 - yhat)
    gU <- s1hat * yhat - s0hat
    I_L <- as.numeric(gL > 0)
    I_U <- as.numeric(gU < 0)
    
    base_nl <- pmax(gL, 0)
    base_nu <- pmin(gU, 0) + s0hat
    
    corr_s0_L <- ((1 - d) / mu0) * (s - s0hat) * I_L
    corr_s0_U <- ((1 - d) / mu0) * (s - s0hat) * (1 - I_U)
    
    corr_s1_L <- (d / mu1) * (s - s1hat) * (-(1 - yhat)) * I_L
    corr_s1_U <- (d / mu1) * (s - s1hat) * ( yhat ) * I_U
    
    corr_y_L  <- (d * s / mu11) * (y - yhat) * ( s1hat ) * I_L
    corr_y_U  <- (d * s / mu11) * (y - yhat) * (-s1hat) * I_U
    
    nl <- weighted.mean(base_nl + corr_s0_L + corr_s1_L + corr_y_L, weights)
    nu <- weighted.mean(base_nu + corr_s0_U + corr_s1_U + corr_y_U, weights)
  } else {
    gL <- s1hat - s0hat * (1 - yhat)
    gU <- s0hat * yhat - s1hat
    I_L <- as.numeric(gL > 0)
    I_U <- as.numeric(gU < 0)
    
    base_nl <- pmax(gL, 0)
    base_nu <- pmin(gU, 0) + s1hat
    
    corr_s1_L <- (d / mu1) * (s - s1hat) * I_L
    corr_s1_U <- (d / mu1) * (s - s1hat) * (1 - I_U)
    
    corr_s0_L <- ((1 - d) / mu0) * (s - s0hat) * (-(1 - yhat)) * I_L
    corr_s0_U <- ((1 - d) / mu0) * (s - s0hat) * ( yhat ) * I_U
    
    corr_y_L  <- (d * s / mu11) * (y - yhat) * ( s0hat ) * I_L
    corr_y_U  <- (d * s / mu11) * (y - yhat) * (-s0hat) * I_U
    
    nl <- weighted.mean(base_nl + corr_s1_L + corr_s0_L + corr_y_L, weights)
    nu <- weighted.mean(base_nu + corr_s1_U + corr_s0_U + corr_y_U, weights)
  }
  list(nl = nl, nu = nu)
}


##### --------------------------------------------
##### Numerator — Discrete Y (plug-in & orthogonal)
##### --------------------------------------------
numerator_discrete <- function(
    s.hat,                # data.frame s.0.hat, s.1.hat
    pi_hat,               # n x K matrix P(Y=beta|D=1,S=1,X)
    Y_levels,             # length K vector of support
    leedata,
    weights   = NULL,
    flag_helps = TRUE,    # TRUE control-anchored; "auto" allowed in orthogonal paths
    inds      = NULL,
    ortho_d   = FALSE,
    mu1       = NULL
){
  if (!is.null(inds)) {
    s.hat   <- s.hat[inds, , drop = FALSE]
    pi_hat  <- pi_hat[inds, , drop = FALSE]
    leedata <- leedata[inds, , drop = FALSE]
    if (!is.null(weights)) weights <- weights[inds]
  }
  stopifnot(all(c("treat","selection","outcome") %in% names(leedata)))
  d <- as.numeric(leedata$treat)
  s <- as.numeric(leedata$selection)
  y <- as.numeric(leedata$outcome)
  
  if (is.null(weights)) {
    weights <- if ("weights" %in% names(leedata)) leedata$weights else rep(1, nrow(leedata))
  }
  
  clip01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)
  s0hat <- clip01(as.numeric(s.hat$s.0.hat))
  s1hat <- clip01(as.numeric(s.hat$s.1.hat))
  
  # normalize pi_hat rows
  rs <- rowSums(pi_hat); rs[rs <= 0] <- 1
  pi_hat <- sweep(pi_hat, 1, rs, "/")
  
  n <- nrow(pi_hat); K <- length(Y_levels)
  stopifnot(nrow(pi_hat) == length(d), ncol(pi_hat) == K)
  
  # plug-in truncated expectations
  Ey_min <- sapply(seq_len(K), function(j) {
    wj <- pmin(Y_levels - Y_levels[j], 0)
    as.numeric(pi_hat %*% wj)
  })
  Ey_max <- sapply(seq_len(K), function(j) {
    wj <- pmax(Y_levels - Y_levels[j], 0)
    as.numeric(pi_hat %*% wj)
  })
  Ey_min <- matrix(Ey_min, nrow = n, ncol = K)
  Ey_max <- matrix(Ey_max, nrow = n, ncol = K)
  
  if (!ortho_d) {
    if (isTRUE(flag_helps)) {
      G_L <- sweep(matrix(rep(s0hat, each = K), nrow = n, byrow = FALSE), 2, Y_levels, `*`) + Ey_min * s1hat
      G_U <- sweep(matrix(rep(s0hat, each = K), nrow = n, byrow = FALSE), 2, Y_levels, `*`) + Ey_max * s1hat
    } else {
      G_L <- sweep(matrix(rep(s1hat, each = K), nrow = n, byrow = FALSE), 2, Y_levels, `*`) + Ey_min * s0hat
      G_U <- sweep(matrix(rep(s1hat, each = K), nrow = n, byrow = FALSE), 2, Y_levels, `*`) + Ey_max * s0hat
    }
    nl <- weighted.mean(apply(G_L, 1, max), w = weights)
    nu <- weighted.mean(apply(G_U, 1, min), w = weights)
    return(list(nl = nl, nu = nu))
  }
  
  # orthogonal: allow "auto" anchor
  if (identical(flag_helps, "auto")) {
    flag_helps <- (weighted.mean(s0hat, weights) <= weighted.mean(s1hat, weights))
  } else {
    flag_helps <- isTRUE(flag_helps)
  }
  
  if (is.null(mu1)) {
    if ("prop1" %in% names(leedata)) mu1 <- as.numeric(leedata$prop1) else mu1 <- rep(mean(d), length(d))
  } else if (length(mu1) == 1L) {
    mu1 <- rep(mu1, length(d))
  }
  mu1  <- clip01(mu1); mu0 <- clip01(1 - mu1)
  mu11 <- clip01(mu1 * s1hat)
  mu00 <- clip01(mu0 * s0hat)
  
  # build G_L, G_U for argmax/argmin beta selection
  if (flag_helps) {
    G_L <- sweep(matrix(rep(s0hat, each = K), nrow = n, byrow = FALSE), 2, Y_levels, `*`) + Ey_min * s1hat
    G_U <- sweep(matrix(rep(s0hat, each = K), nrow = n, byrow = FALSE), 2, Y_levels, `*`) + Ey_max * s1hat
  } else {
    G_L <- sweep(matrix(rep(s1hat, each = K), nrow = n, byrow = FALSE), 2, Y_levels, `*`) + Ey_min * s0hat
    G_U <- sweep(matrix(rep(s1hat, each = K), nrow = n, byrow = FALSE), 2, Y_levels, `*`) + Ey_max * s0hat
  }
  idx_L <- apply(G_L, 1, which.max)
  idx_U <- apply(G_U, 1, which.min)
  beta_star_L <- Y_levels[idx_L]
  beta_star_U <- Y_levels[idx_U]
  
  I_leq_L <- as.numeric(y <= beta_star_L)
  I_geq_U <- as.numeric(y >= beta_star_U)
  
  if (flag_helps) {
    rho_L <- (d * s / mu11) * (y - beta_star_L) * I_leq_L + beta_star_L * ((1 - d) * s / mu00)
    rho_U <- (d * s / mu11) * (y - beta_star_U) * I_geq_U + beta_star_U * ((1 - d) * s / mu00)
  } else {
    rho_L <- ((1 - d) * s / mu00) * (y - beta_star_L) * I_leq_L + beta_star_L * (d * s / mu11)
    rho_U <- ((1 - d) * s / mu00) * (y - beta_star_U) * I_geq_U + beta_star_U * (d * s / mu11)
  }
  
  nl <- weighted.mean(rho_L, w = weights)
  nu <- weighted.mean(rho_U, w = weights)
  list(nl = nl, nu = nu)
}


##### --------------------------------------------
##### Full discrete-Y Lee bounds wrapper
##### --------------------------------------------
lee_bounds_discrete <- function(
    leedata,
    s.hat,                     # data.frame with s.0.hat, s.1.hat
    pi_hat,                    # n x K matrix
    Y_levels,                  # support (ordered as columns of pi_hat)
    inds       = NULL,
    flag_helps = TRUE,         # TRUE / FALSE / "auto" (auto used only for orthogonal)
    ortho_d    = FALSE,
    mu1        = NULL,
    eps        = 1e-6
){
  if (!is.null(inds)) {
    leedata <- leedata[inds, , drop = FALSE]
    s.hat   <- s.hat[inds, , drop = FALSE]
    pi_hat  <- pi_hat[inds, , drop = FALSE]
  }
  stopifnot(all(c("treat","selection","outcome") %in% names(leedata)))
  d  <- as.numeric(leedata$treat)
  s  <- as.numeric(leedata$selection)
  y  <- as.numeric(leedata$outcome)
  w  <- if ("weights" %in% names(leedata)) leedata$weights else rep(1, length(d))
  
  ok <- is.finite(d) & is.finite(s) & is.finite(y) & is.finite(w)
  d <- d[ok]; s <- s[ok]; y <- y[ok]; w <- w[ok]
  s.hat  <- s.hat[ok, , drop = FALSE]
  pi_hat <- pi_hat[ok, , drop = FALSE]
  
  clip01 <- function(x) pmin(pmax(x, eps), 1 - eps)
  K <- length(Y_levels)
  stopifnot(nrow(pi_hat) == length(d), ncol(pi_hat) == K)
  stopifnot(all(c("s.0.hat","s.1.hat") %in% names(s.hat)))
  
  s0hat <- clip01(as.numeric(s.hat$s.0.hat))
  s1hat <- clip01(as.numeric(s.hat$s.1.hat))
  
  # row-normalize pi_hat
  rs <- rowSums(pi_hat); rs[rs <= 0] <- 1
  pi_hat <- sweep(pi_hat, 1, rs, "/")
  
  # build truncated expectations
  Ey_min <- sapply(seq_len(K), function(j) {
    wj <- pmin(Y_levels - Y_levels[j], 0)
    as.numeric(pi_hat %*% wj)
  })
  Ey_max <- sapply(seq_len(K), function(j) {
    wj <- pmax(Y_levels - Y_levels[j], 0)
    as.numeric(pi_hat %*% wj)
  })
  n <- length(d)
  Ey_min <- matrix(Ey_min, nrow = n, ncol = K)
  Ey_max <- matrix(Ey_max, nrow = n, ncol = K)
  
  BetaS0 <- outer(s0hat, Y_levels, `*`)
  BetaS1 <- outer(s1hat, Y_levels, `*`)
  
  if (!ortho_d) {
    if (isTRUE(flag_helps)) {
      G_L <- BetaS0 + (Ey_min * s1hat)
      G_U <- BetaS0 + (Ey_max * s1hat)
      den <- weighted.mean(s[d == 0] == 1, w[d == 0], na.rm = TRUE)  # plug-in: empirical E[S|D=0]
    } else {
      G_L <- BetaS1 + (Ey_min * s0hat)
      G_U <- BetaS1 + (Ey_max * s0hat)
      den <- weighted.mean(s[d == 1] == 1, w[d == 1], na.rm = TRUE)  # plug-in: empirical E[S|D=1]
    }
    nl <- weighted.mean(apply(G_L, 1, max), w, na.rm = TRUE)
    nu <- weighted.mean(apply(G_U, 1, min), w, na.rm = TRUE)
    lower <- nl / max(den, eps); upper <- nu / max(den, eps)
    if (lower > upper) { tmp <- lower; lower <- upper; upper <- tmp }
    return(list(lower_bound = lower, upper_bound = upper, mode = "non-orthogonal",
                diagnostics = list(nl = nl, nu = nu, denominator = den)))
  }
  
  # orthogonal: allow "auto" anchor
  if (identical(flag_helps, "auto")) {
    flag_helps <- (weighted.mean(s0hat, w) <= weighted.mean(s1hat, w))
  } else {
    flag_helps <- isTRUE(flag_helps)
  }
  
  if (flag_helps) {
    G_L <- BetaS0 + (Ey_min * s1hat)
    G_U <- BetaS0 + (Ey_max * s1hat)
  } else {
    G_L <- BetaS1 + (Ey_min * s0hat)
    G_U <- BetaS1 + (Ey_max * s0hat)
  }
  idx_L <- apply(G_L, 1, which.max)
  idx_U <- apply(G_U, 1, which.min)
  beta_star_L <- Y_levels[idx_L]
  beta_star_U <- Y_levels[idx_U]
  
  I_leq_L <- as.numeric(y <= beta_star_L)
  I_geq_U <- as.numeric(y >= beta_star_U)
  
  # propensities
  if (is.null(mu1)) {
    mu1 <- if ("prop1" %in% names(leedata)) as.numeric(leedata$prop1[ok]) else rep(mean(d), length(d))
  } else if (length(mu1) == 1L) {
    mu1 <- rep(mu1, length(d))
  } else {
    mu1 <- as.numeric(mu1)[ok]
  }
  mu1  <- clip01(mu1)
  mu0  <- clip01(1 - mu1)
  mu11 <- clip01(mu1 * s1hat)
  mu00 <- clip01(mu0 * s0hat)
  
  if (flag_helps) {
    rho_L <- (d * s / mu11) * (y - beta_star_L) * I_leq_L + beta_star_L * ((1 - d) * s / mu00)
    rho_U <- (d * s / mu11) * (y - beta_star_U) * I_geq_U + beta_star_U * ((1 - d) * s / mu00)
    den   <- weighted.mean(s0hat + ((1 - d) / mu0) * (s - s0hat), w, na.rm = TRUE)
  } else {
    rho_L <- ((1 - d) * s / mu00) * (y - beta_star_L) * I_leq_L + beta_star_L * (d * s / mu11)
    rho_U <- ((1 - d) * s / mu00) * (y - beta_star_U) * I_geq_U + beta_star_U * (d * s / mu11)
    den   <- weighted.mean(s1hat + (d / mu1) * (s - s1hat), w, na.rm = TRUE)
  }
  
  nl <- weighted.mean(rho_L, w, na.rm = TRUE)
  nu <- weighted.mean(rho_U, w, na.rm = TRUE)
  lower <- nl / max(den, eps); upper <- nu / max(den, eps)
  if (lower > upper) { tmp <- lower; lower <- upper; upper <- tmp }
  
  list(
    lower_bound = lower,
    upper_bound = upper,
    mode        = "orthogonal",
    diagnostics = list(nl = nl, nu = nu, denominator = den,
                       idx_beta_L = idx_L, idx_beta_U = idx_U)
  )
}


##### --------------------------------------------
##### Unified ortho/plug-in wrapper (binary or discrete)
##### --------------------------------------------
ortho_leebounds <- function(
    leedata,
    s.hat,
    # Binary:
    y.hat         = NULL,   # P(Y=1 | D=1,S=1,X)
    # Discrete:
    pi_hat        = NULL,   # n x K matrix
    Y_levels      = NULL,   # length K
    # Mode flags:
    flag_binary   = TRUE,
    flag_discrete = FALSE,
    flag_helps    = TRUE,   # TRUE / FALSE / "auto" (auto only used in orthogonal paths)
    ortho_d       = TRUE,
    # Other:
    weights       = NULL,
    mu1           = NULL,
    eps           = 1e-6,
    ...
) {
  d  <- as.numeric(leedata$treat)
  s  <- as.numeric(leedata$selection)
  sy <- as.numeric(leedata$outcome)
  
  if (is.null(weights)) {
    weights <- if ("weights" %in% names(leedata)) leedata$weights else rep(1, length(d))
  }
  
  if (!flag_binary && !flag_discrete)
    stop("Set either flag_binary=TRUE or flag_discrete=TRUE.")
  if (flag_binary && is.null(y.hat))
    stop("Binary mode needs y.hat = P(Y=1 | D=1,S=1,X).")
  if (flag_discrete && (is.null(pi_hat) || is.null(Y_levels)))
    stop("Discrete mode needs pi_hat (n x K) and Y_levels (length K).")
  
  # numerators
  if (flag_binary) {
    num <- numerator_binary(
      s.hat      = s.hat,
      y.hat      = y.hat,
      leedata    = leedata,
      weights    = weights,
      flag_helps = flag_helps,
      ortho_d    = ortho_d,
      mu1        = mu1,
      ...
    )
  } else {
    num <- numerator_discrete(
      s.hat      = s.hat,
      pi_hat     = pi_hat,
      Y_levels   = Y_levels,
      leedata    = leedata,
      weights    = weights,
      flag_helps = flag_helps,
      ortho_d    = ortho_d,
      mu1        = mu1,
      ...
    )
  }
  
  # denominator
  denom <- denominator(
    s.hat      = s.hat,
    leedata    = leedata,
    weights    = weights,
    flag_helps = flag_helps,
    ortho_d    = ortho_d,
    mu1        = mu1,
    ...
  )
  denom <- max(denom, eps)
  
  # convert to ATE bounds by subtracting/adding the anchor-arm observed mean on S=1
  if (isTRUE(flag_helps)) {
    EY0 <- if (any(d == 0 & s == 1)) weighted.mean(sy[d == 0 & s == 1], weights[d == 0 & s == 1], na.rm = TRUE) else 0
    lower_bound <- (num$nl / denom) - EY0
    upper_bound <- (num$nu / denom) - EY0
  } else if (identical(flag_helps, "auto")) {
    stop("Pass 'auto' only to lower-level calls (numerator/lee_bounds_discrete/basic_lee_bound_binary) or set TRUE/FALSE here explicitly.")
  } else {
    EY1 <- if (any(d == 1 & s == 1)) weighted.mean(sy[d == 1 & s == 1], weights[d == 1 & s == 1], na.rm = TRUE) else 0
    lower_bound <- EY1 - (num$nu / denom)
    upper_bound <- EY1 - (num$nl / denom)
  }
  
  if (lower_bound > upper_bound) { tmp <- lower_bound; lower_bound <- upper_bound; upper_bound <- tmp }
  list(lower_bound = lower_bound, upper_bound = upper_bound,
       diagnostics = list(denominator = denom, nl = num$nl, nu = num$nu))
}


##### --------------------------------------------
##### Simultaneous CI for (L,U) from bootstrap draws
##### --------------------------------------------
compute_confidence_region <- function(ATE_boot, ATE_est, ci_alpha = 0.05, tol = 1e-5) {
  if (sum(is.na(ATE_boot)) + sum(is.na(ATE_est)) > 0) {
    return(c(lower_bound = NA, upper_bound = NA))
  }
  ATE_boot_centered <- cbind(ATE_boot[,1] - ATE_est[1],
                             ATE_boot[,2] - ATE_est[2])
  
  Omega.hat <- matrix(0, 2, 2)
  Omega.hat[1,1] <- var(ATE_boot_centered[,1])
  Omega.hat[2,2] <- var(ATE_boot_centered[,2])
  Omega.hat[1,2] <- Omega.hat[2,1] <- cov(ATE_boot_centered[,1], ATE_boot_centered[,2])
  
  M <- expm::sqrtm(Omega.hat)
  if (max(abs(Im(M))) > tol) stop("Non-trivial imaginary part in sqrtm(Omega.hat)")
  M <- Re(M)
  
  z <- qnorm(sqrt(1 - ci_alpha))
  crit <- as.numeric(M %*% c(-z, +z))
  c(lower_bound = ATE_est[1] + crit[1],
    upper_bound = ATE_est[2] + crit[2])
}


##### --------------------------------------------
##### Bootstrap intersection test H0: theta0 ∈ [L,U]
##### --------------------------------------------
intersection_test_boot <- function(
    leedata,
    bound_fun,           # function(df, ...) -> list(lower_bound, upper_bound)
    ...,                 # forwarded to bound_fun (e.g., s.hat, y.hat, flags)
    theta0    = 0,
    B         = 500,
    seed      = 12345,
    ci_alpha  = 0.05
) {
  stopifnot(is.data.frame(leedata))
  
  # point estimate
  est  <- bound_fun(leedata, ...)
  Lhat <- est$lower_bound
  Uhat <- est$upper_bound
  inside <- (Lhat <= theta0 && theta0 <= Uhat)
  
  # bootstrap
  set.seed(seed)
  n <- nrow(leedata)
  boot_mat <- matrix(NA_real_, B, 2)
  for (b in 1:B) {
    inds <- sample.int(n, n, replace = TRUE)
    resb <- try(bound_fun(leedata[inds, , drop = FALSE], ...), silent = TRUE)
    if (!inherits(resb, "try-error") && all(c("lower_bound","upper_bound") %in% names(resb))) {
      boot_mat[b, ] <- c(resb$lower_bound, resb$upper_bound)
    }
  }
  ok <- stats::complete.cases(boot_mat)
  boot_mat <- boot_mat[ok, , drop = FALSE]
  if (nrow(boot_mat) < max(50, 0.5*B)) warning("Many bootstrap failures; results may be unstable.")
  
  # studentized max one-sided stat
  sd_L <- stats::sd(boot_mat[,1] - Lhat); if (!is.finite(sd_L) || sd_L <= 0) sd_L <- 1e-8
  sd_U <- stats::sd(boot_mat[,2] - Uhat); if (!is.finite(sd_U) || sd_U <= 0) sd_U <- 1e-8
  
  T_obs <- max((Lhat - theta0)/sd_L, (theta0 - Uhat)/sd_U, 0)
  
  L_dev <- (boot_mat[,1] - Lhat)/sd_L
  U_dev <- (Uhat - boot_mat[,2])/sd_U
  T_b   <- pmax(L_dev, U_dev, 0)
  
  pval <- mean(T_b >= T_obs, na.rm = TRUE)
  ci <- compute_confidence_region(ATE_boot = boot_mat, ATE_est = c(Lhat, Uhat), ci_alpha = ci_alpha)
  
  list(
    theta0         = theta0,
    Lhat           = Lhat,
    Uhat           = Uhat,
    inside_bounds  = inside,
    T_obs          = T_obs,
    p_value        = pval,
    sim_CI         = ci,
    boot_draws     = boot_mat
  )
}

GetBounds <- function(x) c(x$lower_bound, x$upper_bound)


##### --------------------------------------------
##### Unified bootstrap across (binary/discrete) × (plug-in/orthogonal)
##### --------------------------------------------
main_bb_unified <- function(
    leedata,
    N_rep       = 500,
    mode        = c("binary","discrete"),
    ortho       = TRUE,
    flag_helps  = TRUE,     # TRUE / FALSE / "auto" (auto honored in orthogonal paths)
    mu1         = NULL,
    Y_levels    = NULL,
    seed        = 12345,
    sel_vars_exclude = c("treat","selection","outcome","weights","prop1"),
    ...
) {
  mode <- match.arg(mode)
  n <- nrow(leedata)
  ATE_bb <- matrix(NA_real_, N_rep, 2, dimnames = list(NULL, c("lower","upper")))
  set.seed(seed)
  
  has_w <- "weights" %in% names(leedata)
  
  for (b in 1:N_rep) {
    inds <- sample.int(n, n, replace = TRUE)
    dfb  <- leedata[inds, , drop = FALSE]
    wb   <- if (has_w) dfb$weights else rep(1, nrow(dfb))
    
    # X's for models
    xvars <- setdiff(colnames(dfb), sel_vars_exclude)
    has_x <- length(xvars) > 0
    
    form_selection <- if (has_x)
      as.formula(paste0("selection ~ treat * (", paste0(xvars, collapse = "+"), ")"))
    else selection ~ treat
    
    form_outcome_bin <- if (has_x)
      as.formula(paste0("outcome ~ ", paste0(xvars, collapse = "+")))
    else outcome ~ 1
    
    # -------- Plug-in --------
    if (!ortho) {
      if (mode == "binary") {
        res <- try(basic_lee_bound_binary(dfb, ortho_d = FALSE), silent = TRUE)
        
      } else {
        if (is.null(Y_levels)) stop("Provide Y_levels for discrete mode (plug-in).")
        
        # selection model
        sel_fit <- try(glm(form_selection, data = dfb, family = binomial(),
                           weights = wb, na.action = na.exclude), silent = TRUE)
        if (inherits(sel_fit, "try-error")) { ATE_bb[b,] <- c(NA,NA); next }
        
        dfb0 <- dfb; dfb0$treat <- 0
        dfb1 <- dfb; dfb1$treat <- 1
        s0hat <- as.numeric(predict(sel_fit, newdata = dfb0, type = "response"))
        s1hat <- as.numeric(predict(sel_fit, newdata = dfb1, type = "response"))
        s.hat.boot <- data.frame(`s.0.hat` = s0hat, `s.1.hat` = s1hat)
        
        # multinomial outcome on D=1,S=1
        dsub <- dfb$treat == 1 & dfb$selection == 1
        if (sum(dsub) < max(5, length(unique(dfb$outcome[dsub])))) { ATE_bb[b,] <- c(NA,NA); next }
        
        sub_df <- dfb[dsub, , drop = FALSE]
        sub_df$y_fac_sub <- factor(sub_df$outcome, levels = Y_levels)
        form_outcome_mult <- if (has_x)
          as.formula(paste0("y_fac_sub ~ ", paste0(xvars, collapse = "+"))) else y_fac_sub ~ 1
        
        mfit <- try(nnet::multinom(form_outcome_mult, data = sub_df,
                                   weights = wb[dsub], trace = FALSE), silent = TRUE)
        if (inherits(mfit, "try-error")) { ATE_bb[b,] <- c(NA,NA); next }
        
        pi_hat.boot <- predict(mfit, newdata = dfb, type = "probs")
        if (is.vector(pi_hat.boot)) pi_hat.boot <- cbind(pi_hat.boot, 1 - pi_hat.boot)
        pi_hat.boot <- as.matrix(pi_hat.boot)
        
        res <- try(
          lee_bounds_discrete(
            leedata    = dfb,
            s.hat      = s.hat.boot,
            pi_hat     = pi_hat.boot,
            Y_levels   = Y_levels,
            flag_helps = flag_helps,
            ortho_d    = FALSE,
            ...
          ),
          silent = TRUE
        )
      }
      
      if (!inherits(res, "try-error") && all(c("lower_bound","upper_bound") %in% names(res))) {
        ATE_bb[b,] <- c(res$lower_bound, res$upper_bound)
      }
      next
    }
    
    # -------- Orthogonal (DML) --------
    sel_fit <- try(glm(form_selection, data = dfb, family = binomial(),
                       weights = wb, na.action = na.exclude), silent = TRUE)
    if (inherits(sel_fit, "try-error")) { ATE_bb[b,] <- c(NA,NA); next }
    
    dfb0 <- dfb; dfb0$treat <- 0
    dfb1 <- dfb; dfb1$treat <- 1
    s0hat <- as.numeric(predict(sel_fit, newdata = dfb0, type = "response"))
    s1hat <- as.numeric(predict(sel_fit, newdata = dfb1, type = "response"))
    s.hat.boot <- data.frame(`s.0.hat` = s0hat, `s.1.hat` = s1hat)
    
    if (mode == "binary") {
      dsub <- dfb$treat == 1 & dfb$selection == 1
      if (sum(dsub) < 5) { ATE_bb[b,] <- c(NA,NA); next }
      
      y_fit <- try(glm(form_outcome_bin, data = dfb[dsub, , drop = FALSE],
                       family = binomial(), weights = wb[dsub], na.action = na.exclude),
                   silent = TRUE)
      if (inherits(y_fit, "try-error")) { ATE_bb[b,] <- c(NA,NA); next }
      y.hat.boot <- as.numeric(predict(y_fit, newdata = dfb, type = "response"))
      
      res <- try(
        ortho_leebounds(
          leedata     = dfb,
          s.hat       = s.hat.boot,
          y.hat       = y.hat.boot,
          flag_binary = TRUE,
          flag_helps  = flag_helps,   # TRUE/FALSE/"auto"
          ortho_d     = TRUE,
          weights     = wb,
          mu1         = mu1,
          ...
        ),
        silent = TRUE
      )
      
    } else {  # discrete
      if (is.null(Y_levels)) stop("Provide Y_levels for discrete mode.")
      dsub <- dfb$treat == 1 & dfb$selection == 1
      if (sum(dsub) < max(5, length(unique(dfb$outcome[dsub])))) { ATE_bb[b,] <- c(NA,NA); next }
      
      sub_df <- dfb[dsub, , drop = FALSE]
      sub_df$y_fac_sub <- factor(sub_df$outcome, levels = Y_levels)
      form_outcome_mult <- if (has_x)
        as.formula(paste0("y_fac_sub ~ ", paste0(xvars, collapse = "+"))) else y_fac_sub ~ 1
      
      mfit <- try(nnet::multinom(form_outcome_mult, data = sub_df,
                                 weights = wb[dsub], trace = FALSE), silent = TRUE)
      if (inherits(mfit, "try-error")) { ATE_bb[b,] <- c(NA,NA); next }
      
      pi_hat.boot <- predict(mfit, newdata = dfb, type = "probs")
      if (is.vector(pi_hat.boot)) pi_hat.boot <- cbind(pi_hat.boot, 1 - pi_hat.boot)
      pi_hat.boot <- as.matrix(pi_hat.boot)
      
      res <- try(
        lee_bounds_discrete(
          leedata    = dfb,
          s.hat      = s.hat.boot,
          pi_hat     = pi_hat.boot,
          Y_levels   = Y_levels,
          flag_helps = flag_helps,  # TRUE/FALSE/"auto"
          ortho_d    = TRUE,
          mu1        = mu1,
          ...
        ),
        silent = TRUE
      )
    }
    
    if (!inherits(res, "try-error") && all(c("lower_bound","upper_bound") %in% names(res))) {
      ATE_bb[b,] <- c(res$lower_bound, res$upper_bound)
    }
  }
  
  ATE_bb
}
