### please explain what this function does
estimate_selection<-function(leedata,form=NULL,selection_function_name,variables_for_selection=NULL,names_to_include=c(),
                             treat_name="treat+",yname="selection",myweights=NULL,...) {
  ### read in data
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  ### set variables for regression
  if (is.null(variables_for_selection)) {
    variables_for_selection<-setdiff(colnames(leedata),c("treat","selection","outcome","X.Intercept.","(Intercept)","weights","prop1","prop0"))
  } else {
    variables_for_selection<-unique(setdiff(variables_for_selection,c("outcome","X.Intercept.","(Intercept)","weights","prop1","prop0")))
  }
  if (is.null(form)) {
    
    form<-as.formula(paste0("selection~(treat)*(", paste0(variables_for_selection,collapse="+"),")"))
  }
  
  #print(form)
  
  if (is.null(myweights)) {
    myweights<-rep(1,dim(leedata)[1])
  }
  
  ## if rlassologit, select covariates
  if (selection_function_name=="rlassologit") {
    glm.fit<-rlassologit(form, leedata[,c("treat", "selection", variables_for_selection)],family="binomial",...)
    # select non-zero coefficients
    non_zero_coefs<-glm.fit$coefficients[glm.fit$coefficients!=0]
    # names whose coefs are non-zero
    selected_names<-setdiff(names(non_zero_coefs),c("(Intercept)","treat"))
    # add manually selected features
    selected_names<-unique(c(selected_names,names_to_include))
    ## (optional): add raw variables after interactions
    #grep("treat:",selected_names,value=TRUE)
    
    # if treatment was dropped, make sure to re-run low-dim analysis with treatment
    ## this corresponds to assumption that Pr (S(1)=S(0)=1)<1 (not everyone is an always-taker)
    if (length(selected_names)>0) {
      form<-as.formula(paste0(yname,"~",treat_name,paste0(selected_names,collapse="+")))
    } else {
      form<-as.formula(paste0(yname,"~",treat_name))
    }
  }
  
  ### final stage is always logistic with low-dim covariates
  leedata$myweights<-myweights
  glm.postlasso<-glm( form,data=leedata[,c("treat", "selection", variables_for_selection)],family="binomial")
  
  return(glm.postlasso)
}

### please explain what this function does
predict_selection<-function(fit,leedata,...) {
  leedata_0treat<-leedata
  leedata_0treat$treat<-0
  
  leedata_1treat<-leedata
  leedata_1treat$treat<-1
  
  s.0.hat<-predict( fit,leedata_0treat,type="response")
  s.1.hat<-predict( fit,leedata_1treat,type="response")
  
  return(list(s.0.hat=s.0.hat,s.1.hat=s.1.hat))
}


basic_lee_bound_binary<-function(leedata,inds=NULL,...) {
  
  if (!is.null(inds)) {
    leedata<-leedata[inds,]
  }
  
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  if ("weights" %in% colnames(leedata)) {
    weights<-leedata$weights
  } else {
    weights<-rep(1,length(d))
  }
  
  treat_size<-sum(d==1)
  control_size<-sum(d==0)
  prop_control_nonmissing<-weighted.mean(s[d==0]==1,weights[d==0])
  prop_treat_nonmissing<-weighted.mean(s[d==1]==1,weights[d==1])
  
  p0<-prop_control_nonmissing/prop_treat_nonmissing
  #print(p0)

 
  if (p0<=1) {
    p0_star = p0
    pY<-weighted.mean(sy[d==1 & s==1]==1,weights[d==1 & s==1])
    
    lower_bound=max(1- (1-pY)/p0,0) - weighted.mean(sy[d==0 & s==1],weights[d==0 & s==1])
    
    upper_bound=min(pY/p0-1,0)+1 - weighted.mean(sy[d==0 & s==1],weights[d==0 & s==1])
    
    
  } else {
    p0_star = 1/p0
   
    pY<-weighted.mean(sy[d==0 & s==1]==1,weights[d==0 & s==1])
    
    lower_bound= weighted.mean(sy[d==1 & s==1],weights[d==1 & s==1]) - (min(pY/p0_star -1,0)+1)
    
    upper_bound=weighted.mean(sy[d==1 & s==1],weights[d==1 & s==1])  - (max(1- (1-pY)/p0_star ,0))
    
  }
  

  return(list( lower_bound= lower_bound,
               upper_bound=upper_bound,
               s0=prop_control_nonmissing,
               s1=prop_treat_nonmissing,
               p0_star= p0_star,
               pY=pY))
}


denominator<-function(s.hat,leedata,weights=NULL,flag_helps=FALSE,inds=NULL,ortho_d=FALSE, ...) {
  
  if (is.null(weights)) {
    weights<-rep(1,dim(s.hat)[1])
  }
  
  
  if (!is.null(inds)) {
    s.hat<-s.hat[inds,]
    leedata<-leedata[inds,]
  }

  
  if (!ortho_d) {
    s<-leedata$selection
    d<-leedata$treat
    weights<-leedata$weights
    
  if (flag_helps) {
    denom<-weighted.mean(s[d==0]==1,weights[d==0])
  } else {
    denom<-weighted.mean(s[d==1]==1,weights[d==1])
  }
  }
  if (ortho_d) {
    
    d<-leedata$treat
    s<-leedata$selection
    if ("weights" %in% colnames(leedata)) {
      weights<-leedata$weights
    } else {
      weights<-rep(1,length(d))
    }
    
    if (flag_helps) {
      denom<-weighted.mean(s.hat$s.0.hat+weights*(1-d)*(s-s.hat$s.0.hat),weights)
    } else {
      denom<-weighted.mean(s.hat$s.1.hat+weights*d*(s-s.hat$s.1.hat),weights)
    }
  }
  return(denom)
}


###  please add orthogonal correction this function.
### 
###  numerator_binary function with orthogonal correction based on envelope score estimator

numerator_binary <- function(s.hat, y.hat, weights = NULL, leedata = NULL, flag_helps = FALSE, 
                             inds = NULL, ortho_d = FALSE, ...) {
  
  if (!is.null(inds)) {
    s.hat <- s.hat[inds, ]
    y.hat <- y.hat[inds]
    leedata <- leedata[inds, ]
  }
  
  if (is.null(weights)) {
    weights <- rep(1, dim(s.hat)[1])
  }
  
  # Extract data components
  d <- leedata$treat
  s <- leedata$selection
  sy <- leedata$outcome
  
  if ("weights" %in% colnames(leedata)) {
    weights <- leedata$weights
  } else {
    weights <- rep(1, length(d))
  }
  
  if (!ortho_d) {
    # Original non-orthogonal version
    if (flag_helps) {
      nl <- weighted.mean(sapply(s.hat$s.0.hat - s.hat$s.1.hat * (1 - y.hat), max, 0), weights)
      nu <- weighted.mean(sapply(s.hat$s.1.hat * y.hat - s.hat$s.0.hat, min, 0) + s.hat$s.0.hat, weights)
    } else {
      nl <- weighted.mean(sapply(s.hat$s.1.hat - s.hat$s.0.hat * (1 - y.hat), max, 0), weights)
      nu <- weighted.mean(sapply(s.hat$s.0.hat * y.hat - s.hat$s.1.hat, min, 0) + s.hat$s.1.hat, weights)
    }
  } else {
    # Orthogonal correction version using envelope score approach
    # This follows the debiased ML framework for aggregated intersection bounds
    
    if (flag_helps) {
      # For the "helps" case (control group bounds)
      
      # Define the envelope functions
      g1 <- s.hat$s.0.hat - s.hat$s.1.hat * (1 - y.hat)  # for lower bound
      g2 <- s.hat$s.1.hat * y.hat - s.hat$s.0.hat        # for upper bound
      
      # Indicators for which function achieves the min/max
      I1_l <- as.numeric(g1 > 0)  # Indicator for lower bound constraint being active
      I2_u <- as.numeric(g2 < 0)  # Indicator for upper bound constraint being active
      
      # Base envelope terms
      base_nl <- pmax(g1, 0)
      base_nu <- pmin(g2, 0) + s.hat$s.0.hat
      
      # Orthogonal corrections - these are the key innovation
      # For selection probabilities s.0.hat
      ortho_s0_nl <- (1 - d) * (s - s.hat$s.0.hat) * I1_l
      ortho_s0_nu <- (1 - d) * (s - s.hat$s.0.hat) * (1 - I2_u)
      
      # For selection probabilities s.1.hat  
      ortho_s1_nl <- d * (s - s.hat$s.1.hat) * I1_l * (-(1 - y.hat))
      ortho_s1_nu <- d * (s - s.hat$s.1.hat) * I2_u * y.hat
      
      # For outcome predictions y.hat - this is the envelope-specific correction
      ortho_y_nl <- (1 - d) * s * (sy - y.hat) * I1_l * s.hat$s.1.hat
      ortho_y_nu <- (1 - d) * s * (sy - y.hat) * I2_u * (-s.hat$s.1.hat)
      
      # Combine base terms with orthogonal corrections
      nl <- weighted.mean(base_nl + ortho_s0_nl + ortho_s1_nl + ortho_y_nl, weights)
      nu <- weighted.mean(base_nu + ortho_s0_nu + ortho_s1_nu + ortho_y_nu, weights)
      
    } else {
      # For the standard case (treatment group bounds)
      
      # Define the envelope functions
      g1 <- s.hat$s.1.hat - s.hat$s.0.hat * (1 - y.hat)  # for lower bound
      g2 <- s.hat$s.0.hat * y.hat - s.hat$s.1.hat        # for upper bound
      
      # Indicators for which function achieves the min/max
      I1_l <- as.numeric(g1 > 0)  # Indicator for lower bound constraint being active
      I2_u <- as.numeric(g2 < 0)  # Indicator for upper bound constraint being active
      
      # Base envelope terms
      base_nl <- pmax(g1, 0)
      base_nu <- pmin(g2, 0) + s.hat$s.1.hat
      
      # Orthogonal corrections following the envelope score approach
      # For selection probabilities s.1.hat
      ortho_s1_nl <- d * (s - s.hat$s.1.hat) * I1_l
      ortho_s1_nu <- d * (s - s.hat$s.1.hat) * (1 - I2_u)
      
      # For selection probabilities s.0.hat
      ortho_s0_nl <- (1 - d) * (s - s.hat$s.0.hat) * I1_l * (-(1 - y.hat))
      ortho_s0_nu <- (1 - d) * (s - s.hat$s.0.hat) * I2_u * y.hat
      
      # For outcome predictions y.hat - envelope-specific correction
      ortho_y_nl <- (1 - d) * s * (sy - y.hat) * I1_l * s.hat$s.0.hat
      ortho_y_nu <- (1 - d) * s * (sy - y.hat) * I2_u * (-s.hat$s.0.hat)
      
      # Combine base terms with orthogonal corrections
      nl <- weighted.mean(base_nl + ortho_s1_nl + ortho_s0_nl + ortho_y_nl, weights)
      nu <- weighted.mean(base_nu + ortho_s1_nu + ortho_s0_nu + ortho_y_nu, weights)
    }
  }
  
  return(list(nl = nl, nu = nu))
}


numerator_discrete<-function(s.hat,y.hat,leedata=NULL,flag_helps=FALSE,inds=NULL,...) {
  ###  please complete this function
  
}



basic_lee_bound_discrete<-function(leedata,inds=NULL,...) {
  ##  please complete this function. when the outcome is binary, the function
  ## should return the same output as basic_lee_bound_binary
  ## it should work for any finitely supported outcome
}

ortho_leebounds <- function(leedata, s.hat, y.hat, flag_binary = TRUE, flag_helps = FALSE, 
                            flag_discrete = FALSE, ortho_d = TRUE, weights = NULL, ...) {
  
  d <- leedata$treat
  s <- leedata$selection
  sy <- leedata$outcome
  
  if (is.null(weights)) {
    if ("weights" %in% colnames(leedata)) {
      weights <- leedata$weights
    } else {
      weights <- rep(1, length(d))
    }
  }
  
  if (flag_binary) {
    num <- numerator_binary(
      s.hat = s.hat,
      y.hat = y.hat,
      leedata = leedata,
      flag_helps = flag_helps,
      weights = weights,
      ortho_d = ortho_d,
      ...
    )
  } else {
    if (flag_discrete) {
      num <- numerator_discrete(
        s.hat = s.hat,
        y.hat = y.hat,
        leedata = leedata,
        flag_helps = flag_helps,
        weights = weights,
        ortho_d = ortho_d,
        ...
      )
    }
  }
  
  denom <- denominator(
    s.hat = s.hat,
    leedata = leedata,
    flag_helps = flag_helps,
    weights = weights,
    ortho_d = ortho_d,
    ...
  )
  
  lower_bound <- num$nl / denom
  upper_bound <- num$nu / denom
  
  if (flag_helps) {
    lower_bound <- num$nl / denom - weighted.mean(sy[d == 0 & s == 1], weights[d == 0 & s == 1])
    upper_bound <- num$nu / denom - weighted.mean(sy[d == 0 & s == 1], weights[d == 0 & s == 1])
  } else {
    lower_bound <- weighted.mean(sy[d == 1 & s == 1], weights[d == 1 & s == 1]) - num$nu / denom
    upper_bound <- weighted.mean(sy[d == 1 & s == 1], weights[d == 1 & s == 1]) - num$nl / denom
  }
  
  return(list(lower_bound = lower_bound, upper_bound = upper_bound))
}


main_bb<-function(function_name,sample_size,N_rep=10,...) {
  
  ATE_bb<-matrix(0,N_rep,2)
  for (b in 1:N_rep) {
    set.seed(b)
    #print(b)
    
    inds<-sample(1:sample_size,sample_size,replace=TRUE)
     resultb = try(function_name(inds,...))
    ATE_bb[b,]<-c(resultb$lower_bound, resultb$upper_bound)
  }
  
  
  return(ATE_bb)
}


compute_confidence_region<-function(ATE_boot, ATE_est, ci_alpha=0.05,tol=1e-5) {
  Omega.hat<-matrix(0,2,2)
  if (sum(is.na(ATE_boot))+sum(is.na(ATE_est))>0) {
    return(c(lower_bound = NA, upper_bound=NA))
  }
  ATE_boot_centered<-matrix(0,dim(ATE_boot)[1],2)
  ## Centered draws of lower bound
  ATE_boot_centered[,1]<-ATE_boot[,1]-ATE_est[1]
  ## Centered draws of upper bound
  ATE_boot_centered[,2]<-ATE_boot[,2]-ATE_est[2]
  
  Omega.hat[1,1]<-var( ATE_boot_centered[,1])
  Omega.hat[2,2]<-var(ATE_boot_centered[,2])
  Omega.hat[1,2]<-cov(ATE_boot_centered[,1],ATE_boot_centered[,2])
  Omega.hat[2,1]<-Omega.hat[1,2]
  
  crit.val<-sqrtm(Omega.hat)%*% c(-qnorm(sqrt(1-ci_alpha)), qnorm(sqrt(1-ci_alpha)) ) 
  if (max(abs(Im(sqrtm(Omega.hat))))>tol) {
    stop ("Non-trivial imaginary part!")
  } else {
    crit.val<-sapply( crit.val,Re)
    
  }
  lower_bound<-ATE_est[1]+ crit.val[1]
  upper_bound<-ATE_est[2] +crit.val[2]
  return(c(lower_bound = lower_bound, upper_bound=upper_bound))
}

GetBounds<-function(x) {
  return(c(x$lower_bound,x$upper_bound))
}


intersection_test <- function(data, weights) {
  # Ensure weights length matches the number of rows in the data
  weights <- weights / sum(weights)
  
  # Calculate weighted proportions
  prop_treat <- sum(weights[data$treat == 1] * data$selection[data$treat == 1]) / sum(weights[data$treat == 1])
  prop_control <- sum(weights[data$treat == 0] * data$selection[data$treat == 0]) / sum(weights[data$treat == 0])
  prop_outcome <- sum(weights[data$treat == 0 & data$selection == 1] * data$outcome[data$treat == 0 & data$selection == 1]) / sum(weights[data$treat == 0 & data$selection == 1])
  
  # Compute sample sizes for each group
  n_treat <- sum(data$treat == 1)
  n_control <- sum(data$treat == 0)
  n_outcome <- sum(data$treat == 0 & data$selection == 1)
  
  # Adjust variances for proportions (Bernoulli distribution)
  var_prop_treat <- (prop_treat * (1 - prop_treat)) / n_treat
  var_prop_control <- (prop_control * (1 - prop_control)) / n_control
  var_prop_outcome <- (prop_outcome * (1 - prop_outcome)) / n_outcome
  
  # Define the hypotheses
  hypothesis_1 <- 1 - prop_treat / prop_control
  hypothesis_2 <- prop_treat / prop_control
  
  # Delta method partial derivatives for hypothesis 1
  g1_treat <- -1 / prop_control
  g1_control <- prop_treat / (prop_control^2)
  g1_outcome <- -1
  
  # Variance of g1(prop_treat, prop_control, prop_outcome)
  var_g1 <- g1_treat^2 * var_prop_treat + g1_control^2 * var_prop_control + g1_outcome^2 * var_prop_outcome
  se_g1 <- sqrt(var_g1)
  
  # Test statistic for hypothesis 1
  test_stat_h1 <- (prop_outcome - hypothesis_1) / se_g1
  
  # Delta method partial derivatives for hypothesis 2
  g2_treat <- 1 / prop_control
  g2_control <- -prop_treat / (prop_control^2)
  g2_outcome <- -1
  
  # Variance of g2(prop_treat, prop_control, prop_outcome)
  var_g2 <- g2_treat^2 * var_prop_treat + g2_control^2 * var_prop_control + g2_outcome^2 * var_prop_outcome
  se_g2 <- sqrt(var_g2)
  
  # Test statistic for hypothesis 2
  test_stat_h2 <- (prop_outcome - hypothesis_2) / se_g2
  
  # Compute p-values for both hypotheses
  p_value_h1 <- 2 * (1 - pnorm(abs(test_stat_h1)))
  p_value_h2 <- 2 * (1 - pnorm(abs(test_stat_h2)))
  
  # Return results
  return(c(p_value_h1 = p_value_h1, p_value_h2 = p_value_h2))
}


main_bb_ortho <- function(sample_size, N_rep = 10, leedata, s.hat, y.hat, flag_binary = TRUE, flag_helps = FALSE, ...) {
  
  ATE_bb <- matrix(0, N_rep, 2)
  
  for (b in 1:N_rep) {
    set.seed(b)
    
    # Bootstrap sample indices
    inds <- sample(1:sample_size, sample_size, replace = TRUE)
    
    # Bootstrap the data
    boot_leedata <- leedata[inds, ]
    boot_s.hat <- s.hat[inds, ]
    boot_y.hat <- y.hat[inds]
    
    # Re-estimate nuisance parameters on bootstrap sample
    # Selection model
    tryCatch({
      form_selection <- as.formula(paste0(
        "selection ~ treat * (", paste0(setdiff(colnames(boot_leedata), c("treat", "selection", "outcome", "weights")), collapse = "+"), ")"
      ))
      
      glm_fit_boot <- estimate_selection(leedata = boot_leedata, form = form_selection, selection_function_name = "glm")
      s.hat.boot <- as.data.frame(predict_selection(glm_fit_boot, leedata = boot_leedata))
      
      # Outcome model
      form_outcome <- as.formula(paste0(
        "outcome ~ (", paste0(setdiff(colnames(boot_leedata), c("treat", "selection", "outcome", "weights")), collapse = "+"), ")"
      ))
      
      glm_y_fit_boot <- glm(
        form_outcome,
        data = boot_leedata[boot_leedata$treat == 0 & boot_leedata$selection == 1, ]
      )
      y.hat.boot <- predict(glm_y_fit_boot, boot_leedata, type = "response")
      
      # Compute orthogonal bounds
      resultb <- ortho_leebounds(
        leedata = boot_leedata,
        s.hat = s.hat.boot,
        y.hat = y.hat.boot,
        flag_binary = flag_binary,
        flag_helps = flag_helps,
        ortho_d = TRUE
      )
      
      ATE_bb[b, ] <- c(resultb$lower_bound, resultb$upper_bound)
      
    }, error = function(e) {
      # If bootstrap iteration fails, use NA
      ATE_bb[b, ] <- c(NA, NA)
    })
  }
  
  return(ATE_bb)
}