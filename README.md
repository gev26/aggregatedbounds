# Aggregated-Intersection Bounds & Welfare (R)

This repository implements debiased machine-learning estimators for:
- **Aggregated intersection bounds** (e.g., Roy model bounds, Horowitz–Manski–Lee (HML) selection bounds),
- **Optimal welfare / policy value** with **orthogonal** scores,
- **Binary and discrete outcomes** (via envelope/min–max over a finite set).

It follows the “envelope score” idea with cross-fitting and multiplier bootstrap, and supports both **binary** and **multi-category (discrete)** outcomes.

---

## File layout

- `utils.R`  
  Core utilities: data prep conventions, orthogonal score construction, selection/propensity modeling, cross-fitting helpers, bootstrap wrappers, and convenience glue for binary vs discrete outcomes.

- `bound.R`  
  High-level estimators for **Roy model** and **HML bounds** (basic & sharp), plus the **envelope** estimator for “min over t∈T” targets. Includes plug-in and weighted bootstrap inference.

---

## Methods (one paragraph each)

- **Envelope score estimator**  
  Targets parameters of the form \( \psi_0 = \mathbb{E}_X[\min_{t\in T}\, \phi(t,\nu_0(X))] \). We replace each \(\phi\) by an **orthogonal/DML** signal \(\rho(W,t,\xi_0)\) and evaluate it at the **estimated argmin** \(\hat t(X)\) with **cross-fitted** nuisances; asymptotically, this achieves an “oracle” property (first-order robust to misclassification of the minimizer).

- **Orthogonality & cross-fitting**  
  Nuisance pieces (e.g., outcome regressions, selection/propensity models, class probabilities for discrete outcomes) are learned on held-out folds and plugged into orthogonal scores; this yields \( \sqrt{N} \)-consistent, asymptotically normal estimators under standard margin and rate conditions.

- **Bootstrap inference**  
  Weighted (multiplier) bootstrap treats the estimated minimizer as fixed inside each resample and perturbs the second-stage score for valid CIs in finite samples.

---

## Data schema

Your analysis data frame (call it `df`) should include:

| Column          | Type / Values                 | Notes |
|-----------------|-------------------------------|-------|
| `treat`         | {0,1}                         | Treatment or instrument (context-dependent). |
| `outcome`       | numeric or **factor**         | If **discrete**, **store as factor** and pass `Y_levels`. |
| `selection`     | {0,1}                         | For selection/HML bounds (survey response, observed outcome flag). |
| Covariates `X*` | numeric / factors             | Any number of baseline covariates; formulas accept standard R syntax. |

---

## Key functions (API at a glance)

> Function names below mirror what’s in the two files (based on the codebase). If your local names differ slightly, the descriptions still map 1-to-1.

### From `utils.R`

- `estimate_selection(leedata, form = NULL, selection_function_name, variables_for_selection = NULL, names_to_include = c(), treat_name = "treat", yname = "selection", myweights = NULL, ...)`
  - Fits **selection/propensity** models (e.g., logistic) possibly after a **post-lasso** step; returns \(\hat s_d(x)\) or \(\hat \mu(x)\) objects used downstream.
  - For selection problems: ensures the **always-taker share** and conditional selection probabilities are available.

- Cross-fitting & bootstrap helpers  
  - Utilities to split folds, refit nuisances out-of-fold, and run **multiplier bootstrap**.

### From `bound.R`

- **Envelope/Welfare (binary & discrete)**  
  - `main_bb_ortho(sample_size, N_rep = 10, leedata, s.hat, y.hat, flag_binary = TRUE, flag_helps = FALSE, ...)`  
    Orthogonal (DML) bootstrap for average welfare and related functionals. Works for **binary** outcomes and (with flags) for **discrete** outcomes.
  - `main_bb_unified(leedata, N_rep = 500, mode = c("binary","discrete"), ...)`  
    Unified entry point combining orthogonal and non-orthogonal plug-ins; supports **discrete outcomes** when `mode = "discrete"` and `Y_levels` supplied.

- **Horowitz–Manski–Lee (HML) bounds with discrete Y**  
  - Functions to compute **basic** (no-covariate) and **sharp** (covariate-adjusted/envelope) bounds for \( \beta_1 = \mathbb{E}[Y(1)\mid S(1)=S(0)=1] \), including numerators/denominator and sorted bounds.
  - Expect arguments like:
    - Fitted selection \( \hat s_0(0,x), \hat s_0(1,x) \),
    - Treated outcome PMF \( \hat \pi_\beta(x) = P(Y=\beta\mid D=1,S=1,X=x) \) for \(\beta \in \mathcal{T}\),
    - Cross-fitting control, bootstrap reps, and weights.

- **Roy model bounds**  
  - Estimators for the **min over instrument values** (or treatment arms) of conditional joint events, with optional aggregation over \(X\).

> Tip: in **discrete** mode, provide `Y_levels` explicitly (sorted unique outcome support) when `outcome` is a **factor**.

---

## Installation

```r
# R >= 4.2 recommended
# Install suggested dependencies
install.packages(c(
  "glmnet",      # lasso/post-lasso
  "nnet",        # multinomial logit for discrete Y
  "Matrix", "data.table", "mltools"
))

# If you use cross-fitting with ML learners, add what you need, e.g.:
# install.packages(c("ranger", "xgboost", "lightgbm"))
