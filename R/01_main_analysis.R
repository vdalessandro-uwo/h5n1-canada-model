############################################################
# Quantifying Human Risk from H5N1 in Canada:
# multi-species SEIR(D) model with an environmental reservoir
#
# Script 1 of 2: main analysis (Technical Appendix 2)
# R version 4.5.1; packages: deSolve, dplyr, tidyr, tibble, ggplot2,
# purrr, stringr, forcats, coda, scales, patchwork
# Figures and tables are written to the working directory.
#
# PART 1  Setup, parameters and ODE system
# PART 2  Synthetic calibration targets, Bayesian calibration (MCMC),
#         posterior predictive checks, human-to-human transmission
# PART 3  Source attribution, sensitivity analysis, PHSM scenarios,
#         deaths and figures
# PART 4  One-year extension: restocking and seasonal wild-bird importation
# PART 5  Supplementary figures
############################################################

############################################################
# PART 1  Setup, parameters and ODE system
############################################################

## ------------------- Packages ---------------------------
suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(forcats)
})

## ------------------- Globals & helpers ------------------
N_H <- 42e6  # Human population (Canada scale, includes workers)
N_W <- 6000  # High-exposure worker subset
N_G <- N_H - N_W  # General human population (excluding workers)

per100k <- function(x) 100000 * x / N_H  # rate per 100 000 total population
week_index <- function(day_vec) floor(day_vec/7)

# Discrete convolution of incidence with a PMF
conv_incidence <- function(inc, pmf){
  inc <- as.numeric(inc); pmf <- as.numeric(pmf)
  n <- length(inc); m <- length(pmf)
  out <- numeric(n)
  for (t in seq_len(n)){
    kmax <- min(t, m)
    out[t] <- sum(inc[t - (1:kmax) + 1] * pmf[1:kmax])
  }
  out
}

# Reporting delay PMF (gamma; shape 4, scale 1.75: mean 7 days, SD 3.5 days)
delay_pmf <- function(max_days=30, shape=4, scale=7/4){
  days <- 0:max_days
  p <- pgamma(days+1, shape=shape, scale=scale) -
    pgamma(days,   shape=shape, scale=scale)
  p / sum(p)
}

# Infection-to-death delay PMF (lognormal; meanlog = log(18) - 0.5, sdlog = 0.4:
# median ~11 days, mean ~12 days)
death_pmf <- function(max_days=45, meanlog=log(18) - 0.5, sdlog=0.4){
  days <- 0:max_days
  p <- plnorm(days+1, meanlog=meanlog, sdlog=sdlog) -
    plnorm(days,   meanlog=meanlog, sdlog=sdlog)
  p / sum(p)
}

# Aggregate daily series to weekly (sum)
aggregate_weekly <- function(daily, times_days){
  stopifnot(length(daily) == length(times_days))
  tibble(day = times_days, val = as.numeric(daily)) %>%
    mutate(week = week_index(day)) %>%
    group_by(week) %>%
    summarise(val = sum(val), .groups="drop")
}

## ------------------- Time & kinetic params ---------------
times <- 0:104  # 15 weeks (0..104 days)

# Incubation & infectious periods
sigma_H <- 1/3; gamma_H <- 1/5
sigma_W <- 1/3; gamma_W <- 1/5

# IFR for those in the modelled clinical compartments
IFR_H <- 0.10
IFR_W <- 0.10

# Mortality hazards so that IFR = mu/(mu+gamma)
total_exit_H <- 1/5
total_exit_W <- 1/5

mu_H    <- IFR_H * total_exit_H
gamma_H <- (1 - IFR_H) * total_exit_H

mu_W    <- IFR_W * total_exit_W
gamma_W <- (1 - IFR_W) * total_exit_W

# Other species
sigma_WB <- 1/3; gamma_WB <- 1/5
sigma_BP <- 1/3; gamma_BP <- 1/5
sigma_CPL <- 1/3; gamma_CPL <- 1/5
sigma_CPS <- 1/3; gamma_CPS <- 1/5
sigma_X   <- 1/3; gamma_X   <- 1/5

# Environment
k_env_sat  <- 2.0    # saturation
delta_env  <- 0.25   # clearance
sh_H <- 0.02; sh_W <- 0.02
sh_WB <- 0.05; sh_BP <- 0.07; sh_CPL <- 0.20; sh_CPS <- 0.15; sh_X <- 0.03

## ------------------- Baseline transmission params --------
params <- c(
  # Human targets (into general humans)
  beta_HH   = 1e-7,
  beta_WH   = 1e-7,
  beta_WBH  = 1e-8,
  beta_BPH  = 1e-6,
  beta_CPLH = 0,
  beta_CPSH = 0,
  beta_XH   = 2e-6,
  beta_EH   = 2e-6,

  # Worker targets
  beta_WW   = 1e-6,
  beta_HW   = 0.0,
  beta_WBW  = 1e-7,
  beta_BPW  = 1e-6,
  beta_CPLW = 5e-4,
  beta_CPSW = 1.5e-4,
  beta_XW   = 2e-5,
  beta_EW   = 1e-6,

  # Wild birds
  beta_WBWB = 0.3,
  beta_WBBP  = 1e-3,
  beta_WBCPL = 5e-4,
  beta_WBCPS = 7e-4,
  beta_EWB  = 0.0,

  # Backyard poultry
  beta_BPBP = 0.6,
  beta_CPLBP= 0.0,
  beta_CPSBP= 0.0,
  beta_EBP  = 0.0,

  # Commercial poultry (Large)
  beta_CPL_CPL = 0.9,
  beta_CPLCPS  = 0.0,
  beta_BPCPL   = 0.0,
  beta_ECPL    = 0.0,

  # Commercial poultry (Small)
  beta_CPS_CPS = 0.7,
  beta_CPSCPL  = 0.0,
  beta_BPCPS   = 0.0,
  beta_ECPS    = 0.0,

  # Reassortment host X
  beta_XX   = 0.15,
  beta_WBX  = 0.0,
  beta_BPX  = 0.0,
  beta_CPLX = 0.0,
  beta_CPSX = 0.0,
  beta_EX   = 0.0,

  # Kinetics (all populations)
  sigma_H = sigma_H, gamma_H = gamma_H, mu_H = mu_H,
  sigma_W = sigma_W, gamma_W = gamma_W, mu_W = mu_W,
  sigma_WB = sigma_WB, gamma_WB = gamma_WB,
  sigma_BP = sigma_BP, gamma_BP = gamma_BP,
  sigma_CPL = sigma_CPL, gamma_CPL = gamma_CPL,
  sigma_CPS = sigma_CPS, gamma_CPS = gamma_CPS,
  sigma_X = sigma_X, gamma_X = gamma_X,

  # Culling rates
  k_cull_CPL = 0.25,
  k_cull_CPS = 0.25,
  # Restocking of depopulated commercial premises and wild-bird turnover (R -> S);
  # 0 in the 15-week analyses, set in stoch_cfg for the one-year extension
  omega_CPL = 0,
  omega_CPS = 0,
  omega_WB  = 0,

  # Environment
  k_env_sat = k_env_sat,
  delta_env = delta_env,
  sh_H = sh_H, sh_W = sh_W, sh_WB = sh_WB, sh_BP = sh_BP,
  sh_CPL = sh_CPL, sh_CPS = sh_CPS, sh_X = sh_X
)

## ------------------- Initial state (proportions) ----------
init_state <- c(
  S_H = 1 - (1 / N_G), E_H = 0, I_H = 1 / N_G, R_H = 0, D_H = 0,
  S_W = 1, E_W = 0, I_W = 0, R_W = 0, D_W = 0,
  S_WB=0.95,   E_WB = 0.02,   I_WB = 0.02,   R_WB = 0.01,
  S_BP=0.97,   E_BP=0.015, I_BP=0.015, R_BP=0.0,
  S_CPL = 0.99995, E_CPL = 0.00003, I_CPL = 0.00002, R_CPL = 0.0,
  S_CPS = 0.99995, E_CPS = 0.00003, I_CPS = 0.00002, R_CPS = 0.0,
  S_X=0.999,   E_X=0.0005, I_X=0.0005, R_X=0.0,
  Env=0.0
)

## ------------------- ODE system (SEIRD + Env) -------------
multi_species_seir_env <- function(t, y, p){
  with(as.list(c(y, p)), {

    Env_eff <- Env / (1 + k_env_sat * Env)

    fracW_inH <- N_W / N_G
    fracH_inW <- 1

    # Force of infection into general humans
    lambda_H <- beta_HH  * I_H  +
      beta_WH   * fracW_inH * I_W   +
      beta_WBH  * I_WB   +
      beta_BPH  * I_BP   +
      beta_CPLH * I_CPL  +
      beta_CPSH * I_CPS  +
      beta_XH   * I_X    +
      beta_EH   * Env_eff

    # Force of infection into workers
    lambda_W <- beta_WW   * I_W   +
      beta_HW    * fracH_inW * I_H  +
      beta_WBW   * I_WB  +
      beta_BPW   * I_BP  +
      beta_CPLW  * I_CPL +
      beta_CPSW  * I_CPS +
      beta_XW    * I_X   +
      beta_EW    * Env_eff

    # Wild birds
    lambda_WB <- beta_WBWB * I_WB +
      beta_WBBP * I_BP +
      beta_WBCPL* I_CPL +
      beta_WBCPS* I_CPS +
      beta_EWB  * Env_eff

    # Backyard poultry
    lambda_BP <- beta_BPBP * I_BP +
      beta_WBBP * I_WB +
      beta_CPLBP* I_CPL +
      beta_CPSBP* I_CPS +
      beta_EBP  * Env_eff

    # Commercial poultry (Large)
    lambda_CPL <- beta_CPL_CPL * I_CPL +
      beta_WBCPL   * I_WB +
      beta_CPLCPS  * I_CPS +
      beta_BPCPL   * I_BP +
      beta_ECPL    * Env_eff

    # Commercial poultry (Small)
    lambda_CPS <- beta_CPS_CPS * I_CPS +
      beta_WBCPS   * I_WB +
      beta_CPSCPL  * I_CPL +
      beta_BPCPS   * I_BP +
      beta_ECPS    * Env_eff

    # Reassortment host X
    lambda_X <- beta_XX * I_X +
      beta_WBX * I_WB +
      beta_BPX * I_BP  +
      beta_CPLX* I_CPL +
      beta_CPSX* I_CPS +
      beta_EX  * Env_eff

    # Humans (SEIRD)
    dS_H <- - lambda_H * S_H
    dE_H <-   lambda_H * S_H - sigma_H * E_H
    dI_H <-   sigma_H * E_H - (gamma_H + mu_H) * I_H
    dR_H <-   gamma_H * I_H
    dD_H <-   mu_H * I_H

    # Workers (SEIRD)
    dS_W <- - lambda_W * S_W
    dE_W <-   lambda_W * S_W - sigma_W * E_W
    dI_W <-   sigma_W * E_W -(gamma_W + mu_W) * I_W
    dR_W <-   gamma_W * I_W
    dD_W <-   mu_W    * I_W

    # Wild birds
    dS_WB <- - lambda_WB * S_WB + omega_WB * R_WB
    dE_WB <-   lambda_WB * S_WB - sigma_WB * E_WB
    dI_WB <-   sigma_WB * E_WB - gamma_WB * I_WB
    dR_WB <-   gamma_WB * I_WB - omega_WB * R_WB

    # Backyard poultry
    dS_BP <- - lambda_BP * S_BP
    dE_BP <-   lambda_BP * S_BP - sigma_BP * E_BP
    dI_BP <-   sigma_BP * E_BP - gamma_BP * I_BP
    dR_BP <-   gamma_BP * I_BP

    # Commercial poultry (Large): removal/depopulation of infected flocks/premises
    dS_CPL <- - lambda_CPL * S_CPL + omega_CPL * R_CPL
    dE_CPL <-   lambda_CPL * S_CPL - sigma_CPL * E_CPL
    dI_CPL <-   sigma_CPL * E_CPL - k_cull_CPL * I_CPL
    dR_CPL <-   k_cull_CPL * I_CPL - omega_CPL * R_CPL

    # Commercial poultry (Small): removal/depopulation of infected flocks/premises
    dS_CPS <- - lambda_CPS * S_CPS + omega_CPS * R_CPS
    dE_CPS <-   lambda_CPS * S_CPS - sigma_CPS * E_CPS
    dI_CPS <-   sigma_CPS * E_CPS - k_cull_CPS * I_CPS
    dR_CPS <-   k_cull_CPS * I_CPS - omega_CPS * R_CPS

    # Reassortment host X
    dS_X <- - lambda_X * S_X
    dE_X <-   lambda_X * S_X - sigma_X * E_X
    dI_X <-   sigma_X * E_X - gamma_X * I_X
    dR_X <-   gamma_X * I_X

    # Environment
    shedding <- sh_H  * I_H + sh_W  * I_W +
      sh_WB * I_WB + sh_BP * I_BP +
      sh_CPL* I_CPL + sh_CPS* I_CPS +
      sh_X  * I_X
    dEnv <- shedding - delta_env * Env

    list(c(
      dS_H, dE_H, dI_H, dR_H, dD_H,
      dS_W, dE_W, dI_W, dR_W, dD_W,
      dS_WB, dE_WB, dI_WB, dR_WB,
      dS_BP, dE_BP, dI_BP, dR_BP,
      dS_CPL, dE_CPL, dI_CPL, dR_CPL,
      dS_CPS, dE_CPS, dI_CPS, dR_CPS,
      dS_X, dE_X, dI_X, dR_X,
      dEnv
    ))
  })
}

############################################################
# PART 2  Synthetic targets, Bayesian calibration, posterior predictive
#         checks and human-to-human transmission
############################################################

## ------------------- Deterministic run -------------------
ode_out <- as_tibble(as.data.frame(
  ode(y=init_state, times=times, func=multi_species_seir_env, parms=params)
))

# Long format with Population/State
seir_long <- ode_out %>%
  pivot_longer(cols = -time, names_to = "compartment", values_to = "value") %>%
  mutate(
    Population = case_when(
      str_ends(compartment, "_H")   ~ "Human",
      str_ends(compartment, "_W")   ~ "Worker",
      str_ends(compartment, "_WB")  ~ "Wild Bird",
      str_ends(compartment, "_BP")  ~ "Backyard Poultry",
      str_ends(compartment, "_CPL") ~ "Commercial (Large)",
      str_ends(compartment, "_CPS") ~ "Commercial (Small)",
      str_ends(compartment, "_X")   ~ "Reassortment Host",
      compartment == "Env"          ~ "Environment",
      TRUE ~ NA_character_
    ),
    State = case_when(
      str_starts(compartment, "S_") ~ "S",
      str_starts(compartment, "E_") ~ "E",
      str_starts(compartment, "I_") ~ "I",
      str_starts(compartment, "R_") ~ "R",
      str_starts(compartment, "D_") ~ "D",
      compartment == "Env"          ~ "Env",
      TRUE ~ NA_character_
    )
  ) %>% filter(!is.na(Population), !is.na(State))

## ------------------- Synthetic weekly counts -------------
DO_DELAY <- TRUE
rho_true <- 0.6
phi_true <- 50

times_days <- ode_out$time
obs_weeks <- 0:14

# Incidence (per-capita) from ODE outputs
incidence_H_daily <- params["sigma_H"] * ode_out$E_H
incidence_W_daily <- params["sigma_W"] * ode_out$E_W

# Apply reporting delay
mu_H_prop_daily <- if (DO_DELAY) conv_incidence(incidence_H_daily, delay_pmf()) else incidence_H_daily
mu_W_prop_daily <- if (DO_DELAY) conv_incidence(incidence_W_daily, delay_pmf()) else incidence_W_daily

# Aggregate to weekly
mu_H_weekly <- aggregate_weekly(mu_H_prop_daily, times_days) %>%
  filter(week %in% obs_weeks) %>% pull(val)

mu_W_weekly <- aggregate_weekly(mu_W_prop_daily, times_days) %>%
  filter(week %in% obs_weeks) %>% pull(val)

# Convert to weekly expected counts
mu_H_counts_weekly <- pmax(mu_H_weekly * N_G * rho_true, 1e-9)
mu_W_counts_weekly <- pmax(mu_W_weekly * N_W * rho_true, 1e-9)
mu_counts_weekly <- mu_H_counts_weekly + mu_W_counts_weekly

set.seed(123)
obs_counts_weekly <- rnbinom(length(mu_counts_weekly), size = phi_true, mu = mu_counts_weekly)

cat("\n=== Synthetic weekly counts (Canada scale) ===\n")
cat("Total observed counts (weeks 0..", max(obs_weeks), "): ",
    format(sum(obs_counts_weekly), big.mark=","), " (",
    round(per100k(sum(obs_counts_weekly)), 3), " per 100k)\n", sep = "")

## ------------------- Bayesian fitting --------------------
par_names <- c(
  "beta_HH","beta_WH","beta_WBH","beta_BPH","beta_EH","beta_XH",
  "beta_WW","beta_WBW","beta_BPW","beta_CPLW","beta_CPSW","beta_EW","beta_XW"
)

# Priors on log10(beta)
prior_mean_log10 <- c(
  beta_HH   = -7, beta_WH   = -5.5, beta_WBH  = -5.5, beta_BPH  = -6.0,
  beta_EH   = -6.0, beta_XH   = -5.7, beta_WW   = -6.5, beta_WBW  = -4.0,
  beta_BPW  = -4.5,beta_CPLW = -3.3, beta_CPSW = -3.8, beta_EW   = -4.7,
  beta_XW   = -4.7
)
prior_sd_log10 <- c(
  beta_HH=0.7, beta_WH=0.6, beta_WBH=0.6, beta_BPH=0.6,
  beta_EH=0.6, beta_XH=0.6,
  beta_WW=0.7, beta_WBW=0.7, beta_BPW=0.7,
  beta_CPLW=0.7, beta_CPSW=0.7, beta_EW=0.7, beta_XW=0.7
)

logprior_betas_core <- function(theta_log10){
  sum(dnorm(theta_log10,
            mean=prior_mean_log10[names(theta_log10)],
            sd  =prior_sd_log10[names(theta_log10)], log=TRUE))
}

logprior_order <- function(theta_log10){
  0
}

logprior_betas <- function(theta_log10){
  lp <- logprior_betas_core(theta_log10) + logprior_order(theta_log10)
  if (!is.finite(lp)) -1e99 else lp
}

logprior_logphi <- function(log_phi){
  lp <- dnorm(log_phi, mean=log(50), sd=0.5, log=TRUE)
  if (!is.finite(lp)) -1e99 else lp
}

logprior_logit_rho <- function(logit_rho){
  rho <- plogis(logit_rho)
  lp <- dbeta(rho, 5, 5, log=TRUE) + log(rho*(1-rho))
  if (!is.finite(lp)) -1e99 else lp
}

loglik_nb <- function(y_counts, mu, phi){
  if (!is.finite(phi) || phi <= 0) return(-1e99)
  mu <- pmax(mu, 1e-12)
  if (any(!is.finite(mu))) return(-1e99)
  ll <- sum(dnbinom(y_counts, size=phi, mu=mu, log=TRUE))
  if (!is.finite(ll)) -1e99 else ll
}

simulate_incidence <- function(par_list, times, init_state){
  sim <- ode(y=init_state, times=times,
             func=multi_species_seir_env, parms=par_list)
  df  <- as_tibble(as.data.frame(sim))

  inc_H <- par_list["sigma_H"] * df$E_H
  inc_W <- par_list["sigma_W"] * df$E_W

  list(df=df, inc_H=inc_H, inc_W=inc_W)
}

posterior_counts_weekly <- function(theta_log10, log_phi, logit_rho,
                                    base_params, times, init_state,
                                    obs_weeks, y_weekly, N_H,
                                    DO_DELAY=TRUE){

  par_try <- base_params
  for (nm in names(theta_log10)) par_try[nm] <- 10^(theta_log10[nm])

  sim <- try(simulate_incidence(par_try, times, init_state), silent=TRUE)
  if (inherits(sim, "try-error")) return(-1e99)

  muH_prop_daily <- if (DO_DELAY) conv_incidence(sim$inc_H, delay_pmf()) else sim$inc_H
  muW_prop_daily <- if (DO_DELAY) conv_incidence(sim$inc_W, delay_pmf()) else sim$inc_W

  muH_weekly <- aggregate_weekly(muH_prop_daily, times) %>%
    filter(week %in% obs_weeks) %>% pull(val)

  muW_weekly <- aggregate_weekly(muW_prop_daily, times) %>%
    filter(week %in% obs_weeks) %>% pull(val)

  if (length(muH_weekly) != length(y_weekly) || length(muW_weekly) != length(y_weekly)) return(-1e99)

  rho <- plogis(logit_rho)
  mu <- (muH_weekly * (N_H - N_W) + muW_weekly * N_W) * rho

  ll <- loglik_nb(y_weekly, mu, phi=exp(log_phi))
  lp <- logprior_betas(theta_log10) +
    logprior_logphi(log_phi) +
    logprior_logit_rho(logit_rho)

  ll + lp
}

# --- MCMC settings ---
n_iter   <- 30000
burn_in  <- 10000
thin     <- 10
n_chains <- 4

keep_idx <- seq(burn_in + 1, n_iter, by = thin)
n_keep   <- length(keep_idx)

theta0_log10 <- log10(as.numeric(params[par_names]))
names(theta0_log10) <- par_names

rho_start <- if (exists("rho_true")) rho_true else 0.5

# Starting values
logphi0    <- log(50)
logit_rho0 <- qlogis(rho_start)

# Proposal settings
base_prop_sd_theta     <- rep(0.15, length(par_names))
names(base_prop_sd_theta) <- par_names
base_prop_sd_logphi    <- 0.15
base_prop_sd_logit_rho <- 0.20

window_size <- 200L

# Function to run one Metropolis-Hastings chain
run_mcmc_chain <- function(chain_id, seed = 2468) {

  set.seed(seed + chain_id)

  # Dispersed starting values across chains
  theta_curr <- theta0_log10 + rnorm(length(theta0_log10), mean = 0, sd = 0.20)
  names(theta_curr) <- par_names

  logphi_curr    <- logphi0 + rnorm(1, mean = 0, sd = 0.20)
  logit_rho_curr <- logit_rho0 + rnorm(1, mean = 0, sd = 0.30)

  # If dispersed starting values produce invalid posterior, fall back to central start
  post_curr <- posterior_counts_weekly(
    theta_curr, logphi_curr, logit_rho_curr,
    params, times, init_state,
    obs_weeks, obs_counts_weekly, N_H, DO_DELAY
  )

  if (!is.finite(post_curr)) {
    theta_curr     <- theta0_log10
    logphi_curr    <- logphi0
    logit_rho_curr <- logit_rho0

    post_curr <- posterior_counts_weekly(
      theta_curr, logphi_curr, logit_rho_curr,
      params, times, init_state,
      obs_weeks, obs_counts_weekly, N_H, DO_DELAY
    )
  }

  prop_sd_theta     <- base_prop_sd_theta
  prop_sd_logphi    <- base_prop_sd_logphi
  prop_sd_logit_rho <- base_prop_sd_logit_rho

  # Draw storage on natural scale for downstream analysis
  draws_theta <- matrix(NA_real_, nrow = n_keep, ncol = length(par_names))
  colnames(draws_theta) <- par_names

  phi_draws <- numeric(n_keep)
  rho_draws <- numeric(n_keep)

  # Draw storage on transformed scale for MCMC diagnostics
  diag_draws <- matrix(
    NA_real_,
    nrow = n_keep,
    ncol = length(par_names) + 2
  )
  colnames(diag_draws) <- c(par_names, "log_phi", "logit_rho")

  acc <- 0L
  row_i <- 0L
  accept_window <- 0L

  for (it in 1:n_iter) {

    theta_prop     <- theta_curr + rnorm(length(theta_curr), 0, prop_sd_theta)
    names(theta_prop) <- par_names

    logphi_prop    <- logphi_curr + rnorm(1, 0, prop_sd_logphi)
    logit_rho_prop <- logit_rho_curr + rnorm(1, 0, prop_sd_logit_rho)

    post_prop <- posterior_counts_weekly(
      theta_prop, logphi_prop, logit_rho_prop,
      params, times, init_state,
      obs_weeks, obs_counts_weekly, N_H, DO_DELAY
    )

    log_alpha <- post_prop - post_curr

    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      theta_curr     <- theta_prop
      logphi_curr    <- logphi_prop
      logit_rho_curr <- logit_rho_prop
      post_curr      <- post_prop

      acc <- acc + 1L
      accept_window <- accept_window + 1L
    }

    # Adaptive proposal scaling during burn-in only
    if (it <= burn_in && it %% window_size == 0) {
      ar <- accept_window / window_size

      scale <- if (ar < 0.18) {
        0.8
      } else if (ar > 0.40) {
        1.2
      } else {
        1.0
      }

      prop_sd_theta     <- prop_sd_theta * scale
      prop_sd_logphi    <- prop_sd_logphi * scale
      prop_sd_logit_rho <- prop_sd_logit_rho * scale

      accept_window <- 0L
    }

    # Save thinned post-burn-in draws
    if (it %in% keep_idx) {
      row_i <- row_i + 1L

      draws_theta[row_i, ] <- 10^(theta_curr)
      phi_draws[row_i]     <- exp(logphi_curr)
      rho_draws[row_i]     <- plogis(logit_rho_curr)

      diag_draws[row_i, ]  <- c(theta_curr, logphi_curr, logit_rho_curr)
    }
  }

  draws_df <- as_tibble(draws_theta) %>%
    mutate(
      phi = phi_draws,
      rho = rho_draws,
      chain = chain_id,
      draw = seq_len(n_keep),
      iteration = keep_idx
    )

  list(
    draws = draws_df,
    diag_draws = diag_draws,
    acceptance = acc / n_iter
  )
}

# Run chains
chain_results <- vector("list", n_chains)

for (ch in seq_len(n_chains)) {
  cat("\nRunning MCMC chain ", ch, " of ", n_chains, "...\n", sep = "")
  chain_results[[ch]] <- run_mcmc_chain(chain_id = ch, seed = 2468)
  cat(
    "Chain ", ch, " acceptance: ",
    sprintf("%.1f%%", 100 * chain_results[[ch]]$acceptance),
    "\n",
    sep = ""
  )
}

# Combine posterior draws
posterior_draws_with_meta <- bind_rows(lapply(chain_results, `[[`, "draws"))

# Use this object for all downstream analyses
posterior_draws <- posterior_draws_with_meta %>%
  select(all_of(par_names), phi, rho)

cat("\nMCMC retained posterior samples:\n")
cat("Chains: ", n_chains, "\n", sep = "")
cat("Iterations per chain: ", n_iter, "\n", sep = "")
cat("Burn-in per chain: ", burn_in, "\n", sep = "")
cat("Thinning interval: ", thin, "\n", sep = "")
cat("Retained samples per chain: ", n_keep, "\n", sep = "")
cat("Total retained posterior samples: ", nrow(posterior_draws), "\n", sep = "")

cat("\nMCMC acceptance rates:\n")
print(
  tibble(
    chain = seq_len(n_chains),
    acceptance_percent = 100 * sapply(chain_results, `[[`, "acceptance")
  )
)

# Posterior summaries
post_sum <- posterior_draws %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  summarise(
    mean   = mean(value),
    median = median(value),
    p2.5   = quantile(value, 0.025),
    p25    = quantile(value, 0.25),
    p75    = quantile(value, 0.75),
    p97.5  = quantile(value, 0.975),
    .groups = "drop"
  )

cat("\nPosterior summaries (Canada scale):\n")
print(post_sum)

# Posterior-median parameter set (deterministic final size, weekly attribution
# figure and one-year extension)
par_med <- params
for (nm in par_names) par_med[nm] <- median(posterior_draws[[nm]])

# Convergence diagnostics (Gelman-Rubin R-hat, effective sample size)
if (requireNamespace("coda", quietly = TRUE)) {

  mcmc_list <- coda::mcmc.list(
    lapply(chain_results, function(x) {
      coda::mcmc(x$diag_draws)
    })
  )

  gelman_diag <- coda::gelman.diag(
    mcmc_list,
    autoburnin = FALSE,
    multivariate = FALSE
  )

  ess <- coda::effectiveSize(mcmc_list)

  cat("\nGelman-Rubin diagnostics / R-hat:\n")
  print(gelman_diag)

  cat("\nEffective sample sizes:\n")
  print(ess)

} else {
  message(
    "\nPackage 'coda' is not installed. ",
    "Install it with install.packages('coda') to calculate R-hat and effective sample sizes."
  )
}

# Technical Appendix Figure 2. MCMC trace plots
trace_df <- posterior_draws_with_meta %>%
  select(chain, iteration, all_of(par_names), phi, rho) %>%
  pivot_longer(
    cols = -c(chain, iteration),
    names_to = "parameter",
    values_to = "value"
  )

p_trace <- ggplot(trace_df, aes(x = iteration, y = value, colour = factor(chain))) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(
    title = "MCMC trace plots",
    x = "Iteration",
    y = "Parameter value",
    colour = "Chain"
  ) +
  theme_minimal()

print(p_trace)
ggsave("MCMC_trace_plots.png", p_trace, width = 14, height = 10, dpi = 300)

## ------------------- Posterior predictive counts ----------
## Separate reported cases from total modelled infections

# Helper to make sure all requested weeks are returned
get_weekly_vals <- function(x_daily, time, obs_weeks) {
  wk_tbl <- aggregate_weekly(x_daily, time)

  out <- tibble(week = obs_weeks) %>%
    left_join(wk_tbl, by = "week") %>%
    arrange(week)

  out$val[is.na(out$val)] <- 0
  out$val
}

# Unified H + W function
simulate_weekly_human_counts <- function(par, times, init_state, obs_weeks,
                                         N_H, N_W,
                                         rho = 1,
                                         outcome = c("infections", "reported"),
                                         DO_DELAY = TRUE) {

  outcome <- match.arg(outcome)
  N_G <- N_H - N_W

  sim <- try(
    ode(y = init_state, times = times,
        func = multi_species_seir_env, parms = par),
    silent = TRUE
  )

  if (inherits(sim, "try-error")) {
    return(tibble(
      week = obs_weeks,
      H = NA_real_,
      W = NA_real_,
      total = NA_real_
    ))
  }

  df <- as_tibble(as.data.frame(sim))

  if (!all(c("time", "E_H", "E_W") %in% names(df))) {
    return(tibble(
      week = obs_weeks,
      H = NA_real_,
      W = NA_real_,
      total = NA_real_
    ))
  }

  # Incident modelled infections before reporting adjustment
  incH_prop_daily <- par["sigma_H"] * df$E_H
  incW_prop_daily <- par["sigma_W"] * df$E_W

  # Reported cases include reporting delay and reporting probability
  if (outcome == "reported") {
    if (DO_DELAY) {
      incH_prop_daily <- conv_incidence(incH_prop_daily, delay_pmf())
      incW_prop_daily <- conv_incidence(incW_prop_daily, delay_pmf())
    }
    multiplier <- rho
  } else {
    # Total modelled infections should not be multiplied by rho
    multiplier <- 1
  }

  H_week <- get_weekly_vals(incH_prop_daily, df$time, obs_weeks) * N_G * multiplier
  W_week <- get_weekly_vals(incW_prop_daily, df$time, obs_weeks) * N_W * multiplier

  tibble(
    week = obs_weeks,
    H = H_week,
    W = W_week,
    total = H_week + W_week
  )
}

# --- Posterior predictive loop ---
K <- min(400, nrow(posterior_draws))
set.seed(765)
idx <- sample(seq_len(nrow(posterior_draws)), K)

pred_mat_reported   <- matrix(NA_real_, nrow = K, ncol = length(obs_weeks))
pred_mat_infections <- matrix(NA_real_, nrow = K, ncol = length(obs_weeks))

for (i in seq_along(idx)) {
  ii <- idx[i]

  par_i <- params
  for (nm in par_names) {
    par_i[nm] <- posterior_draws[[nm]][ii]
  }

  # Expected reported cases: includes reporting delay and rho
  wk_reported <- simulate_weekly_human_counts(
    par = par_i,
    times = times,
    init_state = init_state,
    obs_weeks = obs_weeks,
    N_H = N_H,
    N_W = N_W,
    rho = posterior_draws$rho[ii],
    outcome = "reported",
    DO_DELAY = DO_DELAY
  )

  pred_mat_reported[i, ] <- wk_reported$total

  # Total modelled infections: no reporting delay and no rho adjustment
  wk_infections <- simulate_weekly_human_counts(
    par = par_i,
    times = times,
    init_state = init_state,
    obs_weeks = obs_weeks,
    N_H = N_H,
    N_W = N_W,
    rho = 1,
    outcome = "infections",
    DO_DELAY = FALSE
  )

  pred_mat_infections[i, ] <- wk_infections$total
}

# Data for Figure 2: posterior predictive reported cases vs synthetic targets
pred_df <- tibble(
  week = obs_weeks,
  med  = apply(pred_mat_reported, 2, median),
  lo50 = apply(pred_mat_reported, 2, quantile, 0.25),
  hi50 = apply(pred_mat_reported, 2, quantile, 0.75),
  lo95 = apply(pred_mat_reported, 2, quantile, 0.025),
  hi95 = apply(pred_mat_reported, 2, quantile, 0.975)
) %>%
  mutate(med_per100k = per100k(med))

# Weekly modelled infections
pred_df_infections <- tibble(
  week = obs_weeks,
  med  = apply(pred_mat_infections, 2, median),
  lo50 = apply(pred_mat_infections, 2, quantile, 0.25),
  hi50 = apply(pred_mat_infections, 2, quantile, 0.75),
  lo95 = apply(pred_mat_infections, 2, quantile, 0.025),
  hi95 = apply(pred_mat_infections, 2, quantile, 0.975)
) %>%
  mutate(med_per100k = per100k(med))

# Reported case summaries
peak_week_reported <- pred_df$week[which.max(pred_df$med)]
peak_reported_med  <- max(pred_df$med)

cum_reported_draws  <- rowSums(pred_mat_reported)
cum_infection_draws <- rowSums(pred_mat_infections)

cat("\nPeak weekly reported cases (median predictive): ",
    format(round(peak_reported_med), big.mark = ","),
    " (", round(per100k(peak_reported_med), 3),
    " per 100k) on week ", peak_week_reported, "\n", sep = "")

cat("Cumulative reported cases (weeks 0..", max(obs_weeks), "): median ",
    format(round(median(cum_reported_draws)), big.mark = ","),
    " [", format(round(quantile(cum_reported_draws, 0.025)), big.mark = ","),
    ", ", format(round(quantile(cum_reported_draws, 0.975)), big.mark = ","),
    "] ~ ", round(per100k(median(cum_reported_draws)), 3), " per 100k\n", sep = "")

cat("Cumulative modelled infections (weeks 0..", max(obs_weeks), "): median ",
    format(round(median(cum_infection_draws)), big.mark = ","),
    " [", format(round(quantile(cum_infection_draws, 0.025)), big.mark = ","),
    ", ", format(round(quantile(cum_infection_draws, 0.975)), big.mark = ","),
    "] ~ ", round(per100k(median(cum_infection_draws)), 3), " per 100k\n", sep = "")

## ------------------- Sustained H2H probability -----------
h2h_thresh_main <- 0.05
h2h_thresh_alt  <- 0.02

sustained_from_draw <- function(draw_row, thresh) {
  par_try <- params
  for (nm in par_names) par_try[nm] <- posterior_draws[[nm]][draw_row]
  par_try[c("beta_WH","beta_WBH","beta_BPH","beta_CPLH","beta_CPSH","beta_EH","beta_XH")] <- 0
  init_k <- c(S_H=0.999, E_H=0, I_H=0.001, R_H=0, D_H=0,
              S_W=1, E_W=0, I_W=0, R_W=0, D_W=0,
              S_WB=1, E_WB=0, I_WB=0, R_WB=0,
              S_BP=1, E_BP=0, I_BP=0, R_BP=0,
              S_CPL=1,E_CPL=0,I_CPL=0,R_CPL=0,
              S_CPS=1,E_CPS=0,I_CPS=0,R_CPS=0,
              S_X=1, E_X=0, I_X=0, R_X=0,
              Env=0)
  df <- as_tibble(as.data.frame(
    ode(y=init_k, times=0:180,
        func=multi_species_seir_env, parms=par_try)
  ))
  any(df$I_H > thresh & par_try["beta_HH"] > 0)
}

set.seed(101)
idx_check <- sample(seq_len(nrow(posterior_draws)), min(1000, nrow(posterior_draws)))
prob_sustained_main <- mean(sapply(idx_check, sustained_from_draw, thresh = h2h_thresh_main))
prob_sustained_alt  <- mean(sapply(idx_check, sustained_from_draw, thresh = h2h_thresh_alt))

# Reproduction numbers; infectious period = 1 / (gamma + mu)
R0_HH_post <- quantile(posterior_draws$beta_HH / (params["gamma_H"] + params["mu_H"]), c(.5, .025, .975))
R0_WW_post <- quantile(posterior_draws$beta_WW / (params["gamma_W"] + params["mu_W"]), c(.5, .025, .975))
cat(sprintf("Posterior R0_WW: median %.2e [%.2e, %.2e]\n", R0_WW_post[1], R0_WW_post[2], R0_WW_post[3]))
cat(sprintf("Posterior R0_HH: median %.2e [%.2e, %.2e]\n",
            R0_HH_post[1], R0_HH_post[2], R0_HH_post[3]))
cat("\nPosterior Prob(sustained H2H): ", round(prob_sustained_main,5),
    " (threshold I_H>", h2h_thresh_main, "); sensitivity @ ",
    h2h_thresh_alt, ": ", round(prob_sustained_alt,5), "\n", sep="")

############################################################
# PART 3  Source attribution, sensitivity analysis, PHSM scenarios,
#         deaths and figures
############################################################

## ------------------- Attribution (posterior-weighted) ----
## Outcome = modelled human infections, H + W, not reported cases

sources <- c(
  "Human/Worker→Human/Worker",
  "WildBird→Human/Worker",
  "Backyard→Human/Worker",
  "CPL→Worker",
  "CPS→Worker",
  "Environment→Human/Worker",
  "X→Human/Worker"
)

get_par <- function(par, nm, default = 0) {
  if (nm %in% names(par)) return(as.numeric(par[nm]))
  default
}

attrib_weekly_from_par <- function(par) {

  if (!exists("N_G")) N_G <- N_H - N_W

  sim <- simulate_incidence(par, times, init_state)
  df  <- sim$df

  Env_eff <- df$Env / (1 + par["k_env_sat"] * df$Env)

  ## ---------------- Source-specific forces into general humans ----------------

  lambda_HH_human <- get_par(par, "beta_HH") * df$I_H
  lambda_WH_human <- get_par(par, "beta_WH") * (N_W / N_G) * df$I_W

  lambda_WB_human  <- get_par(par, "beta_WBH")  * df$I_WB
  lambda_BP_human  <- get_par(par, "beta_BPH")  * df$I_BP
  lambda_X_human   <- get_par(par, "beta_XH")   * df$I_X
  lambda_E_human   <- get_par(par, "beta_EH")   * Env_eff

  lambda_total_H <- lambda_HH_human + lambda_WH_human +
    lambda_WB_human + lambda_BP_human + lambda_X_human + lambda_E_human

  ## ---------------- Source-specific forces into workers ----------------

  lambda_HW_worker <- get_par(par, "beta_HW") * df$I_H
  lambda_WW_worker <- get_par(par, "beta_WW") * df$I_W

  lambda_WB_worker  <- get_par(par, "beta_WBW")  * df$I_WB
  lambda_BP_worker  <- get_par(par, "beta_BPW")  * df$I_BP
  lambda_CPL_worker <- get_par(par, "beta_CPLW") * df$I_CPL
  lambda_CPS_worker <- get_par(par, "beta_CPSW") * df$I_CPS
  lambda_X_worker   <- get_par(par, "beta_XW")   * df$I_X
  lambda_E_worker   <- get_par(par, "beta_EW")   * Env_eff

  lambda_total_W <- lambda_HW_worker + lambda_WW_worker +
    lambda_WB_worker + lambda_BP_worker + lambda_CPL_worker +
    lambda_CPS_worker + lambda_X_worker + lambda_E_worker

  eps <- 1e-12

  ## ---------------- Daily modelled infections, not reported cases ----------------
  incH_counts_daily <- par["sigma_H"] * df$E_H * N_G
  incW_counts_daily <- par["sigma_W"] * df$E_W * N_W

  ## Attribution weights for general humans
  w_H_human <- (lambda_HH_human + lambda_WH_human) / pmax(lambda_total_H, eps)
  w_H_WB    <- lambda_WB_human  / pmax(lambda_total_H, eps)
  w_H_BP    <- lambda_BP_human  / pmax(lambda_total_H, eps)
  w_H_E     <- lambda_E_human   / pmax(lambda_total_H, eps)
  w_H_X     <- lambda_X_human   / pmax(lambda_total_H, eps)

  ## Attribution weights for workers
  w_W_human <- (lambda_HW_worker + lambda_WW_worker) / pmax(lambda_total_W, eps)
  w_W_WB    <- lambda_WB_worker  / pmax(lambda_total_W, eps)
  w_W_BP    <- lambda_BP_worker  / pmax(lambda_total_W, eps)
  w_W_CPL   <- lambda_CPL_worker / pmax(lambda_total_W, eps)
  w_W_CPS   <- lambda_CPS_worker / pmax(lambda_total_W, eps)
  w_W_E     <- lambda_E_worker   / pmax(lambda_total_W, eps)
  w_W_X     <- lambda_X_worker   / pmax(lambda_total_W, eps)

  attrib_daily <- tibble(
    day = df$time,
    `Human/Worker→Human/Worker` =
      incH_counts_daily * w_H_human + incW_counts_daily * w_W_human,
    `WildBird→Human/Worker` =
      incH_counts_daily * w_H_WB + incW_counts_daily * w_W_WB,
    `Backyard→Human/Worker` =
      incH_counts_daily * w_H_BP + incW_counts_daily * w_W_BP,
    `CPL→Worker` =
      incW_counts_daily * w_W_CPL,
    `CPS→Worker` =
      incW_counts_daily * w_W_CPS,
    `Environment→Human/Worker` =
      incH_counts_daily * w_H_E + incW_counts_daily * w_W_E,
    `X→Human/Worker` =
      incH_counts_daily * w_H_X + incW_counts_daily * w_W_X
  )

  attrib_week <- attrib_daily %>%
    pivot_longer(-day, names_to = "source", values_to = "counts") %>%
    mutate(week = week_index(day)) %>%
    group_by(week, source) %>%
    summarise(counts = sum(counts, na.rm = TRUE), .groups = "drop") %>%
    filter(week %in% obs_weeks)

  attrib_week
}

## Get median attribution for plotting
attrib_week_med <- attrib_weekly_from_par(par_med)

## ------------------- Deterministic baseline size & peak ----
if (!exists("N_G")) N_G <- N_H - N_W

times_det <- 0:365

ode_det <- as_tibble(as.data.frame(
  ode(y = init_state,
      times = times_det,
      func  = multi_species_seir_env,
      parms = par_med)
))

S0_H <- init_state["S_H"]
S0_W <- init_state["S_W"]

final_S_H <- tail(ode_det$S_H, 1)
final_S_W <- tail(ode_det$S_W, 1)

final_size_prop_H <- S0_H - final_S_H
final_size_prop_W <- S0_W - final_S_W

final_size_n_H <- final_size_prop_H * N_G
final_size_n_W <- final_size_prop_W * N_W

final_size_n_total <- final_size_n_H + final_size_n_W
final_size_per100k <- per100k(final_size_n_total)

peak_IH <- max(ode_det$I_H)
peak_IW <- max(ode_det$I_W)

peak_day_H <- ode_det$time[which.max(ode_det$I_H)]
peak_day_W <- ode_det$time[which.max(ode_det$I_W)]

peak_prev_pct_H <- 100 * peak_IH
peak_prev_pct_W <- 100 * peak_IW

peak_prev_n_H <- peak_IH * N_G
peak_prev_n_W <- peak_IW * N_W
peak_prev_n_total <- peak_prev_n_H + peak_prev_n_W

cat("\n=== Deterministic baseline (posterior median params) ===\n")

cat("Final epidemic size (general humans): ",
    format(round(final_size_n_H), big.mark=","), " cases (",
    round(100 * final_size_prop_H, 3), "% of general humans)\n", sep="")

cat("Final epidemic size (workers): ",
    format(round(final_size_n_W), big.mark=","), " cases (",
    round(100 * final_size_prop_W, 3), "% of workers)\n", sep="")

cat("Final epidemic size (Canada total): ",
    format(round(final_size_n_total), big.mark=","), " cases (",
    round(final_size_per100k, 3), " per 100,000 of Canada)\n", sep="")

cat("Peak concurrently infectious (general humans): ",
    round(peak_prev_pct_H, 3), "% on day ", peak_day_H,
    " (~", format(round(peak_prev_n_H), big.mark=","), ")\n", sep="")

cat("Peak concurrently infectious (workers): ",
    round(peak_prev_pct_W, 3), "% on day ", peak_day_W,
    " (~", format(round(peak_prev_n_W), big.mark=","), ")\n", sep="")

cat("Peak concurrently infectious (Canada total, H+W): ~",
    format(round(peak_prev_n_total), big.mark=","), "\n", sep="")

# Cumulative attribution with CrI
Kattr <- min(300, nrow(posterior_draws)); set.seed(909)
idx_attr <- sample(seq_len(nrow(posterior_draws)), Kattr)
cum_attr_mat <- matrix(0, nrow = Kattr, ncol = length(sources))
colnames(cum_attr_mat) <- sources

for (i in seq_along(idx_attr)) {
  ii <- idx_attr[i]
  par_i <- params
  for (nm in par_names) par_i[nm] <- posterior_draws[[nm]][ii]
  aw <- attrib_weekly_from_par(par_i)
  vec <- setNames(numeric(length(sources)), sources)
  cum_attr <- aw %>%
    group_by(source) %>%
    summarise(cases = sum(counts), .groups = "drop")
  vec[cum_attr$source] <- cum_attr$cases
  cum_attr_mat[i, ] <- vec
}

attr_sum <- apply(
  cum_attr_mat, 2,
  function(x) c(
    med = median(x, na.rm = TRUE),
    lo  = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    hi  = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
)

attr_df <- t(attr_sum) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("source") %>%
  as_tibble() %>%
  mutate(
    median_cases   = round(med, 0),
    lo             = round(lo, 0),
    hi             = round(hi, 0),
    median_per100k = round(per100k(med), 3)
  ) %>%
  select(source, median_cases, lo, hi, median_per100k) %>%
  arrange(desc(median_cases))

cat("\n=== Posterior attribution (cumulative, weeks 0..", max(obs_weeks), ", Canada) ===\n", sep = "")
print(attr_df)

## ------------------- Sensitivity (Spearman) ---------------
## Outcome = total modelled human infections, H + W, not reported cases

post_to_cum_infections <- function(ii) {
  par_i <- params
  for (nm in par_names) {
    par_i[nm] <- posterior_draws[[nm]][ii]
  }

  wk <- simulate_weekly_human_counts(
    par = par_i,
    times = times,
    init_state = init_state,
    obs_weeks = obs_weeks,
    N_H = N_H,
    N_W = N_W,
    rho = 1,
    outcome = "infections",
    DO_DELAY = FALSE
  )

  sum(wk$total)
}

Ksens <- min(300, nrow(posterior_draws))
set.seed(404)
idx_sens <- sample(seq_len(nrow(posterior_draws)), Ksens)

cum_cases_vec <- purrr::map_dbl(idx_sens, post_to_cum_infections)

sens_df <- as_tibble(posterior_draws[idx_sens, par_names])
Ksens
spearman <- map_dbl(names(sens_df), ~ {
  x <- sens_df[[.x]]
  if (sd(x) == 0 || sd(cum_cases_vec) == 0) return(NA_real_)
  cor(x, cum_cases_vec, method = "spearman")
})

spearman <- set_names(spearman, names(sens_df))

sens_tbl <- tibble::enframe(spearman, name = "parameter", value = "rho") %>%
  arrange(desc(abs(rho)))

cat("\n=== Sensitivity (Spearman ρ with cumulative modelled human infections, H+W) ===\n")
print(sens_tbl)

## ------------------- Interventions (ARR/RRR) --------------
scenarios <- list(
  list(name = "PPE_50", label = "Enhanced PPE & Hygiene Among Poultry Workers"),
  list(name = "Backyard_25",  label = "Reduced Community Level Exposure"),
  list(name = "CP_biosec_30", label = "Enhanced Biosecurity"),
  list(name = "Env_decon_2x", label = "Environmental Management"),
  list(name = "Depop_Emergency (3days)", label = "Depopulation of Infected Commercial Farms")
)

apply_intervention <- function(par, scen) {
  nm <- as.character(scen$name); p <- par
  switch(
    nm,
    "PPE_50" = {
      keys <- c(
        "beta_WBW", "beta_BPW", "beta_CPLW", "beta_CPSW", "beta_EW", "beta_XW"
      )

      keys <- intersect(keys, names(p))
      p[keys] <- p[keys] * 0.5
      p
    },
    "Backyard_25" = {
      keys <- c("beta_BPBP","beta_BPH","beta_BPW","beta_BPX","beta_WBBP")
      keys <- intersect(keys, names(p))
      p[keys] <- p[keys] * 0.75
      p
    },
    "CP_biosec_30" = {
      for (k in c("beta_CPL_CPL", "beta_CPLW", "beta_WBCPL"))
        if (k %in% names(p)) p[k] <- p[k] * 0.70
      for (k in c("beta_CPS_CPS", "beta_CPSW", "beta_WBCPS"))
        if (k %in% names(p)) p[k] <- p[k] * 0.85
      p
    },
    "Env_decon_2x" = {
      if ("delta_env" %in% names(p)) p["delta_env"] <- p["delta_env"] * 2
      p
    },
    "Depop_Emergency (3days)" = {
      p["k_cull_CPL"] <- 0.33
      p["k_cull_CPS"] <- 0.33
      p
    },
    p
  )
}

cum_infections_from_par <- function(par) {
  wk <- simulate_weekly_human_counts(
    par = par,
    times = times,
    init_state = init_state,
    obs_weeks = obs_weeks,
    N_H = N_H,
    N_W = N_W,
    rho = 1,
    outcome = "infections",
    DO_DELAY = FALSE
  )

  sum(wk$total)
}

# PHSMs are evaluated on the same 400 posterior draws as Table 4 (idx, PART 2),
# so Tables 4 and 5 share one baseline
idx_int <- idx
Kint    <- length(idx_int)
baseline_cum <- numeric(Kint)
int_cum      <- array(NA_real_, dim = c(Kint, length(scenarios)))
colnames(int_cum) <- vapply(scenarios, function(s) s$name, character(1))

for (i in seq_along(idx_int)) {
  ii   <- idx_int[i]
  par0 <- params
  for (nm in par_names) par0[nm] <- posterior_draws[[nm]][ii]
  rho0 <- posterior_draws$rho[ii]
  baseline_cum[i] <- cum_infections_from_par(par0)

  for (j in seq_along(scenarios)) {
    par_j <- apply_intervention(par0, scenarios[[j]])
    int_cum[i, j] <- cum_infections_from_par(par_j)
  }
}

reduct_mat <- 1 - sweep(int_cum, 1, baseline_cum, "/")
arr_mat    <- matrix(baseline_cum, nrow = length(baseline_cum), ncol = ncol(int_cum)) - int_cum
colnames(arr_mat) <- colnames(int_cum)
valid_rows <- is.finite(baseline_cum) & baseline_cum > 0
reduct_mat <- reduct_mat[valid_rows, , drop = FALSE]
arr_mat    <- arr_mat[valid_rows, , drop = FALSE]

summarise_draws <- function(x) c(
  med = median(x, na.rm = TRUE),
  lo  = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
  hi  = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
)

rrr_sum <- apply(reduct_mat, 2, summarise_draws)
arr_sum <- apply(arr_mat,    2, summarise_draws)

lab_map <- tibble(
  scenario = vapply(scenarios, function(s) s$name, character(1)),
  label    = vapply(scenarios, function(s) s$label, character(1))
)

rrr_df <- tibble(
  scenario = colnames(reduct_mat),
  med = as.numeric(rrr_sum["med", ]),
  lo  = as.numeric(rrr_sum["lo",  ]),
  hi  = as.numeric(rrr_sum["hi",  ])
) %>% left_join(lab_map, by = "scenario") %>%
  select(label, med, lo, hi)

arr_df <- tibble(
  scenario  = colnames(arr_mat),
  med_cases = round(as.numeric(arr_sum["med", ])),
  lo_cases  = round(as.numeric(arr_sum["lo",  ])),
  hi_cases  = round(as.numeric(arr_sum["hi",  ]))
) %>% left_join(lab_map, by = "scenario") %>%
  mutate(
    med_per100k = round(per100k(med_cases), 3),
    lo_per100k  = round(per100k(lo_cases),  3),
    hi_per100k  = round(per100k(hi_cases),  3)
  ) %>%
  select(label, med_cases, lo_cases, hi_cases, med_per100k, lo_per100k, hi_per100k)

base_med <- median(baseline_cum[valid_rows])
base_lo  <- quantile(baseline_cum[valid_rows], 0.025, names = FALSE)
base_hi  <- quantile(baseline_cum[valid_rows], 0.975, names = FALSE)

cat("\n=== Decision outcomes (Canada) ===\n")
cat("Baseline cumulative modelled human infections, H+W (weeks 0..", max(obs_weeks), "): ",
    format(round(base_med), big.mark=","), " [",
    format(round(base_lo),  big.mark=","), ", ",
    format(round(base_hi),  big.mark=","), "]  ~ ",
    round(per100k(base_med), 3), " per 100k\n", sep = "")

cat("\nRelative risk reduction (percent) — median [95% CrI]:\n")
print(rrr_df %>%
        mutate(across(c(med, lo, hi), ~ round(100*.x, 2))) %>%
        arrange(desc(med)))

cat("\nAbsolute cases averted — median [95% CrI] (counts; per 100k):\n")
print(
  arr_df %>%
    arrange(desc(med_cases)) %>%
    mutate(display = paste0(
      format(med_cases, big.mark=",")," [",
      format(lo_cases,  big.mark=","),", ",
      format(hi_cases,  big.mark=","),"] ; ",
      med_per100k," [",lo_per100k,", ",hi_per100k,"]"
    )) %>%
    select(label, display)
)

## ------------------- Deaths (TOTAL only + 95% CrI) --------
simulate_total_deaths_counts <- function(par, times, init_state, obs_weeks,
                                         N_H, N_W, IFR_H, IFR_W,
                                         DO_DELAY = TRUE){

  N_G <- N_H - N_W

  sim <- try(
    ode(y = init_state, times = times, func = multi_species_seir_env, parms = par),
    silent = TRUE
  )
  if (inherits(sim, "try-error")) return(c(Human = NA_real_, Worker = NA_real_, Total = NA_real_))

  df <- as_tibble(as.data.frame(sim))

  if (!all(c("time","E_H","E_W") %in% names(df))) {
    return(c(Human = NA_real_, Worker = NA_real_, Total = NA_real_))
  }

  incH_prop <- par["sigma_H"] * df$E_H
  incW_prop <- par["sigma_W"] * df$E_W

  dpmf <- death_pmf()

  deathsH_prop_daily <- if (DO_DELAY) conv_incidence(incH_prop * IFR_H, dpmf) else incH_prop * IFR_H
  deathsW_prop_daily <- if (DO_DELAY) conv_incidence(incW_prop * IFR_W, dpmf) else incW_prop * IFR_W

  deathsH_daily <- deathsH_prop_daily * N_G
  deathsW_daily <- deathsW_prop_daily * N_W

  wkH <- aggregate_weekly(deathsH_daily, df$time) %>%
    filter(week %in% obs_weeks) %>% pull(val)

  wkW <- aggregate_weekly(deathsW_daily, df$time) %>%
    filter(week %in% obs_weeks) %>% pull(val)

  totalH <- sum(wkH, na.rm = TRUE)
  totalW <- sum(wkW, na.rm = TRUE)

  c(Human = totalH, Worker = totalW, Total = totalH + totalW)
}

summarise_crit <- function(x){
  x <- x[is.finite(x)]
  if (length(x) < 10) return(c(med = NA_real_, lo = NA_real_, hi = NA_real_))
  c(
    med = median(x),
    lo  = as.numeric(quantile(x, 0.025, names = FALSE)),
    hi  = as.numeric(quantile(x, 0.975, names = FALSE))
  )
}

K_death <- min(300, nrow(posterior_draws))
set.seed(4219)
idx_death <- sample(seq_len(nrow(posterior_draws)), K_death)

death_tot_mat <- matrix(NA_real_, nrow = K_death, ncol = 3)
colnames(death_tot_mat) <- c("Human","Worker","Total")

for (i in seq_along(idx_death)) {
  ii <- idx_death[i]

  par_i <- params
  for (nm in par_names) par_i[nm] <- posterior_draws[[nm]][ii]

  death_tot_mat[i, ] <- simulate_total_deaths_counts(
    par = par_i,
    times = times,
    init_state = init_state,
    obs_weeks = obs_weeks,
    N_H = N_H,
    N_W = N_W,
    IFR_H = IFR_H,
    IFR_W = IFR_W,
    DO_DELAY = DO_DELAY
  )
}

death_sum <- rbind(
  Human  = summarise_crit(death_tot_mat[, "Human"]),
  Worker = summarise_crit(death_tot_mat[, "Worker"]),
  Total  = summarise_crit(death_tot_mat[, "Total"])
)

death_summary <- tibble(
  Group  = rownames(death_sum),
  median = round(death_sum[, "med"]),
  lo     = round(death_sum[, "lo"]),
  hi     = round(death_sum[, "hi"])
)

cat("\n=== TOTAL deaths (weeks 0..", max(obs_weeks), ") — median [95% CrI] ===\n", sep="")
print(death_summary)

## ------------------- PLOTS -------------------------------
# Colour-blind friendly palette: Okabe-Ito style
cb_palette <- c(
  "#0072B2", # blue
  "#E69F00", # orange
  "#009E73", # bluish green
  "#CC79A7", # reddish purple
  "#56B4E9", # sky blue
  "#D55E00", # vermillion
  "#F0E442", # yellow
  "#000000"  # black
)

line_types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "solid", "dashed")

make_cols <- function(x) {
  setNames(rep_len(cb_palette, length(x)), x)
}

make_ltys <- function(x) {
  setNames(rep_len(line_types, length(x)), x)
}

state_cols <- c(
  "S" = "#0072B2",
  "E" = "#E69F00",
  "I" = "#D55E00",
  "R" = "#009E73",
  "D" = "#000000"
)

state_ltys <- c(
  "S" = "solid",
  "E" = "dashed",
  "I" = "dotted",
  "R" = "dotdash",
  "D" = "longdash"
)

# Add per-100k columns directly to pred_df
pred_df <- pred_df %>%
  mutate(
    lo95_per100k = per100k(lo95),
    hi95_per100k = per100k(hi95),
    lo50_per100k = per100k(lo50),
    hi50_per100k = per100k(hi50),
    med_per100k  = per100k(med)
  )

# Figure 3. SEIR(D) trajectories by population
p1 <- seir_long %>%
  filter(Population != "Environment") %>%
  mutate(State = factor(State, levels = c("S", "E", "I", "R", "D"))) %>%
  ggplot(aes(x = time, y = value, colour = State, linetype = State)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = state_cols, drop = FALSE) +
  scale_linetype_manual(values = state_ltys, drop = FALSE) +
  facet_wrap(~ Population, scales = "free_y", ncol = 3) +
  labs(
    title = "SEIR(D) trajectories by population (daily)",
    x = "Day",
    y = "Proportion",
    colour = "State",
    linetype = "State"
  ) +
  theme_minimal()

print(p1)
ggsave("P1_SEIRD_by_population.png", p1, width = 12, height = 8, dpi = 300)

# Supplemental Figure 1. Environmental reservoir
p1b <- seir_long %>%
  filter(State == "Env") %>%
  ggplot(aes(x = time, y = value)) +
  geom_line(linewidth = 0.9, colour = "#0072B2") +
  labs(
    title = "Environmental reservoir (daily)",
    x = "Day",
    y = "Level (arbitrary units)"
  ) +
  theme_minimal()

print(p1b)
ggsave("P1b_Environment.png", p1b, width = 8, height = 4, dpi = 300)

# Figure 2. Weekly reported cases vs synthetic calibration targets
# Posterior predictive intervals include negative-binomial sampling variation
obs_df <- tibble(
  week = obs_weeks,
  y = as.numeric(obs_counts_weekly),
  y_per100k = per100k(as.numeric(obs_counts_weekly))
)

phi_i <- posterior_draws$phi[idx]
set.seed(20261020)
y_sim <- matrix(
  rnbinom(length(pred_mat_reported),
          size = rep(phi_i, times = ncol(pred_mat_reported)),
          mu   = as.vector(pred_mat_reported)),
  nrow = nrow(pred_mat_reported)
)
ppi_df <- tibble(
  week = obs_weeks,
  med  = apply(pred_mat_reported, 2, median),
  lo50 = apply(y_sim, 2, quantile, 0.25),
  hi50 = apply(y_sim, 2, quantile, 0.75),
  lo95 = apply(y_sim, 2, quantile, 0.025),
  hi95 = apply(y_sim, 2, quantile, 0.975)
)
cov95 <- sum(obs_df$y >= ppi_df$lo95 & obs_df$y <= ppi_df$hi95)
cov50 <- sum(obs_df$y >= ppi_df$lo50 & obs_df$y <= ppi_df$hi50)
cat("\nPosterior predictive coverage: ", cov95, " of ", nrow(obs_df),
    " weekly targets inside the 95% interval; ", cov50, " inside the 50% interval\n", sep = "")

p2 <- ggplot(data = ppi_df, mapping = aes(x = week)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = "95% posterior predictive interval")) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50, fill = "50% posterior predictive interval")) +
  geom_line(aes(y = med, colour = "Posterior predictive median"), linewidth = 1.0) +
  geom_point(data = obs_df, aes(x = week, y = y, shape = "Synthetic calibration target"), size = 2, colour = "#000000") +
  scale_fill_manual(values = c(
    "95% posterior predictive interval" = "grey80",
    "50% posterior predictive interval" = "grey60"
  )) +
  scale_colour_manual(values = c(
    "Posterior predictive median" = "#0072B2"
  )) +
  scale_shape_manual(values = c(
    "Synthetic calibration target" = 16
  )) +
  labs(
    x = "Week",
    y = "Reported cases",
    fill = "",
    colour = "",
    shape = ""
  ) +
  theme_minimal()

print(p2)
ggsave("Figure2_posterior_predictive.png", p2, width = 10, height = 6, dpi = 300, bg = "white")

# Weekly reported cases per 100 000
p2b <- ggplot(data = pred_df, mapping = aes(x = week)) +
  geom_ribbon(aes(ymin = lo95_per100k, ymax = hi95_per100k, fill = "95% credible interval"), alpha = 0.20) +
  geom_ribbon(aes(ymin = lo50_per100k, ymax = hi50_per100k, fill = "50% credible interval"), alpha = 0.35) +
  geom_line(aes(y = med_per100k, colour = "Posterior predictive median"), linewidth = 1.0) +
  geom_point(data = obs_df, aes(x = week, y = y_per100k, shape = "Synthetic calibration target"), size = 2, colour = "#000000") +
  scale_fill_manual(values = c(
    "95% credible interval" = "grey75",
    "50% credible interval" = "grey55"
  )) +
  scale_colour_manual(values = c(
    "Posterior predictive median" = "#0072B2"
  )) +
  scale_shape_manual(values = c(
    "Synthetic calibration target" = 16
  )) +
  labs(
    title = "Weekly reported human cases per 100k — posterior predictive",
    subtitle = "Compared with synthetic calibration targets",
    x = "Week",
    y = "Reported cases per 100k",
    fill = "",
    colour = "",
    shape = ""
  ) +
  theme_minimal()

print(p2b)
ggsave("P2b_weekly_reported_cases_per100k.png", p2b, width = 10, height = 6, dpi = 300)

# Total deaths with 95% CrI
p2c <- ggplot(death_summary, aes(x = Group, y = median, fill = Group)) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  scale_fill_manual(values = make_cols(unique(death_summary$Group))) +
  coord_flip() +
  labs(
    title = paste0("Total deaths, weeks 0–", max(obs_weeks)),
    subtitle = "Median with 95% credible interval",
    x = "",
    y = "Deaths"
  ) +
  theme_minimal()

print(p2c)
ggsave("P2c_total_deaths_CrI.png", p2c, width = 8, height = 5, dpi = 300)

# Figure 4. Weekly source attribution
p3 <- attrib_week_med %>%
  mutate(source = factor(source, levels = sources)) %>%
  ggplot(aes(x = factor(week), y = counts, fill = source)) +
  geom_col(position = position_dodge2(width = 0.85, preserve = "single"), width = 0.75) +
  scale_fill_manual(values = make_cols(sources), drop = FALSE) +
  labs(
    title = "Attribution of weekly modelled human infections",
    x = "Week",
    y = "Modelled infections",
    fill = "Source"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(p3)
ggsave("P3_attribution_grouped_bars.png", p3, width = 12, height = 6, dpi = 300)

# Cumulative attribution (Supplemental Table 2)
p3b <- attr_df %>%
  mutate(source = factor(source, levels = attr_df$source[order(attr_df$median_cases)])) %>%
  ggplot(aes(x = source, y = median_cases, fill = source)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = make_cols(levels(factor(attr_df$source)))) +
  coord_flip() +
  labs(
    title = "Cumulative modelled human infections by source",
    x = "",
    y = paste0("Modelled infections, cumulative over weeks 0–", max(obs_weeks))
  ) +
  theme_minimal()

print(p3b)
ggsave("P3b_cumulative_attribution.png", p3b, width = 8, height = 6, dpi = 300)

# Sensitivity analysis (Supplemental Table 3)
if (all(is.na(sens_tbl$rho))) {
  message("Sensitivity coefficients all NA; skipping P4 plot.")
} else {
  p4 <- sens_tbl %>%
    mutate(
      abs_rho = abs(rho),
      parameter = forcats::fct_reorder(parameter, abs_rho, .na_rm = TRUE)
    ) %>%
    ggplot(aes(x = parameter, y = abs_rho, fill = parameter)) +
    geom_col(na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(values = make_cols(sens_tbl$parameter)) +
    coord_flip() +
    labs(
      title = "Sensitivity: |Spearman rho| with cumulative modelled human infections",
      x = "Parameter",
      y = "|Spearman rho|"
    ) +
    theme_minimal()

  print(p4)
  ggsave("P4_sensitivity_spearman.png", p4, width = 8, height = 6, dpi = 300)
}

############################################################
# PART 4  One-year extension: restocking and seasonal wild-bird importation
############################################################

## ------------------- One-year extension (restocking, wild-bird importation) -------
## Main analysis: Eq. 18 burden-linked spillover is switched off (kappa_link = 0).
## The Eq. 18 sensitivity analysis (kappa_link = 1e-6) is in script 2.

DO_STOCH  <- TRUE
N_REP     <- 200      # replicates per scenario
SEED_BASE <- 20261020 # replicate r uses seed SEED_BASE + r in every scenario (common random numbers)
YEAR_DAYS <- 365

# Configuration for spillover hazard and seasonal WB import
stoch_cfg <- list(
  eta0       = 0.0,
  kappa_link = 0,        # Eq. 18 off in the main analysis; 1e-06 in the sensitivity analysis

  w_WB  = 1.0,
  w_BP  = 1.0,
  w_CPL = 1.5,
  w_CPS = 0.9,
  w_X   = 0.6,
  w_E   = 0.3,

  imp_WB_base   = 5e-05,
  imp_WB_amp    = 0.8,
  imp_WB_phase  = 0,
  imp_pop_scale = 1e6,
  # Restocking of depopulated commercial premises (mean 60 days after depopulation)
  omega_restock = 1/60,
  # optional wild-bird population turnover (0 = off)
  omega_WB = 0
)

# Expected number of additional primary zoonotic infections per day
spillover_rate_cases <- function(par, state_row, t, cfg = stoch_cfg) {
  with(as.list(c(state_row, par, cfg)), {
    Env_eff <- Env / (1 + k_env_sat * Env)

    linked <- kappa_link * (
      w_WB  * I_WB +
        w_BP  * I_BP +
        w_CPL * I_CPL +
        w_CPS * I_CPS +
        w_X   * I_X +
        w_E   * Env_eff
    )

    rate <- N_H * pmax(eta0 + linked, 0)
    rate
  })
}

# One-year simulation (daily ODE steps with Poisson importation)
# Outcome = modelled human infections, H + W, not reported cases
simulate_one_year_stoch <- function(par, init,
                                    n_days = YEAR_DAYS,
                                    cfg = stoch_cfg) {

  N_G <- N_H - N_W
  # Restocking and wild-bird turnover apply to the one-year horizon only
  par["omega_CPL"] <- cfg$omega_restock
  par["omega_CPS"] <- cfg$omega_restock
  par["omega_WB"]  <- cfg$omega_WB

  states <- matrix(NA_real_, nrow = n_days + 1, ncol = length(init))
  colnames(states) <- names(init)
  states[1, ] <- init

  daily_model_infections <- numeric(n_days)
  daily_spill_infections <- numeric(n_days)

  for (d in 1:n_days) {

    # 1. Step ODE forward by one day
    sol <- ode(
      y = states[d, ],
      times = c(0, 1),
      func = multi_species_seir_env,
      parms = par
    )

    out <- as.numeric(sol[2, -1])
    names(out) <- colnames(states)

    # 2. Within-model human infections, H + W
    inc_prop_H <- par["sigma_H"] * out["E_H"]
    inc_prop_W <- par["sigma_W"] * out["E_W"]

    daily_model_infections[d] <- inc_prop_H * N_G + inc_prop_W * N_W

    # 3. Seasonal wild-bird importation (Poisson)
    # Imported birds move from S_WB to E_WB
    season_factor <- 1 + cfg$imp_WB_amp *
      cos(2 * pi * ((d - 1 - cfg$imp_WB_phase) / 365))
    season_factor <- max(season_factor, 0)
    lambda_imp <- cfg$imp_WB_base * season_factor * cfg$imp_pop_scale
    imp_prop   <- rpois(1, lambda = lambda_imp) / cfg$imp_pop_scale
    imp_prop   <- min(imp_prop, out[["S_WB"]])
    out["S_WB"] <- out["S_WB"] - imp_prop
    out["E_WB"] <- out["E_WB"] + imp_prop

    # 4. Additional primary zoonotic spillovers into humans
    rate <- spillover_rate_cases(par, out, t = d - 1, cfg = cfg)
    extra_infections <- rpois(1, lambda = rate)

    daily_spill_infections[d] <- extra_infections

    # 5. Store updated state
    states[d + 1, ] <- out
  }

  list(
    states = states,
    model_infections = daily_model_infections,
    spill_infections = daily_spill_infections
  )
}

# Helper for weekly aggregation in the one-year simulations
weekly_vec_fixed <- function(daily, wk) {
  out <- aggregate_weekly(daily, 0:(length(daily) - 1)) %>%
    right_join(tibble(week = wk), by = "week") %>%
    arrange(week)

  out$val[is.na(out$val)] <- 0
  out$val
}

if (DO_STOCH) {

  par_med_all <- par_med

  stoch_scenarios <- c(
    "Baseline",
    "PPE_50",
    "Backyard_25",
    "CP_biosec_30",
    "Env_decon_2x",
    "Depop_Emergency (3days)"
  )

  get_scen_label <- function(s) {
    if (s == "Baseline") return("Baseline (no intervention)")

    if (exists("lab_map")) {
      lab <- lab_map$label[match(s, lab_map$scenario)]
      if (!is.na(lab)) return(lab)
    }

    s
  }

  wk <- 0:floor((YEAR_DAYS - 1) / 7)

  stoch_week_list   <- list()
  annual_cases_list <- list()
  annual_reps       <- list()

  for (sc in stoch_scenarios) {

    cat("\nRunning one-year scenario:", sc, "\n")

    par_s <- par_med_all

    if (sc != "Baseline") {
      par_s <- apply_intervention(par_s, list(name = sc))
    }

    wk_mat     <- matrix(0, nrow = N_REP, ncol = length(wk))
    annual_vec <- numeric(N_REP)

    for (r in 1:N_REP) {

      if (r %% 50 == 0) {
        cat("  replicate", r, "of", N_REP, "\n")
      }
      set.seed(SEED_BASE + r)

      st <- simulate_one_year_stoch(
        par = par_s,
        init = init_state,
        n_days = YEAR_DAYS,
        cfg = stoch_cfg
      )

      model_week <- weekly_vec_fixed(st$model_infections, wk)
      spill_week <- weekly_vec_fixed(st$spill_infections, wk)

      total_week <- model_week + spill_week

      wk_mat[r, ] <- total_week
      annual_vec[r] <- sum(total_week)
    }
    annual_reps[[sc]] <- annual_vec

    stoch_week_list[[sc]] <- tibble(
      week = wk,
      med  = apply(wk_mat, 2, median),
      lo95 = apply(wk_mat, 2, quantile, 0.025),
      hi95 = apply(wk_mat, 2, quantile, 0.975),
      scenario = get_scen_label(sc)
    )

    annual_cases_list[[sc]] <- tibble(
      scenario = get_scen_label(sc),
      med  = median(annual_vec),
      lo95 = quantile(annual_vec, 0.025),
      hi95 = quantile(annual_vec, 0.975),
      med_per100k  = per100k(median(annual_vec)),
      lo95_per100k = per100k(quantile(annual_vec, 0.025)),
      hi95_per100k = per100k(quantile(annual_vec, 0.975))
    )
  }

  stoch_week_df <- bind_rows(stoch_week_list) %>%
    mutate(scenario = factor(scenario))

  stoch_annual_df <- bind_rows(annual_cases_list) %>%
    mutate(
      med  = round(med),
      lo95 = round(lo95),
      hi95 = round(hi95),
      med_per100k  = round(med_per100k, 3),
      lo95_per100k = round(lo95_per100k, 3),
      hi95_per100k = round(hi95_per100k, 3)
    )

  # Table 6. Annual modelled human infections. Replicates differ only through the
  # Poisson importation draws, which change annual totals by less than one infection,
  # so the table reports the median without intervals.
  base_reps <- annual_reps[["Baseline"]]
  table6 <- map_dfr(stoch_scenarios, function(sc) {
    a <- annual_reps[[sc]]
    w <- stoch_week_list[[sc]]
    tibble(
      scenario = get_scen_label(sc),
      annual   = median(a),
      per100k  = per100k(median(a)),
      averted  = if (sc == "Baseline") NA_real_ else median(base_reps - a),
      relative_reduction_pct = if (sc == "Baseline") NA_real_ else median(100 * (1 - a / base_reps)),
      plateau_per_week = mean(w$med[w$week >= 20 & w$week <= 51]),
      first_wave_peak  = max(w$med),
      peak_week        = w$week[which.max(w$med)],
      interval_width_95 = diff(quantile(a, c(0.025, 0.975), names = FALSE))
    )
  }) %>% arrange(is.na(averted) == FALSE, desc(relative_reduction_pct))

  cat("\n=== Table 6. One-year extension (Eq. 18 off), modelled human infections, H+W ===\n")
  print(table6 %>% mutate(across(c(annual, averted), round),
                          across(c(per100k), ~ round(.x, 3)),
                          across(c(relative_reduction_pct, plateau_per_week, first_wave_peak), ~ round(.x, 1)),
                          interval_width_95 = round(interval_width_95, 2)))
  cat("Largest 95% interval width across scenarios: ",
      round(max(table6$interval_width_95), 2), " infections per year\n", sep = "")
  write.csv(table6, "Table6_one_year_extension.csv", row.names = FALSE)

  # Figure 6. Weekly modelled infections, weeks 0-51 (intervals narrower than one infection)
  fig6_levels <- c("Baseline (no intervention)",
                   "Depopulation of Infected Commercial Farms",
                   "Enhanced PPE & Hygiene Among Poultry Workers",
                   "Enhanced Biosecurity",
                   "Reduced Community Level Exposure",
                   "Environmental Management")
  fig6_cols <- c("#000000", "#D55E00", "#E69F00", "#CC79A7", "#009E73", "#0072B2")
  fig6_lty  <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  names(fig6_cols) <- names(fig6_lty) <- fig6_levels

  p_stoch <- stoch_week_df %>%
    filter(week <= 51) %>%
    mutate(scenario = factor(as.character(scenario), levels = fig6_levels)) %>%
    ggplot(aes(x = week, y = med, colour = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = fig6_cols) +
    scale_linetype_manual(values = fig6_lty) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Week", y = "Modelled human infections per week",
         colour = NULL, linetype = NULL) +
    guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  print(p_stoch)
  ggsave("Figure6_one_year_extension.png", p_stoch, width = 9, height = 6, dpi = 300, bg = "white")
}

############################################################
# PART 5  Supplementary figures
############################################################

## ------------------- Technical Appendix Figure 3: posterior distributions ----------
# Calibrated parameters (Table 2; Supplemental Table 1)
param_order <- c(
  "beta_BPH",
  "beta_BPW",
  "beta_CPLW",
  "beta_CPSW",
  "beta_EH",
  "beta_EW",
  "beta_HH",
  "beta_WBH",
  "beta_WBW",
  "beta_WH",
  "beta_WW",
  "beta_XH",
  "beta_XW",
  "rho",
  "phi"
)

param_cols <- intersect(param_order, names(posterior_draws))

if (length(param_cols) == 0) {
  stop("No calibrated posterior parameter columns were found in posterior_draws.")
}

# Convert to long format
post_long <- posterior_draws %>%
  select(all_of(param_cols)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    scale_type = if_else(str_detect(parameter, "^beta_"), "log10", "natural"),
    value_plot = if_else(scale_type == "log10", log10(value), value),
    parameter_label = case_when(
      parameter == "beta_BPH"  ~ "β_BPH: Backyard poultry → human",
      parameter == "beta_BPW"  ~ "β_BPW: Backyard poultry → worker",
      parameter == "beta_CPLW" ~ "β_CPLW: Large commercial poultry → worker",
      parameter == "beta_CPSW" ~ "β_CPSW: Small commercial poultry → worker",
      parameter == "beta_EH"   ~ "β_EH: Environment → human",
      parameter == "beta_EW"   ~ "β_EW: Environment → worker",
      parameter == "beta_HH"   ~ "β_HH: Human → human",
      parameter == "beta_WBH"  ~ "β_WBH: Wild bird → human",
      parameter == "beta_WBW"  ~ "β_WBW: Wild bird → worker",
      parameter == "beta_WH"   ~ "β_WH: Worker → human",
      parameter == "beta_WW"   ~ "β_WW: Worker → worker",
      parameter == "beta_XH"   ~ "β_XH: Reassortment host → human",
      parameter == "beta_XW"   ~ "β_XW: Reassortment host → worker",
      parameter == "rho"       ~ "ρ: Reporting probability",
      parameter == "phi"       ~ "φ: Negative binomial dispersion",
      TRUE ~ parameter
    )
  ) %>%
  filter(is.finite(value_plot))

# Summary statistics for vertical lines and shaded 95% CrI
post_summary <- post_long %>%
  group_by(parameter, parameter_label, scale_type) %>%
  summarise(
    q025 = quantile(value_plot, 0.025, na.rm = TRUE),
    median = median(value_plot, na.rm = TRUE),
    q975 = quantile(value_plot, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Build density data manually so the 95% CrI region can be shaded
density_df <- post_long %>%
  group_by(parameter, parameter_label, scale_type) %>%
  group_modify(~{
    d <- density(.x$value_plot, na.rm = TRUE)
    tibble(x = d$x, y = d$y)
  }) %>%
  ungroup() %>%
  left_join(post_summary, by = c("parameter", "parameter_label", "scale_type")) %>%
  mutate(in_95cri = x >= q025 & x <= q975)

# Plot
p_posterior_distributions <- ggplot(density_df, aes(x = x, y = y)) +
  geom_area(fill = "grey85", colour = NA) +
  geom_area(
    data = density_df %>% filter(in_95cri),
    fill = "grey55",
    colour = NA
  ) +
  geom_vline(
    data = post_summary,
    aes(xintercept = median),
    linewidth = 0.6
  ) +
  geom_vline(
    data = post_summary,
    aes(xintercept = q025),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    data = post_summary,
    aes(xintercept = q975),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  facet_wrap(~ parameter_label, scales = "free", ncol = 3) +
  labs(
    title = "Posterior distributions of calibrated model parameters",
    subtitle = "Dark shading indicates the 95% credible interval; solid vertical line indicates the posterior median",
    x = "Posterior value; transmission coefficients shown on log10 scale",
    y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(size = 8),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

print(p_posterior_distributions)

ggsave(
  filename = "TA_Figure3_posterior_distributions.png",
  plot = p_posterior_distributions,
  width = 13,
  height = 10,
  dpi = 300,
  bg = "white"
)

## ------------------- Figure 5: relative impact of PHSMs ----------------------------
## Uses rrr_df and arr_df (15-week PHSM analysis on the 400 draws of Table 4, PART 3)
## and table6 (one-year extension, PART 4).
library(patchwork)

lab_key <- tibble::tribble(
  ~label,                                          ~short,
  "Depopulation of Infected Commercial Farms",     "A. Depopulation of infected\ncommercial farms",
  "Enhanced PPE & Hygiene Among Poultry Workers",  "B. Enhanced PPE & hygiene\namong poultry workers",
  "Enhanced Biosecurity",                          "C. Enhanced biosecurity",
  "Reduced Community Level Exposure",              "D. Reduced community-level\nexposure",
  "Environmental Management",                      "E. Environmental management"
)
lab_levels <- rev(lab_key$short)

fig5_rrr <- rrr_df %>%
  inner_join(lab_key, by = "label") %>%
  mutate(across(c(med, lo, hi), ~ 100 * .x),
         short = factor(short, levels = lab_levels))

fig5_arr <- arr_df %>%
  inner_join(lab_key, by = "label") %>%
  mutate(short = factor(short, levels = lab_levels))

fig5_base <- table6 %>% filter(grepl("^Baseline", scenario))
fig5_ann <- table6 %>%
  inner_join(lab_key, by = c("scenario" = "label")) %>%
  mutate(lab_txt = sprintf("%s (\u2212%.1f%%)", format(round(annual), big.mark = ","), relative_reduction_pct),
         short   = factor(short, levels = lab_levels))

col_pt <- "#2b6cb0"; ink <- "#1f2328"; muted <- "#6e7781"
theme_fig5 <- theme_classic(base_size = 10) +
  theme(
    panel.grid.major.x  = element_line(colour = "#e5e7eb", linewidth = 0.4),
    axis.line           = element_line(colour = "#9ca3af"),
    axis.ticks          = element_line(colour = "#9ca3af"),
    axis.text           = element_text(colour = muted),
    axis.text.y         = element_text(colour = ink, size = 10),
    axis.title.x        = element_text(colour = ink, size = 9.5),
    plot.title          = element_text(face = "bold", size = 10.5, colour = ink),
    plot.title.position = "plot"
  )

pA <- ggplot(fig5_rrr, aes(y = short, x = med)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.9, colour = col_pt) +
  geom_point(size = 3, shape = 21, fill = col_pt, colour = "white", stroke = 0.8) +
  geom_text(aes(x = hi, label = sprintf("%.1f%%", med)), nudge_y = 0.28, hjust = 0, size = 3, colour = ink) +
  scale_x_continuous(limits = c(0, 32), expand = c(0, 0)) +
  labs(title = "A  15 weeks: relative reduction",
       x = "Relative reduction (%, median, 95% CrI)", y = NULL) +
  theme_fig5

pB <- ggplot(fig5_arr, aes(y = short, x = med_cases)) +
  geom_errorbarh(aes(xmin = lo_cases, xmax = hi_cases), height = 0, linewidth = 0.9, colour = col_pt) +
  geom_point(size = 3, shape = 21, fill = col_pt, colour = "white", stroke = 0.8) +
  geom_text(aes(x = hi_cases, label = med_cases), nudge_y = 0.28, hjust = 0, size = 3, colour = ink) +
  scale_x_continuous(limits = c(0, 340), expand = c(0, 0)) +
  labs(title = "B  15 weeks: infections averted",
       x = "Infections averted (median, 95% CrI)", y = NULL) +
  theme_fig5 + theme(axis.text.y = element_blank())

pC <- ggplot(fig5_ann, aes(y = short, x = annual)) +
  geom_vline(xintercept = fig5_base$annual, linetype = "dashed", colour = muted, linewidth = 0.6) +
  annotate("text", x = fig5_base$annual + 15, y = 5.45, hjust = 0, vjust = 1, size = 2.9, colour = muted,
           label = paste0("Baseline\n", format(round(fig5_base$annual), big.mark = ","))) +
  geom_point(size = 2.8, shape = 22, fill = col_pt, colour = "white", stroke = 0.8) +
  geom_text(aes(label = lab_txt), nudge_y = 0.28, size = 3, colour = ink) +
  scale_x_continuous(limits = c(1200, 1950), expand = c(0, 0), labels = scales::comma) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(title = "C  One-year extension: annual infections",
       x = "Annual modelled human infections", y = NULL) +
  theme_fig5 + theme(axis.text.y = element_blank())

fig5 <- pA + pB + pC + plot_layout(widths = c(1, 1, 1.05))
print(fig5)
ggsave("Figure5_PHSM_impact.png", fig5, width = 12.5, height = 3.9, dpi = 300, bg = "white")
