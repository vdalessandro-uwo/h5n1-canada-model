############################################################
# Quantifying Human Risk from H5N1 in Canada:
# multi-species SEIR(D) model with an environmental reservoir
#
# Script 2 of 2: sensitivity and supporting analyses
#
# Run after 01_main_analysis.R, in the same R session. Results are written
# to the folder "outputs" and summarised in outputs/sensitivity_summary.txt.
#
#   S1  Scale of the scenario: share of each animal population infected by
#       day 104, and plug-in 15-week totals (Results, scale paragraph)
#   S2  Decomposition of the one-year burden into restocking, wild-bird
#       importation and Eq. 18 spillover (Results, one-year extension)
#   S3  One-year extension with Eq. 18 burden-linked spillover switched on
#       (Supplemental Table 4)
#   S4  Wild-bird importation rate x1, x10, x100, with and without restocking
#       (Supplemental Table 6)
#   S5  Prior sensitivity of source attribution: three alternative prior sets,
#       each re-calibrated by MCMC (Table 7; Supplemental Table 5)
#
# Run time: S2-S4 about 1 hour; S5 about 2 hours (three MCMC runs; set (iii)
# uses 60,000 iterations per chain). Switch sections off below if needed.
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(parallel)
})

needed <- c("posterior_draws", "params", "par_names", "par_med", "init_state", "times",
            "obs_weeks", "stoch_cfg", "simulate_one_year_stoch", "weekly_vec_fixed",
            "apply_intervention", "attrib_weekly_from_par", "simulate_weekly_human_counts",
            "multi_species_seir_env", "run_mcmc_chain", "prior_mean_log10", "prior_sd_log10",
            "lab_map", "YEAR_DAYS", "SEED_BASE")
missing <- needed[!vapply(needed, exists, logical(1))]
if (length(missing) > 0) stop("Run 01_main_analysis.R first. Missing: ", paste(missing, collapse = ", "))

RUN_S1 <- TRUE
RUN_S2 <- TRUE
RUN_S3 <- TRUE
RUN_S4 <- TRUE
RUN_S5 <- TRUE

N_REP_S    <- 200        # replicates for the one-year analyses
KAPPA_EQ18 <- 1e-06      # Eq. 18 coupling constant used in the sensitivity analyses
N_CORES    <- max(1, min(4, parallel::detectCores() - 1, na.rm = TRUE))

out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

summary_lines <- character(0)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  summary_lines <<- c(summary_lines, txt)
}
fmt <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")
q95 <- function(x) quantile(x, c(0.5, 0.025, 0.975), names = FALSE, na.rm = TRUE)

say("H5N1 model - sensitivity and supporting analyses (", format(Sys.time(), "%Y-%m-%d %H:%M"), ")")

## =====================================================================
## Helper: one-year runs with common random numbers
## =====================================================================
wk_all <- 0:floor((YEAR_DAYS - 1) / 7)

run_one_year <- function(scen_names, cfg, n_rep = N_REP_S, par_base = par_med) {
  annual <- matrix(NA_real_, n_rep, length(scen_names), dimnames = list(NULL, scen_names))
  weekly <- array(NA_real_, c(n_rep, length(wk_all), length(scen_names)),
                  dimnames = list(NULL, NULL, scen_names))
  for (r in seq_len(n_rep)) {
    if (r %% 50 == 0) message("  replicate ", r, " of ", n_rep)
    for (sc in scen_names) {
      par_s <- if (sc == "Baseline") par_base else apply_intervention(par_base, list(name = sc))
      set.seed(SEED_BASE + r)
      st <- simulate_one_year_stoch(par = par_s, init = init_state, n_days = YEAR_DAYS, cfg = cfg)
      w <- weekly_vec_fixed(st$model_infections, wk_all) + weekly_vec_fixed(st$spill_infections, wk_all)
      weekly[r, , sc] <- w
      annual[r, sc] <- sum(w)
    }
  }
  list(annual = annual, weekly = weekly)
}

baseline_summary <- function(res) {
  a <- q95(res$annual[, "Baseline"])
  wmed <- apply(res$weekly[, , "Baseline", drop = TRUE], 2, median)
  tibble(annual_med = a[1], annual_lo = a[2], annual_hi = a[3],
         plateau_wk20_51 = mean(wmed[wk_all >= 20 & wk_all <= 51]),
         peak_weekly = max(wmed), peak_week = wk_all[which.max(wmed)])
}

scen_all <- c("Baseline", "Depop_Emergency (3days)", "PPE_50", "CP_biosec_30", "Backyard_25", "Env_decon_2x")
scen_label <- function(s) if (s == "Baseline") "Baseline (no intervention)" else lab_map$label[match(s, lab_map$scenario)]

## =====================================================================
## S1. Scale of the scenario at the posterior medians
## =====================================================================
if (RUN_S1) {
  wk15 <- simulate_weekly_human_counts(par_med, times, init_state, obs_weeks, N_H, N_W,
                                       rho = 1, outcome = "infections", DO_DELAY = FALSE)
  say("\n[S1] Posterior-median (plug-in) cumulative modelled infections, days 0-104: ",
      fmt(sum(wk15$total)), " (general population ", fmt(sum(wk15$H)), ", workers ", fmt(sum(wk15$W)), ")")
  say("     Plug-in worker cumulative risk: ", sprintf("%.1f%%", 100 * sum(wk15$W) / N_W))

  det104 <- as.data.frame(ode(y = init_state, times = 0:104, func = multi_species_seir_env, parms = par_med))
  last <- det104[nrow(det104), ]
  infected_share <- c(
    "Large-scale commercial poultry" = init_state[["S_CPL"]] - last$S_CPL,
    "Small-scale commercial poultry" = init_state[["S_CPS"]] - last$S_CPS,
    "Backyard poultry"               = init_state[["S_BP"]]  - last$S_BP,
    "Wild birds"                     = init_state[["S_WB"]]  - last$S_WB,
    "Reassortment hosts"             = init_state[["S_X"]]   - last$S_X
  )
  say("     Share of each population infected by day 104:")
  for (nm in names(infected_share)) say(sprintf("       %-32s %5.1f%%", nm, 100 * infected_share[[nm]]))
}

## =====================================================================
## S2. Decomposition of the one-year burden (baseline, posterior medians)
## =====================================================================
if (RUN_S2) {
  message("S2: decomposition of the one-year burden ...")
  mk <- function(omega, imp, kappa) {
    cf <- stoch_cfg; cf$omega_restock <- omega; cf$imp_WB_base <- imp; cf$kappa_link <- kappa; cf
  }
  om <- stoch_cfg$omega_restock; im <- stoch_cfg$imp_WB_base
  decomp_cfgs <- list(
    "Closed system (no restocking, importation or Eq. 18)" = mk(0,  0,  0),
    "Restocking only"                                      = mk(om, 0,  0),
    "Wild-bird importation only"                           = mk(0,  im, 0),
    "Restocking + importation (main analysis)"             = mk(om, im, 0),
    "Eq. 18 spillover only"                                = mk(0,  0,  KAPPA_EQ18),
    "Restocking + importation + Eq. 18"                    = mk(om, im, KAPPA_EQ18)
  )
  decomp <- imap_dfr(decomp_cfgs, function(cfg, nm) {
    message("  ", nm)
    baseline_summary(run_one_year("Baseline", cfg)) %>% mutate(component = nm, .before = 1)
  })
  write.csv(decomp, file.path(out_dir, "S2_decomposition.csv"), row.names = FALSE)
  say("\n[S2] Decomposition of the one-year burden (baseline, posterior medians)")
  for (k in seq_len(nrow(decomp))) {
    d <- decomp[k, ]
    say(sprintf("     %-54s annual %s; plateau weeks 20-51 %.1f/week; peak %.0f/week (week %d)",
                d$component, fmt(d$annual_med), d$plateau_wk20_51, d$peak_weekly, d$peak_week))
  }
}

## =====================================================================
## S3. Supplemental Table 4: one-year extension with Eq. 18 on
## =====================================================================
if (RUN_S3) {
  message("S3: one-year extension with Eq. 18 on ...")
  cfg_eq18 <- stoch_cfg; cfg_eq18$kappa_link <- KAPPA_EQ18
  res <- run_one_year(scen_all, cfg_eq18)
  a <- res$annual
  supp4 <- map_dfr(colnames(a), function(sc) {
    s <- q95(a[, sc])
    red <- if (sc == "Baseline") rep(NA, 3) else q95(a[, "Baseline"] - a[, sc])
    rel <- if (sc == "Baseline") rep(NA, 3) else q95(100 * (1 - a[, sc] / a[, "Baseline"]))
    tibble(scenario = scen_label(sc), annual_med = s[1], annual_lo = s[2], annual_hi = s[3],
           per100k = per100k(s[1]),
           averted_med = red[1], averted_lo = red[2], averted_hi = red[3],
           rel_reduction_med = rel[1], rel_reduction_lo = rel[2], rel_reduction_hi = rel[3])
  })
  write.csv(supp4, file.path(out_dir, "SuppTable4_Eq18_sensitivity.csv"), row.names = FALSE)
  say("\n[S3] Supplemental Table 4 - annual modelled infections with Eq. 18 on (kappa_link = ", KAPPA_EQ18, ")")
  for (k in seq_len(nrow(supp4))) {
    t <- supp4[k, ]
    line <- sprintf("     %-46s %s [%s-%s]", t$scenario, fmt(t$annual_med), fmt(t$annual_lo), fmt(t$annual_hi))
    if (!is.na(t$rel_reduction_med))
      line <- paste0(line, sprintf("; averted %s [%s-%s]; %.1f%% [%.1f-%.1f]",
                                   fmt(t$averted_med), fmt(t$averted_lo), fmt(t$averted_hi),
                                   t$rel_reduction_med, t$rel_reduction_lo, t$rel_reduction_hi))
    say(line)
  }
}

## =====================================================================
## S4. Supplemental Table 6: wild-bird importation rate
## =====================================================================
if (RUN_S4) {
  imp_rows <- list()
  for (restock in c(TRUE, FALSE)) {
    for (mult in c(1, 10, 100)) {
      cfg <- stoch_cfg
      cfg$kappa_link  <- 0
      cfg$imp_WB_base <- stoch_cfg$imp_WB_base * mult
      if (!restock) cfg$omega_restock <- 0
      lab <- sprintf("%s; importation rate x %d", if (restock) "With restocking" else "Without restocking", mult)
      message("S4: ", lab)
      imp_rows[[lab]] <- baseline_summary(run_one_year("Baseline", cfg)) %>%
        mutate(setting = lab, restocking = restock, multiplier = mult, nu0_per_day = cfg$imp_WB_base, .before = 1)
    }
  }
  imp_tab <- bind_rows(imp_rows) %>%
    group_by(restocking) %>%
    mutate(change_vs_x1_pct = 100 * (annual_med / annual_med[multiplier == 1] - 1)) %>%
    ungroup()
  write.csv(imp_tab, file.path(out_dir, "SuppTable6_importation_sensitivity.csv"), row.names = FALSE)
  say("\n[S4] Supplemental Table 6 - importation rate (baseline, Eq. 18 off)")
  for (k in seq_len(nrow(imp_tab))) {
    d <- imp_tab[k, ]
    say(sprintf("     %-42s annual %s (%+.1f%%); plateau weeks 20-51 %.1f/week; peak %.0f/week (week %d)",
                d$setting, fmt(d$annual_med), d$change_vs_x1_pct, d$plateau_wk20_51, d$peak_weekly, d$peak_week))
  }
}

## =====================================================================
## S5. Prior sensitivity of source attribution (Table 7; Supplemental Table 5)
## =====================================================================
attribution_draws <- function(draws, K = 300, seed = 909) {
  set.seed(seed)
  ii <- sample(seq_len(nrow(draws)), min(K, nrow(draws)))
  map_dfr(ii, function(i) {
    p <- params; for (nm in par_names) p[nm] <- draws[[nm]][i]
    aw <- attrib_weekly_from_par(p) %>% group_by(source) %>% summarise(n = sum(counts), .groups = "drop")
    tot <- sum(aw$n)
    wk <- simulate_weekly_human_counts(p, times, init_state, obs_weeks, N_H, N_W,
                                       rho = 1, outcome = "infections", DO_DELAY = FALSE)
    tibble(env  = aw$n[aw$source == "Environment→Human/Worker"] / tot,
           wb   = aw$n[aw$source == "WildBird→Human/Worker"]    / tot,
           cplw = aw$n[aw$source == "CPL→Worker"]               / tot,
           bp   = aw$n[aw$source == "Backyard→Human/Worker"]    / tot,
           worker_risk  = sum(wk$W) / N_W,
           general_risk = sum(wk$H) / N_G)
  })
}
row_table7 <- function(rows, set_name) {
  f <- function(x) 100 * q95(x)
  e <- f(rows$env); w <- f(rows$wb); cp <- f(rows$cplw); b <- f(rows$bp)
  wr <- f(rows$worker_risk); rr <- q95(rows$worker_risk / rows$general_risk)
  tibble(prior_set = set_name,
         environment = sprintf("%.0f%% (%.0f-%.0f)", e[1], e[2], e[3]),
         wild_bird   = sprintf("%.0f%% (%.0f-%.0f)", w[1], w[2], w[3]),
         cpl_worker  = sprintf("%.0f%% (%.0f-%.0f)", cp[1], cp[2], cp[3]),
         backyard    = sprintf("%.0f%% (%.0f-%.0f)", b[1], b[2], b[3]),
         worker_cum_risk = sprintf("%.2f%% (%.2f-%.2f)", wr[1], wr[2], wr[3]),
         worker_vs_general_risk_ratio = sprintf("%.0f (%.0f-%.0f)", rr[1], rr[2], rr[3]))
}
row_largest <- function(rows, set_name) {
  lv <- c("Environment", "Wild bird", "Large commercial poultry to worker", "Backyard poultry")
  top <- lv[max.col(as.matrix(rows[, c("env", "wb", "cplw", "bp")]), ties.method = "first")]
  pct <- round(100 * as.numeric(table(factor(top, levels = lv))) / length(top))
  tibble(prior_set = set_name, !!!setNames(as.list(pct), lv))
}

fit_prior_set <- function(mean_l10, sd_l10, iter, burn) {
  mcmc_orig <- list(n_iter = n_iter, burn_in = burn_in, keep_idx = keep_idx, n_keep = n_keep)
  pm_orig <- prior_mean_log10; ps_orig <- prior_sd_log10
  on.exit({
    n_iter   <<- mcmc_orig$n_iter;   burn_in <<- mcmc_orig$burn_in
    keep_idx <<- mcmc_orig$keep_idx; n_keep  <<- mcmc_orig$n_keep
    prior_mean_log10 <<- pm_orig;    prior_sd_log10 <<- ps_orig
  })
  n_iter   <<- iter; burn_in <<- burn
  keep_idx <<- seq(burn_in + 1, n_iter, by = thin); n_keep <<- length(keep_idx)
  prior_mean_log10 <<- mean_l10; prior_sd_log10 <<- sd_l10

  run_chains <- function(k) run_mcmc_chain(chain_id = k, seed = 2468)
  ch <- NULL
  if (.Platform$OS.type == "unix" && N_CORES > 1) {
    ch <- tryCatch(parallel::mclapply(seq_len(n_chains), run_chains, mc.cores = N_CORES),
                   error = function(e) NULL)
    if (!is.null(ch) && any(vapply(ch, function(x) inherits(x, "try-error") || is.null(x$draws), logical(1)))) ch <- NULL
  }
  if (is.null(ch)) ch <- lapply(seq_len(n_chains), run_chains)

  diag_txt <- ""
  if (requireNamespace("coda", quietly = TRUE)) {
    ml <- coda::mcmc.list(lapply(ch, function(x) coda::mcmc(x$diag_draws)))
    g <- coda::gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)$psrf
    ess <- coda::effectiveSize(ml)
    diag_txt <- sprintf("R-hat max %.2f (upper CI %.2f); ESS %.0f-%.0f", max(g[, 1]), max(g[, 2]), min(ess), max(ess))
  }
  list(draws = bind_rows(lapply(ch, `[[`, "draws")) %>% select(all_of(par_names), phi, rho),
       acceptance = round(100 * sapply(ch, `[[`, "acceptance"), 1),
       diag = diag_txt)
}

if (RUN_S5) {
  pm <- prior_mean_log10; ps <- prior_sd_log10
  prior_sets <- list(
    "(i) Wild-bird and environmental prior means halved" = list(
      mean = within(as.list(pm), { beta_WBH <- beta_WBH - 0.3; beta_EH <- beta_EH - 0.3 }),
      sd = as.list(ps), iter = n_iter, burn = burn_in),
    "(ii) Occupational prior means doubled" = list(
      mean = within(as.list(pm), {
        beta_CPLW <- beta_CPLW + 0.3; beta_CPSW <- beta_CPSW + 0.3; beta_WBW <- beta_WBW + 0.3
        beta_BPW <- beta_BPW + 0.3; beta_EW <- beta_EW + 0.3 }),
      sd = as.list(ps), iter = n_iter, burn = burn_in),
    "(iii) All prior SDs on log10 beta set to 1.0" = list(
      mean = as.list(pm),
      sd = as.list(setNames(rep(1.0, length(ps)), names(ps))),
      iter = 60000, burn = 20000)       # longer chains needed for convergence
  )

  rows_main <- attribution_draws(posterior_draws)
  table7 <- row_table7(rows_main, "Main priors")
  largest <- row_largest(rows_main, "Main priors")
  say("\n[S5] Prior sensitivity (each set re-calibrated; 4 chains, thinning ", thin, ")")
  for (nm in names(prior_sets)) {
    s <- prior_sets[[nm]]
    message("S5: MCMC under prior set ", nm, " (", s$iter, " iterations per chain) ...")
    fit <- fit_prior_set(unlist(s$mean)[names(pm)], unlist(s$sd)[names(ps)], s$iter, s$burn)
    saveRDS(fit$draws, file.path(out_dir, paste0("posterior_draws_prior_set_", sub("^\\((i+)\\).*", "\\1", nm), ".rds")))
    say("     ", nm, ": ", s$iter, " iterations, burn-in ", s$burn, "; acceptance ",
        paste(fit$acceptance, collapse = "/"), "%; ", fit$diag)
    rows <- attribution_draws(fit$draws)
    table7  <- bind_rows(table7,  row_table7(rows, nm))
    largest <- bind_rows(largest, row_largest(rows, nm))
  }
  write.csv(table7,  file.path(out_dir, "Table7_prior_sensitivity.csv"), row.names = FALSE)
  write.csv(largest, file.path(out_dir, "SuppTable5_largest_pathway_by_draw.csv"), row.names = FALSE)

  say("\n[S5] Table 7 - attributed share of modelled infections, median (95% CrI), 300 draws per set")
  for (k in seq_len(nrow(table7))) {
    t <- table7[k, ]
    say("     ", t$prior_set, ": environment ", t$environment, "; wild bird ", t$wild_bird,
        "; CPL-to-worker ", t$cpl_worker, "; backyard ", t$backyard,
        "; worker risk ", t$worker_cum_risk, "; risk ratio ", t$worker_vs_general_risk_ratio)
  }
  say("\n[S5] Supplemental Table 5 - largest pathway in each posterior draw (% of draws)")
  for (k in seq_len(nrow(largest))) {
    d <- largest[k, ]
    say(sprintf("     %-52s environment %d%%; wild bird %d%%; CPL-to-worker %d%%; backyard %d%%",
                d$prior_set, d$Environment, d$`Wild bird`, d$`Large commercial poultry to worker`, d$`Backyard poultry`))
  }
}

writeLines(summary_lines, file.path(out_dir, "sensitivity_summary.txt"))
message("Done. Summary in ", file.path(out_dir, "sensitivity_summary.txt"))
