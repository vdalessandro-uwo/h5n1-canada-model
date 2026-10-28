############################################################
# Quantifying Human Risk from H5N1 in Canada:
# multi-species SEIR(D) model with an environmental reservoir
#
# Script 3 of 3: the simulation study (parameter recovery and identifiability)
#
# Run after 01_main_analysis.R, in the same R session (script 2 is not needed).
# Results go to outputs/recovery and are summarised in
# outputs/recovery/recovery_summary.txt.
#
#   R0  The "true world": every scenario output at the data-generating
#       parameter values (attribution, worker risk, Table 4 quantities,
#       15-week PHSM effects, one-year extension, burden decomposition),
#       with Figures 4-6 redrawn at these values
#   R1  Recovery table (new Table 2): prior, data-generating value, posterior,
#       width ratio, recovered yes/no  - uses the fit from script 1
#   R2  beta_HH experiment: targets generated with beta_HH = 0, 1e-8, 1e-7,
#       1e-6, 1e-2, 1e-1; each re-fitted with the main priors
#   R3  Recovery experiment: beta_WBH, beta_EH, beta_BPH, beta_CPLW each at
#       0.1x and 10x the generating value, plus three replicate datasets at
#       the generating values with new noise seeds; new Figure 7
#   R4  One fit with weak priors (SD 2 on every log10 beta) on the main
#       synthetic targets (4 chains x 120,000 iterations, burn-in 40,000)
#
# Every fit is saved as it finishes (outputs/recovery/fits/*.rds) and is
# re-used if the script is run again, so an interrupted run can resume.
# Run time: about 5 hours (18 fits of 4 chains), plus about 1 hour for R0.
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ggplot2); library(parallel); library(patchwork)
})

needed <- c("params", "par_names", "posterior_draws", "init_state", "times", "obs_weeks",
            "obs_counts_weekly", "rho_true", "phi_true", "multi_species_seir_env",
            "conv_incidence", "delay_pmf", "aggregate_weekly", "run_mcmc_chain",
            "prior_mean_log10", "prior_sd_log10", "n_iter", "burn_in", "thin", "n_chains",
            "attrib_weekly_from_par", "simulate_weekly_human_counts", "simulate_total_deaths_counts",
            "cum_infections_from_par", "apply_intervention", "scenarios", "lab_map",
            "simulate_one_year_stoch", "weekly_vec_fixed", "stoch_cfg", "YEAR_DAYS", "SEED_BASE")
missing <- needed[!vapply(needed, exists, logical(1))]
if (length(missing) > 0) stop("Run 01_main_analysis.R first. Missing: ", paste(missing, collapse = ", "))

if (!exists("RUN_R0")) RUN_R0 <- TRUE   # set to FALSE before sourcing to skip
if (!exists("RUN_R1")) RUN_R1 <- TRUE   # set to FALSE before sourcing to skip
if (!exists("RUN_R2")) RUN_R2 <- TRUE   # set to FALSE before sourcing to skip
if (!exists("RUN_R3")) RUN_R3 <- TRUE   # set to FALSE before sourcing to skip
if (!exists("RUN_R4")) RUN_R4 <- TRUE   # set to FALSE before sourcing to skip

N_REP_R  <- 200          # replicates for the one-year extension
K_ATTR   <- 200          # posterior draws used for attribution shares of each fit
N_CORES  <- max(1, min(4, parallel::detectCores() - 1, na.rm = TRUE))

out_dir <- file.path("outputs", "recovery")
fit_dir <- file.path(out_dir, "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

summary_lines <- character(0)
summary_file <- if (all(RUN_R0, RUN_R1, RUN_R2, RUN_R3, RUN_R4)) "recovery_summary.txt" else "recovery_summary_partial.txt"
WEAK_ITER <- 120000; WEAK_BURN <- 40000
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  summary_lines <<- c(summary_lines, txt)
  writeLines(summary_lines, file.path(out_dir, summary_file))
}
fmt  <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")
sci  <- function(x) ifelse(x == 0, "0", formatC(x, format = "e", digits = 2))
q95  <- function(x) quantile(x, c(0.5, 0.025, 0.975), names = FALSE, na.rm = TRUE)
N_G  <- N_H - N_W
gen  <- params                                  # the data-generating parameter set
R0_den <- as.numeric(params["gamma_H"] + params["mu_H"])

say("H5N1 model - simulation study: parameter recovery and identifiability (",
    format(Sys.time(), "%Y-%m-%d %H:%M"), ")")

## =====================================================================
## Helpers
## =====================================================================

# Synthetic weekly targets, generated exactly as in script 1
gen_targets <- function(par, seed = 123) {
  out <- as.data.frame(ode(y = init_state, times = times, func = multi_species_seir_env, parms = par))
  incH <- conv_incidence(par["sigma_H"] * out$E_H, delay_pmf())
  incW <- conv_incidence(par["sigma_W"] * out$E_W, delay_pmf())
  mH <- aggregate_weekly(incH, out$time) %>% filter(week %in% obs_weeks) %>% pull(val)
  mW <- aggregate_weekly(incW, out$time) %>% filter(week %in% obs_weeks) %>% pull(val)
  mu <- pmax(mH * N_G * rho_true, 1e-9) + pmax(mW * N_W * rho_true, 1e-9)
  set.seed(seed)
  rnbinom(length(mu), size = phi_true, mu = mu)
}
stopifnot(identical(as.numeric(gen_targets(gen, 123)), as.numeric(obs_counts_weekly)))

# Fit one synthetic dataset with the script-1 sampler; saved and re-used
fit_dataset <- function(y, tag, sd_override = NULL, iter = n_iter, burn = burn_in) {
  f <- file.path(fit_dir, paste0(tag, ".rds"))
  if (file.exists(f)) { message("  using saved fit ", tag); return(readRDS(f)) }
  keep <- list(obs = obs_counts_weekly, n_iter = n_iter, burn_in = burn_in,
               keep_idx = keep_idx, n_keep = n_keep, sd = prior_sd_log10)
  on.exit({
    obs_counts_weekly <<- keep$obs; n_iter <<- keep$n_iter; burn_in <<- keep$burn_in
    keep_idx <<- keep$keep_idx; n_keep <<- keep$n_keep; prior_sd_log10 <<- keep$sd
  })
  obs_counts_weekly <<- y
  n_iter <<- iter; burn_in <<- burn
  keep_idx <<- seq(burn_in + 1, n_iter, by = thin); n_keep <<- length(keep_idx)
  if (!is.null(sd_override)) prior_sd_log10 <<- setNames(rep(sd_override, length(keep$sd)), names(keep$sd))
  t0 <- Sys.time()
  message("  fitting ", tag, " (", n_chains, " chains x ", iter, " iterations) ...")
  run_chains <- function(k) run_mcmc_chain(chain_id = k, seed = 2468)
  ch <- NULL
  if (.Platform$OS.type == "unix" && N_CORES > 1) {
    ch <- tryCatch(parallel::mclapply(seq_len(n_chains), run_chains, mc.cores = N_CORES),
                   error = function(e) NULL)
    if (!is.null(ch) && any(vapply(ch, function(x) inherits(x, "try-error") || is.null(x$draws), logical(1)))) ch <- NULL
  }
  if (is.null(ch)) ch <- lapply(seq_len(n_chains), run_chains)
  rhat <- NA_real_
  if (requireNamespace("coda", quietly = TRUE)) {
    ml <- coda::mcmc.list(lapply(ch, function(x) coda::mcmc(x$diag_draws)))
    rhat <- max(coda::gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1])
  }
  res <- list(tag = tag, y = y,
              draws = bind_rows(lapply(ch, `[[`, "draws")) %>% select(all_of(par_names), phi, rho),
              acceptance = sapply(ch, `[[`, "acceptance"), rhat = rhat,
              prior_sd = prior_sd_log10,
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(res, f)
  res
}

# Attribution shares and worker risk at one parameter set
src_key <- c("Environment→Human/Worker" = "environment", "WildBird→Human/Worker" = "wild_bird",
             "CPL→Worker" = "cpl_worker", "Backyard→Human/Worker" = "backyard",
             "CPS→Worker" = "cps_worker", "X→Human/Worker" = "reassortment",
             "Human/Worker→Human/Worker" = "human_to_human")
world_at <- function(par) {
  aw <- attrib_weekly_from_par(par) %>% group_by(source) %>% summarise(n = sum(counts), .groups = "drop")
  wk <- simulate_weekly_human_counts(par, times, init_state, obs_weeks, N_H, N_W,
                                     rho = 1, outcome = "infections", DO_DELAY = FALSE)
  sh <- setNames(aw$n / sum(aw$n), src_key[aw$source])
  c(sh[unname(src_key)], total = sum(wk$total), general = sum(wk$H), workers = sum(wk$W),
    worker_risk = sum(wk$W) / N_W, general_risk = sum(wk$H) / N_G)
}
par_from_draw <- function(d, i) { p <- params; for (nm in par_names) p[nm] <- d[[nm]][i]; p }
par_median <- function(d) { p <- params; for (nm in par_names) p[nm] <- median(d[[nm]]); p }
shares_draws <- function(d, K = K_ATTR, seed = 909) {
  set.seed(seed); ii <- sample(seq_len(nrow(d)), min(K, nrow(d)))
  do.call(rbind, lapply(ii, function(i) world_at(par_from_draw(d, i))))
}

# Recovery summary of one fit against the parameter set that generated its data
recovery_rows <- function(fit, truth, label) {
  map_dfr(c(par_names, "rho", "phi"), function(nm) {
    v <- fit$draws[[nm]]; q <- q95(v)
    tv <- switch(nm, rho = rho_true, phi = phi_true, as.numeric(truth[nm]))
    if (nm %in% par_names) {
      prior_w <- 2 * qnorm(0.975) * fit$prior_sd[[nm]]
      wr <- diff(quantile(log10(v), c(0.025, 0.975), names = FALSE)) / prior_w
    } else wr <- NA_real_
    tibble(fit = label, parameter = nm, true_value = tv, post_median = q[1], post_lo = q[2], post_hi = q[3],
           covered = tv >= q[2] & tv <= q[3], ratio_median_to_true = ifelse(tv > 0, q[1] / tv, NA_real_),
           width_ratio = wr)
  })
}

## =====================================================================
## R0. The true world: every output at the data-generating values
## =====================================================================
wk_all <- 0:floor((YEAR_DAYS - 1) / 7)
run_one_year <- function(scen_names, cfg, par_base, n_rep = N_REP_R) {
  annual <- matrix(NA_real_, n_rep, length(scen_names), dimnames = list(NULL, scen_names))
  weekly <- array(NA_real_, c(n_rep, length(wk_all), length(scen_names)), dimnames = list(NULL, NULL, scen_names))
  for (r in seq_len(n_rep)) {
    if (r %% 50 == 0) message("  replicate ", r, " of ", n_rep)
    for (sc in scen_names) {
      par_s <- if (sc == "Baseline") par_base else apply_intervention(par_base, list(name = sc))
      set.seed(SEED_BASE + r)
      st <- simulate_one_year_stoch(par = par_s, init = init_state, n_days = YEAR_DAYS, cfg = cfg)
      w <- weekly_vec_fixed(st$model_infections, wk_all) + weekly_vec_fixed(st$spill_infections, wk_all)
      weekly[r, , sc] <- w; annual[r, sc] <- sum(w)
    }
  }
  list(annual = annual, weekly = weekly)
}
scen_all  <- c("Baseline", "Depop_Emergency (3days)", "PPE_50", "CP_biosec_30", "Backyard_25", "Env_decon_2x")
scen_label <- function(s) if (s == "Baseline") "Baseline (no intervention)" else lab_map$label[match(s, lab_map$scenario)]
pal6 <- c("Baseline (no intervention)" = "#000000", "Depopulation of Infected Commercial Farms" = "#D55E00",
          "Enhanced PPE & Hygiene Among Poultry Workers" = "#E69F00", "Enhanced Biosecurity" = "#CC79A7",
          "Reduced Community Level Exposure" = "#009E73", "Environmental Management" = "#0072B2")

if (RUN_R0) {
  message("R0: outputs at the data-generating values ...")
  tw <- world_at(gen); pm <- world_at(par_median(posterior_draws))
  say("\n[R0] The true world (data-generating values) against the posterior medians, 15-week window")
  say(sprintf("     Total modelled infections: %s true; %s posterior median", fmt(tw["total"]), fmt(pm["total"])))
  say(sprintf("     Workers: %s infections, risk %.2f%% true; %s, %.2f%% posterior median",
              fmt(tw["workers"]), 100 * tw["worker_risk"], fmt(pm["workers"]), 100 * pm["worker_risk"]))
  say(sprintf("     Worker/general risk ratio: %.0f true; %.0f posterior median",
              tw["worker_risk"] / tw["general_risk"], pm["worker_risk"] / pm["general_risk"]))
  attr_tab <- tibble(pathway = unname(src_key), true_share = as.numeric(tw[unname(src_key)]),
                     posterior_median_share = as.numeric(pm[unname(src_key)]),
                     true_infections = as.numeric(tw[unname(src_key)]) * tw["total"])
  write.csv(attr_tab, file.path(out_dir, "R0_attribution_true_vs_posterior.csv"), row.names = FALSE)
  for (k in seq_len(nrow(attr_tab)))
    say(sprintf("     %-16s true %5.1f%% (%s infections); posterior median %5.1f%%", attr_tab$pathway[k],
                100 * attr_tab$true_share[k], fmt(attr_tab$true_infections[k]), 100 * attr_tab$posterior_median_share[k]))

  # Table 4 quantities at the generating values
  sim <- as.data.frame(ode(y = init_state, times = times, func = multi_species_seir_env, parms = gen))
  inf_now <- sim$I_H * N_G + sim$I_W * N_W
  rep_wk <- simulate_weekly_human_counts(gen, times, init_state, obs_weeks, N_H, N_W,
                                         rho = rho_true, outcome = "reported", DO_DELAY = TRUE)
  deaths <- simulate_total_deaths_counts(gen, times, init_state, obs_weeks, N_H, N_W, IFR_H, IFR_W, DO_DELAY = TRUE)
  fs <- as.data.frame(ode(y = init_state, times = 0:365, func = multi_species_seir_env, parms = gen))
  final_size <- N_G * (init_state[["S_H"]] - tail(fs$S_H, 1)) + N_W * (init_state[["S_W"]] - tail(fs$S_W, 1))
  t4 <- tibble(
    outcome = c("Peak concurrently infectious (general + workers)", "Day of peak, general population", "Day of peak, workers",
                "Peak weekly expected reported cases", "Week of peak reported cases", "Cumulative expected reported cases",
                "Cumulative modelled infections", "Total deaths (weeks 0-14)", "Final size, 365 days"),
    value = c(max(inf_now), sim$time[which.max(sim$I_H)], sim$time[which.max(sim$I_W)],
              max(rep_wk$total), rep_wk$week[which.max(rep_wk$total)], sum(rep_wk$total),
              tw["total"], deaths[["Total"]], final_size))
  write.csv(t4, file.path(out_dir, "R0_Table4_true_world.csv"), row.names = FALSE)
  say("\n[R0] Table 4 quantities at the data-generating values")
  for (k in seq_len(nrow(t4))) say(sprintf("     %-50s %s", t4$outcome[k], fmt(t4$value[k], ifelse(t4$value[k] < 100, 1, 0))))
  say(sprintf("     Synthetic targets (observed): %s reported cases in total", fmt(sum(obs_counts_weekly))))

  # 15-week PHSM effects at the generating values (one parameter set, so no interval)
  b0 <- cum_infections_from_par(gen)
  t5 <- map_dfr(scenarios, function(s) {
    v <- cum_infections_from_par(apply_intervention(gen, s))
    tibble(scenario = s$label, infections = v, averted = b0 - v, rrr_pct = 100 * (1 - v / b0),
           averted_per100k = per100k(b0 - v))
  }) %>% arrange(desc(rrr_pct))
  write.csv(t5, file.path(out_dir, "R0_Table5_PHSM_true_world.csv"), row.names = FALSE)
  say("\n[R0] Table 5 at the data-generating values (baseline ", fmt(b0), " infections)")
  for (k in seq_len(nrow(t5))) say(sprintf("     %-46s RRR %5.1f%%; averted %s", t5$scenario[k], t5$rrr_pct[k], fmt(t5$averted[k])))

  # Weekly attribution (Figure 4 at the generating values)
  aw <- attrib_weekly_from_par(gen)
  f4 <- ggplot(aw, aes(x = week, y = counts, fill = source)) + geom_col(width = 0.8) +
    labs(x = "Week", y = "Modelled human infections", fill = "Source pathway") + theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, "Figure4_true_world.png"), f4, width = 10, height = 6, dpi = 300, bg = "white")

  # One-year extension at the generating values (Eq. 18 off), Table 6 and Figure 6
  message("R0: one-year extension at the data-generating values ...")
  res <- run_one_year(scen_all, stoch_cfg, gen)
  a <- res$annual
  t6 <- map_dfr(scen_all, function(sc) {
    wmed <- apply(res$weekly[, , sc, drop = TRUE], 2, median)
    tibble(scenario = scen_label(sc), annual = median(a[, sc]), per100k = per100k(median(a[, sc])),
           averted = if (sc == "Baseline") NA_real_ else median(a[, "Baseline"] - a[, sc]),
           rrr_pct = if (sc == "Baseline") NA_real_ else median(100 * (1 - a[, sc] / a[, "Baseline"])),
           plateau_wk20_51 = mean(wmed[wk_all >= 20 & wk_all <= 51]), peak = max(wmed), peak_week = wk_all[which.max(wmed)],
           interval_width_95 = diff(quantile(a[, sc], c(0.025, 0.975), names = FALSE)))
  })
  write.csv(t6, file.path(out_dir, "R0_Table6_one_year_true_world.csv"), row.names = FALSE)
  say("\n[R0] Table 6 at the data-generating values (one-year extension, Eq. 18 off)")
  for (k in seq_len(nrow(t6))) {
    d <- t6[k, ]
    say(sprintf("     %-46s annual %s%s; plateau %.1f/week; peak %.0f (week %d)", d$scenario, fmt(d$annual),
                ifelse(is.na(d$rrr_pct), "", sprintf("; averted %s; %.1f%%", fmt(d$averted), d$rrr_pct)),
                d$plateau_wk20_51, d$peak, d$peak_week))
  }
  say(sprintf("     Ratio of annual burden to the closed-system final size: %.2f", t6$annual[1] / final_size))
  wk6 <- map_dfr(scen_all, function(sc) tibble(scenario = scen_label(sc), week = wk_all,
                                              med = apply(res$weekly[, , sc, drop = TRUE], 2, median)))
  write.csv(wk6, file.path(out_dir, "R0_Figure6_weekly_true_world.csv"), row.names = FALSE)
  f6 <- wk6 %>% filter(week <= 51) %>% mutate(scenario = factor(scenario, levels = names(pal6))) %>%
    ggplot(aes(x = week, y = med, colour = scenario, linetype = scenario)) + geom_line(linewidth = 0.9) +
    scale_colour_manual(values = pal6) +
    scale_linetype_manual(values = setNames(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), names(pal6))) +
    scale_y_continuous(limits = c(0, NA)) + labs(x = "Week", y = "Modelled human infections per week", colour = NULL, linetype = NULL) +
    guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "Figure6_true_world.png"), f6, width = 9, height = 6, dpi = 300, bg = "white")

  # Decomposition of the annual burden at the generating values
  mk <- function(omega, imp, kappa) { cf <- stoch_cfg; cf$omega_restock <- omega; cf$imp_WB_base <- imp; cf$kappa_link <- kappa; cf }
  om <- stoch_cfg$omega_restock; im <- stoch_cfg$imp_WB_base
  dec <- imap_dfr(list("Closed system" = mk(0, 0, 0), "Restocking only" = mk(om, 0, 0),
                       "Importation only" = mk(0, im, 0), "Restocking + importation (main)" = mk(om, im, 0),
                       "Eq. 18 only" = mk(0, 0, 1e-6), "All three" = mk(om, im, 1e-6)), function(cfg, nm) {
    message("  decomposition: ", nm)
    r <- run_one_year("Baseline", cfg, gen)
    wmed <- apply(r$weekly[, , "Baseline", drop = TRUE], 2, median)
    tibble(component = nm, annual = median(r$annual[, "Baseline"]),
           plateau_wk20_51 = mean(wmed[wk_all >= 20 & wk_all <= 51]), peak = max(wmed))
  })
  write.csv(dec, file.path(out_dir, "R0_decomposition_true_world.csv"), row.names = FALSE)
  say("\n[R0] Decomposition of the annual burden at the data-generating values")
  for (k in seq_len(nrow(dec))) say(sprintf("     %-34s annual %s; plateau %.1f/week; peak %.0f/week",
                                            dec$component[k], fmt(dec$annual[k]), dec$plateau_wk20_51[k], dec$peak[k]))

  # Figure 5 at the generating values (one parameter set: points, no intervals)
  key <- c("Depopulation of Infected Commercial Farms" = "A. Depopulation of infected\ncommercial farms",
           "Enhanced PPE & Hygiene Among Poultry Workers" = "B. Enhanced PPE & hygiene\namong poultry workers",
           "Enhanced Biosecurity" = "C. Enhanced biosecurity",
           "Reduced Community Level Exposure" = "D. Reduced community-level\nexposure",
           "Environmental Management" = "E. Environmental management")
  lv <- rev(unname(key))
  d5 <- t5 %>% mutate(short = factor(key[scenario], levels = lv))
  d6 <- t6 %>% filter(!is.na(rrr_pct)) %>% mutate(short = factor(key[scenario], levels = lv))
  th <- theme_classic(base_size = 10) + theme(panel.grid.major.x = element_line(colour = "#e5e7eb"),
                                            plot.title = element_text(face = "bold", size = 10.5), plot.title.position = "plot")
  pA <- ggplot(d5, aes(y = short, x = rrr_pct)) + geom_point(size = 3, colour = "#2b6cb0") +
    geom_text(aes(label = sprintf("%.1f%%", rrr_pct)), nudge_y = 0.3, size = 3) +
    labs(title = "A  15 weeks: relative reduction", x = "Relative reduction (%)", y = NULL) + th
  pB <- ggplot(d5, aes(y = short, x = averted)) + geom_point(size = 3, colour = "#2b6cb0") +
    geom_text(aes(label = round(averted)), nudge_y = 0.3, size = 3) +
    labs(title = "B  15 weeks: infections averted", x = "Infections averted", y = NULL) + th + theme(axis.text.y = element_blank())
  pC <- ggplot(d6, aes(y = short, x = annual)) + geom_vline(xintercept = t6$annual[1], linetype = "dashed", colour = "#6e7781") +
    geom_point(size = 3, shape = 15, colour = "#2b6cb0") +
    geom_text(aes(label = sprintf("%s (−%.1f%%)", fmt(annual), rrr_pct)), nudge_y = 0.3, size = 3) +
    labs(title = "C  One-year extension: annual infections", x = "Annual modelled human infections", y = NULL) +
    th + theme(axis.text.y = element_blank())
  ggsave(file.path(out_dir, "Figure5_true_world.png"), pA + pB + pC + plot_layout(widths = c(1, 1, 1.05)),
         width = 12.5, height = 3.9, dpi = 300, bg = "white")
}

## =====================================================================
## R1. Recovery table for the main fit (new Table 2)
## =====================================================================
if (RUN_R1) {
  main_fit <- list(draws = posterior_draws, prior_sd = as.list(prior_sd_log10))
  t2 <- recovery_rows(main_fit, gen, "main") %>%
    mutate(prior_median = ifelse(parameter %in% par_names, 10^prior_mean_log10[parameter],
                                 ifelse(parameter == "rho", 0.5, 50)),
           prior_lo = ifelse(parameter %in% par_names, 10^(prior_mean_log10[parameter] - qnorm(0.975) * prior_sd_log10[parameter]),
                             ifelse(parameter == "rho", qbeta(0.025, 5, 5), exp(log(50) - qnorm(0.975) * 0.5))),
           prior_hi = ifelse(parameter %in% par_names, 10^(prior_mean_log10[parameter] + qnorm(0.975) * prior_sd_log10[parameter]),
                             ifelse(parameter == "rho", qbeta(0.975, 5, 5), exp(log(50) + qnorm(0.975) * 0.5))),
           prior_centred_on_truth = abs(log10(prior_median / true_value)) < 0.1)
  write.csv(t2, file.path(out_dir, "R1_Table2_recovery.csv"), row.names = FALSE)
  say("\n[R1] Recovery of the data-generating values by the main fit (new Table 2)")
  for (k in seq_len(nrow(t2))) {
    d <- t2[k, ]
    say(sprintf("     %-10s true %s; prior median %s; posterior %s (%s-%s); width ratio %s; %s",
                d$parameter, sci(d$true_value), sci(d$prior_median), sci(d$post_median), sci(d$post_lo), sci(d$post_hi),
                ifelse(is.na(d$width_ratio), "-", sprintf("%.2f", d$width_ratio)),
                ifelse(d$covered, "recovered (inside 95% CrI)", "NOT recovered")))
  }
}

## =====================================================================
## R2. beta_HH experiment
## =====================================================================
if (RUN_R2) {
  message("R2: beta_HH experiment ...")
  bhh <- c(0, 1e-8, 1e-7, 1e-6, 1e-2, 1e-1)
  r2 <- map_dfr(bhh, function(b) {
    p <- gen; p["beta_HH"] <- b
    y <- gen_targets(p, 123)
    fit <- fit_dataset(y, paste0("betaHH_", format(b, scientific = TRUE)))
    v <- fit$draws$beta_HH; q <- q95(v); r0 <- q95(v / R0_den)
    wr <- diff(quantile(log10(v), c(0.025, 0.975), names = FALSE)) / (2 * qnorm(0.975) * prior_sd_log10[["beta_HH"]])
    tibble(true_beta_HH = b, true_R0_HH = b / R0_den, total_cases = sum(y), modelled_infections = cum_infections_from_par(p),
           post_median = q[1], post_lo = q[2], post_hi = q[3], R0_median = r0[1], R0_lo = r0[2], R0_hi = r0[3],
           width_ratio = wr, covered = b >= q[2] & b <= q[3], max_rhat = fit$rhat)
  })
  write.csv(r2, file.path(out_dir, "R2_betaHH_experiment.csv"), row.names = FALSE)
  say("\n[R2] beta_HH experiment (main priors; prior median 1e-7)")
  for (k in seq_len(nrow(r2))) {
    d <- r2[k, ]
    say(sprintf("     true beta_HH %-8s (R0 %s): %s reported cases, %s modelled infections; posterior %s (%s-%s); R0 %s (%s-%s); width ratio %.2f; %s",
                sci(d$true_beta_HH), sci(d$true_R0_HH), fmt(d$total_cases), fmt(d$modelled_infections, 1),
                sci(d$post_median), sci(d$post_lo), sci(d$post_hi), sci(d$R0_median), sci(d$R0_lo), sci(d$R0_hi),
                d$width_ratio, ifelse(d$covered, "covered", "not covered")))
  }
}

## =====================================================================
## R3. Recovery experiment for four key coefficients + replicate seeds
## =====================================================================
if (RUN_R3) {
  message("R3: recovery experiment ...")
  base_fit <- list(draws = posterior_draws, prior_sd = as.list(prior_sd_log10))
  designs <- bind_rows(
    tibble(label = "baseline seed 123", coef = NA_character_, mult = 1, seed = 123),
    expand_grid(coef = c("beta_WBH", "beta_EH", "beta_BPH", "beta_CPLW"), mult = c(0.1, 10)) %>%
      mutate(label = paste0(coef, " x", mult), seed = 123),
    tibble(label = paste0("baseline seed ", 124:126), coef = NA_character_, mult = 1, seed = 124:126)
  )
  r3_rec <- list(); r3_share <- list()
  for (k in seq_len(nrow(designs))) {
    d <- designs[k, ]
    p <- gen; if (!is.na(d$coef)) p[d$coef] <- gen[d$coef] * d$mult
    fit <- if (d$label == "baseline seed 123") base_fit else
      fit_dataset(gen_targets(p, d$seed), gsub("[^A-Za-z0-9_.]", "_", d$label))
    r3_rec[[d$label]] <- recovery_rows(fit, p, d$label) %>% mutate(varied = d$coef, multiplier = d$mult, seed = d$seed)
    tw <- world_at(p); sd_ <- shares_draws(fit$draws)
    r3_share[[d$label]] <- map_dfr(c("environment", "wild_bird", "cpl_worker", "backyard"), function(pw) {
      q <- q95(sd_[, pw])
      tibble(fit = d$label, pathway = pw, true_share = tw[[pw]], est_median = q[1], est_lo = q[2], est_hi = q[3],
             covered = tw[[pw]] >= q[2] & tw[[pw]] <= q[3])
    }) %>% bind_rows(tibble(fit = d$label, pathway = "worker_risk", true_share = tw[["worker_risk"]],
                            est_median = median(sd_[, "worker_risk"]), est_lo = quantile(sd_[, "worker_risk"], 0.025),
                            est_hi = quantile(sd_[, "worker_risk"], 0.975),
                            covered = tw[["worker_risk"]] >= quantile(sd_[, "worker_risk"], 0.025) &
                              tw[["worker_risk"]] <= quantile(sd_[, "worker_risk"], 0.975)))
  }
  rec <- bind_rows(r3_rec); shr <- bind_rows(r3_share)
  write.csv(rec, file.path(out_dir, "R3_recovery_parameters.csv"), row.names = FALSE)
  write.csv(shr, file.path(out_dir, "R3_recovery_shares.csv"), row.names = FALSE)
  say("\n[R3] Recovery experiment: varied coefficient, estimate vs truth (main priors)")
  focus <- rec %>% filter(!is.na(varied) & parameter == varied | is.na(varied) & parameter %in% c("beta_WBH", "beta_EH", "beta_BPH", "beta_CPLW"))
  for (k in seq_len(nrow(focus))) {
    d <- focus[k, ]
    say(sprintf("     %-22s %-10s true %s; posterior %s (%s-%s); %s", d$fit, d$parameter, sci(d$true_value),
                sci(d$post_median), sci(d$post_lo), sci(d$post_hi), ifelse(d$covered, "covered", "NOT covered")))
  }
  cov_all <- rec %>% filter(parameter %in% par_names) %>% group_by(fit) %>% summarise(covered = sum(covered), n = n(), .groups = "drop")
  say("     Coefficients whose true value lies inside the 95% CrI, per dataset:")
  for (k in seq_len(nrow(cov_all))) say(sprintf("       %-22s %d of %d", cov_all$fit[k], cov_all$covered[k], cov_all$n[k]))
  say("\n[R3] Attributed shares and worker risk: estimate (95% CrI) vs truth")
  for (k in seq_len(nrow(shr))) {
    d <- shr[k, ]
    pct <- function(x) sprintf("%.1f%%", 100 * x)
    say(sprintf("     %-22s %-12s true %s; estimate %s (%s-%s); %s", d$fit, d$pathway, pct(d$true_share), pct(d$est_median),
                pct(d$est_lo), pct(d$est_hi), ifelse(d$covered, "covered", "NOT covered")))
  }

  # Figure 7: true value against estimate for the four varied coefficients
  f7d <- rec %>% filter(parameter %in% c("beta_WBH", "beta_EH", "beta_BPH", "beta_CPLW"),
                        is.na(varied) | varied == parameter) %>%
    mutate(parameter = factor(parameter, levels = c("beta_WBH", "beta_EH", "beta_BPH", "beta_CPLW"),
                              labels = c("Wild bird → human (β_WBH)", "Environment → human (β_EH)",
                                         "Backyard poultry → human (β_BPH)", "Large commercial poultry → worker (β_CPLW)")),
           dataset = ifelse(is.na(varied), "Generating values (4 noise seeds)", "Coefficient varied 0.1× or 10×"))
  pri <- tibble(parameter = levels(f7d$parameter),
                prior_median = 10^prior_mean_log10[c("beta_WBH", "beta_EH", "beta_BPH", "beta_CPLW")]) %>%
    mutate(parameter = factor(parameter, levels = levels(f7d$parameter)))
  f7 <- ggplot(f7d, aes(x = true_value, y = post_median, colour = dataset, group = fit)) +
    geom_abline(slope = 1, intercept = 0, colour = "#9ca3af") +
    geom_hline(data = pri, aes(yintercept = prior_median), linetype = "dashed", colour = "#6e7781") +
    geom_errorbar(aes(ymin = post_lo, ymax = post_hi), width = 0.08, position = position_dodge(width = 0.15)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.15)) +
    scale_x_log10() + scale_y_log10() + facet_wrap(~parameter, scales = "free") +
    scale_colour_manual(values = c("#2b6cb0", "#D55E00")) +
    labs(x = "Data-generating (true) value", y = "Posterior median and 95% CrI", colour = NULL,
         caption = "Grey line: estimate equals truth. Dashed line: prior median.") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "Figure7_recovery.png"), f7, width = 10, height = 7.5, dpi = 300, bg = "white")
}

## =====================================================================
## R4. Weak priors on every coefficient
## =====================================================================
if (RUN_R4) {
  message("R4: weak-prior fit ...")
  fit <- fit_dataset(obs_counts_weekly, paste0("weak_priors_sd2_", format(WEAK_ITER, scientific = FALSE)), sd_override = 2, iter = WEAK_ITER, burn = WEAK_BURN)
  rec4 <- recovery_rows(fit, gen, "weak priors (SD 2)")
  sd4 <- shares_draws(fit$draws); tw <- world_at(gen)
  sh4 <- map_dfr(c("environment", "wild_bird", "cpl_worker", "backyard", "worker_risk"), function(pw) {
    q <- q95(sd4[, pw]); tibble(pathway = pw, true = tw[[pw]], est_median = q[1], est_lo = q[2], est_hi = q[3])
  })
  write.csv(rec4, file.path(out_dir, "R4_weak_prior_parameters.csv"), row.names = FALSE)
  write.csv(sh4, file.path(out_dir, "R4_weak_prior_shares.csv"), row.names = FALSE)
  say("\n[R4] Weak priors (SD 2 on every log10 beta), main synthetic targets; ", format(WEAK_ITER, big.mark = ","), " iterations; max R-hat ",
      sprintf("%.2f", fit$rhat))
  for (k in seq_len(nrow(rec4))) {
    d <- rec4[k, ]
    say(sprintf("     %-10s true %s; posterior %s (%s-%s); width ratio %s; %s", d$parameter, sci(d$true_value),
                sci(d$post_median), sci(d$post_lo), sci(d$post_hi),
                ifelse(is.na(d$width_ratio), "-", sprintf("%.2f", d$width_ratio)), ifelse(d$covered, "covered", "NOT covered")))
  }
  for (k in seq_len(nrow(sh4))) say(sprintf("     share %-12s true %.1f%%; estimate %.1f%% (%.1f-%.1f%%)", sh4$pathway[k],
                                            100 * sh4$true[k], 100 * sh4$est_median[k], 100 * sh4$est_lo[k], 100 * sh4$est_hi[k]))
}

say("\nDone ", format(Sys.time(), "%Y-%m-%d %H:%M"))
message("Done. Summary in ", file.path(out_dir, summary_file))
