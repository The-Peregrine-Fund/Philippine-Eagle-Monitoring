# =============================================================================
# FULL POWER ANALYSIS WITH TWO-LEVEL COVARIATE MODEL — p ESTIMATED + COVERAGE
# =============================================================================
# Uses BFGS optimizer (faster + provides Hessian efficiently)
# Parallel across 54 parameter combinations on 6 cores
# Delta method for SE(N_hat) and 95% CI coverage
# =============================================================================

library(tidyverse)
library(terra)
library(parallel)

cat("================================================================\n")
cat("POWER ANALYSIS: Two-Level Covariate Model (p estimated + coverage)\n")
cat("================================================================\n\n")

# Load suitability data
r <- rast("C:/Users/chrim/Dropbox/Current projects/PHEA/PHEA Surveys/PHEA H3/PHEA_Continuous.tif")
suit_sample <- spatSample(r, size = 50000, method = "random", na.rm = TRUE)
suit_values <- suit_sample[[1]]

cat("Suitability: min =", round(min(suit_values), 3),
    ", mean =", round(mean(suit_values), 3),
    ", max =", round(max(suit_values), 3), "\n\n")

# Progress directory
prog_dir <- "C:/Users/chrim/Dropbox/Current projects/PHEA/PHEA Surveys/PHEA H3/sim_progress"
dir.create(prog_dir, showWarnings = FALSE)
# Clear old progress files
file.remove(list.files(prog_dir, full.names = TRUE))

# Build parameter grid
scenarios <- data.frame(
  scenario = c("Low occupancy", "Moderate occupancy", "Higher occupancy"),
  beta0 = c(-3.5, -2.5, -1.8),
  beta1 = c(4.5, 4.0, 3.5),
  stringsAsFactors = FALSE
)
gamma0 <- 0
gamma1 <- 2.5
p_values <- c(0.05, 0.10, 0.15)
n_values <- c(150, 200, 250, 300, 350, 400)
M <- 1500
n_sims <- 1000
J <- 8
K_total <- 7
K_sampled <- 4

param_grid <- expand.grid(
  s = 1:nrow(scenarios),
  p_true = p_values,
  n_res6 = n_values,
  stringsAsFactors = FALSE
)

cat("Total combinations:", nrow(param_grid), "\n")
cat("Simulations per combo:", n_sims, "\n")
cat("Total iterations:", nrow(param_grid) * n_sims, "\n\n")

# Worker function — runs all sims for one parameter combo
run_combo <- function(combo_idx, param_grid, scenarios, suit_values,
                      M, K_total, K_sampled, J, n_sims,
                      gamma0, gamma1, prog_dir) {

  row <- param_grid[combo_idx, ]
  s <- row$s
  p_true <- row$p_true
  n_res6 <- row$n_res6
  beta0 <- scenarios$beta0[s]
  beta1 <- scenarios$beta1[s]
  t_max <- J                                     # max survey duration (visits = time units)
  lambda_true <- -log(1 - p_true)               # hazard per visit; P(found in J visits) = 1-(1-p)^J

  N_estimates      <- numeric(n_sims)
  N_estimates_null <- numeric(n_sims)
  lambda_estimates <- numeric(n_sims)
  true_N <- numeric(n_sims)
  converged <- logical(n_sims)
  ci_covers <- logical(n_sims)

  for (sim in 1:n_sims) {

    # --- Generate landscape ---
    suit_res7_all <- matrix(
      sample(suit_values, M * K_total, replace = TRUE),
      nrow = M, ncol = K_total
    )
    suit_res6_all <- apply(suit_res7_all, 1, max)

    psi_all <- plogis(beta0 + beta1 * suit_res6_all)
    z_res6_all <- rbinom(M, 1, psi_all)

    theta_all <- plogis(gamma0 + gamma1 * suit_res7_all)
    z_res7_all <- matrix(0L, nrow = M, ncol = K_total)
    occ_idx <- which(z_res6_all == 1)
    for (i in occ_idx) {
      z_res7_all[i, ] <- rbinom(K_total, 1, theta_all[i, ])
    }

    true_N[sim] <- sum(z_res6_all) * 2

    # --- PPS sampling at Res 6 ---
    pi_res6 <- suit_res6_all / sum(suit_res6_all)
    sampled_idx <- sample.int(M, n_res6, replace = FALSE, prob = pi_res6)

    # --- Deterministic top-K at Res 7, vectorized ---
    # Pre-allocate observation arrays
    n_obs <- n_res6 * K_sampled
    obs_suit_res6 <- numeric(n_obs)
    obs_suit_res7 <- numeric(n_obs)
    obs_found   <- integer(n_obs)
    obs_t_det   <- numeric(n_obs)
    obs_cell_id <- integer(n_obs)

    for (ii in seq_along(sampled_idx)) {
      i <- sampled_idx[ii]
      top_k <- order(suit_res7_all[i, ], decreasing = TRUE)[1:K_sampled]
      rng <- ((ii - 1) * K_sampled + 1):(ii * K_sampled)

      obs_cell_id[rng] <- ii
      obs_suit_res6[rng] <- suit_res6_all[i]
      obs_suit_res7[rng] <- suit_res7_all[i, top_k]

      for (jj in seq_along(top_k)) {
        j <- top_k[jj]
        if (z_res6_all[i] == 1 && z_res7_all[i, j] == 1) {
          t_det <- rexp(1, rate = lambda_true)
          obs_found[rng[jj]] <- as.integer(t_det <= t_max)
          obs_t_det[rng[jj]] <- min(t_det, t_max)
        } else {
          obs_t_det[rng[jj]] <- t_max   # right-censored
        }
      }
    }

    # --- Organize by cell (using sequential cell IDs) ---
    cell_suit     <- tapply(obs_suit_res6, obs_cell_id, `[`, 1)
    cell_suit_res7 <- split(obs_suit_res7, obs_cell_id)
    cell_found    <- split(obs_found, obs_cell_id)
    cell_t        <- split(obs_t_det,  obs_cell_id)

    n_cells <- length(cell_suit)

    # --- Negative log-likelihood: TTD exponential hazard ---
    # params[5] = log(lambda); lambda = exponential hazard rate
    # Found sub-cell at time t:   log(theta) + log(lam) - lam*t
    # Censored sub-cell at t_max: log((1-theta) + theta*exp(-lam*t_max))
    neg_log_lik <- function(params) {
      b0 <- params[1]; b1 <- params[2]
      g0 <- params[3]; g1 <- params[4]
      log_lam  <- params[5]
      lam      <- exp(log_lam)
      surv_max <- exp(-lam * t_max)   # P(T > t_max) under Exp(lam)

      ll <- 0
      for (ci in 1:n_cells) {
        psi_i   <- plogis(b0 + b1 * cell_suit[ci])
        theta_j <- plogis(g0 + g1 * cell_suit_res7[[ci]])
        f_j     <- cell_found[[ci]]
        t_j     <- cell_t[[ci]]

        if (any(f_j > 0)) {
          ll_sub <- sum(ifelse(f_j > 0,
            log(pmax(theta_j, 1e-15)) + log_lam - lam * t_j,
            log(pmax((1 - theta_j) + theta_j * surv_max, 1e-15))
          ))
          ll <- ll + log(max(psi_i, 1e-15)) + ll_sub
        } else {
          ll_sub_zero <- sum(log(pmax((1 - theta_j) + theta_j * surv_max, 1e-15)))
          ll <- ll + log(max((1 - psi_i) + psi_i * exp(ll_sub_zero), 1e-15))
        }
      }
      return(-ll)
    }

    # --- Gradient (analytical for BFGS, TTD parameterisation) ---
    # params[5] = log(lambda); chain rule: d/d(log_lam) = d/d(lam) * lam
    neg_grad <- function(params) {
      b0 <- params[1]; b1 <- params[2]
      g0 <- params[3]; g1 <- params[4]
      log_lam  <- params[5]
      lam      <- exp(log_lam)
      surv_max <- exp(-lam * t_max)

      db0 <- 0; db1 <- 0; dg0 <- 0; dg1 <- 0; dlam <- 0

      for (ci in 1:n_cells) {
        psi_i   <- plogis(b0 + b1 * cell_suit[ci])
        theta_j <- plogis(g0 + g1 * cell_suit_res7[[ci]])
        f_j     <- cell_found[[ci]]
        t_j     <- cell_t[[ci]]
        s7      <- cell_suit_res7[[ci]]

        if (any(f_j > 0)) {
          # Cell is certainly occupied at Res 6
          db0 <- db0 + (1 - psi_i)
          db1 <- db1 + (1 - psi_i) * cell_suit[ci]

          for (k in seq_along(f_j)) {
            if (f_j[k] > 0) {
              # Sub-cell: nest found at t_j[k]
              # log-lik contribution: log(theta) + log_lam - lam*t_j
              dg0  <- dg0  + (1 - theta_j[k])
              dg1  <- dg1  + (1 - theta_j[k]) * s7[k]
              dlam <- dlam + (1 - lam * t_j[k])   # d/d(log_lam) = 1 - lam*t
            } else {
              # Sub-cell: right-censored at t_max
              # log-lik: log((1-theta) + theta*surv_max)
              p_cens   <- (1 - theta_j[k]) + theta_j[k] * surv_max
              dtheta_k <- theta_j[k] * (1 - theta_j[k])
              dg0  <- dg0  + (surv_max - 1) * dtheta_k / p_cens
              dg1  <- dg1  + (surv_max - 1) * dtheta_k * s7[k] / p_cens
              # d/d(log_lam): d(surv_max)/d(log_lam) = surv_max * (-lam*t_max)
              dlam <- dlam + theta_j[k] * surv_max * (-lam * t_max) / p_cens
            }
          }
        } else {
          # All sub-cells censored — integrate out psi
          p_cens_j  <- (1 - theta_j) + theta_j * surv_max
          ll_sub_zero <- sum(log(pmax(p_cens_j, 1e-15)))
          A     <- exp(ll_sub_zero)
          denom <- (1 - psi_i) + psi_i * A

          ratio_psi <- (A - 1) * psi_i * (1 - psi_i) / denom
          db0 <- db0 + ratio_psi
          db1 <- db1 + ratio_psi * cell_suit[ci]

          for (k in seq_along(f_j)) {
            dtheta_k    <- theta_j[k] * (1 - theta_j[k])
            dp_cens_dg  <- (surv_max - 1) * dtheta_k
            contrib_g   <- psi_i * A / p_cens_j[k] * dp_cens_dg / denom
            dg0  <- dg0  + contrib_g
            dg1  <- dg1  + contrib_g * s7[k]

            dp_cens_dlam <- theta_j[k] * surv_max * (-lam * t_max)
            dlam <- dlam + psi_i * A / p_cens_j[k] * dp_cens_dlam / denom
          }
        }
      }
      return(-c(db0, db1, dg0, dg1, dlam))
    }

    # --- Fit with BFGS ---
    start_log_lam <- log(lambda_true)   # params[5] = log(lambda)

    fit <- tryCatch({
      optim(c(-2, 3, 0, 2, start_log_lam), neg_log_lik, gr = neg_grad,
            method = "BFGS", control = list(maxit = 500), hessian = TRUE)
    }, error = function(e) {
      # Fallback to Nelder-Mead without hessian
      tryCatch({
        optim(c(-2, 3, 0, 2, start_log_lam), neg_log_lik,
              method = "Nelder-Mead", control = list(maxit = 5000), hessian = TRUE)
      }, error = function(e2) {
        list(par = c(-2, 3, 0, 2, start_log_lam), convergence = 1, hessian = NULL)
      })
    })

    converged[sim] <- (fit$convergence == 0)
    lambda_estimates[sim] <- exp(fit$par[5])

    b0_hat <- fit$par[1]
    b1_hat <- fit$par[2]
    psi_pred_all <- plogis(b0_hat + b1_hat * suit_res6_all)
    N_hat <- sum(psi_pred_all) * 2
    N_estimates[sim] <- N_hat

    # --- Delta method for SE(N_hat) and coverage ---
    ci_covers[sim] <- FALSE
    if (!is.null(fit$hessian) && fit$convergence == 0) {
      vcov <- tryCatch(solve(fit$hessian), error = function(e) NULL)
      if (!is.null(vcov) && all(is.finite(diag(vcov))) && all(diag(vcov)[1:2] > 0)) {
        dpsi_db0 <- psi_pred_all * (1 - psi_pred_all)
        dpsi_db1 <- psi_pred_all * (1 - psi_pred_all) * suit_res6_all
        grad_N <- 2 * c(sum(dpsi_db0), sum(dpsi_db1), 0, 0, 0)
        var_N <- as.numeric(t(grad_N) %*% vcov %*% grad_N)
        if (is.finite(var_N) && var_N > 0) {
          se_N <- sqrt(var_N)
          ci_lo <- N_hat - 1.96 * se_N
          ci_hi <- N_hat + 1.96 * se_N
          ci_covers[sim] <- (true_N[sim] >= ci_lo && true_N[sim] <= ci_hi)
        }
      }
    }

    # --- Null model: psi~1, theta~1, lambda~1 (TTD, intercept-only) ---
    neg_log_lik_null <- function(params) {
      a0      <- params[1]   # logit(psi)
      g0      <- params[2]   # logit(theta)
      log_lam <- params[3]   # log(lambda)
      psi      <- plogis(a0)
      theta    <- plogis(g0)
      lam      <- exp(log_lam)
      surv_max <- exp(-lam * t_max)

      log_psi  <- log(max(psi,  1e-15))
      log_th   <- log(max(theta, 1e-15))
      log_cens <- log(max((1 - theta) + theta * surv_max, 1e-15))

      ll <- 0
      for (ci in 1:n_cells) {
        f_j <- cell_found[[ci]]
        t_j <- cell_t[[ci]]
        if (any(f_j > 0)) {
          ll_sub <- sum(ifelse(f_j > 0,
            log_th + log_lam - lam * t_j,
            log_cens))
          ll <- ll + log_psi + ll_sub
        } else {
          ll <- ll + log(max((1 - psi) + psi * exp(K_sampled * log_cens), 1e-15))
        }
      }
      return(-ll)
    }

    fit_null <- tryCatch({
      optim(c(-2, 0, start_log_lam), neg_log_lik_null,
            method = "BFGS", control = list(maxit = 500))
    }, error = function(e) {
      tryCatch({
        optim(c(-2, 0, start_log_lam), neg_log_lik_null,
              method = "Nelder-Mead", control = list(maxit = 5000))
      }, error = function(e2) {
        list(par = c(-2, 0, start_log_lam), convergence = 1)
      })
    })

    N_estimates_null[sim] <- M * plogis(fit_null$par[1]) * 2

    # Write progress every 100 iterations
    if (sim %% 100 == 0) {
      writeLines(as.character(sim),
                 file.path(prog_dir, paste0("combo_", combo_idx, ".txt")))
    }
  }

  # Final progress
  writeLines(as.character(n_sims),
             file.path(prog_dir, paste0("combo_", combo_idx, ".txt")))

  # Convert estimated lambda back to per-visit detection probability
  # lambda is hazard per visit, so p_daily = 1 - exp(-lambda * 1)
  p_hat_vec <- 1 - exp(-lambda_estimates)
  data.frame(
    scenario = scenarios$scenario[s],
    p_true = p_true,
    n = n_res6,
    true_N = mean(true_N),
    mean_N_hat = mean(N_estimates),
    bias_pct = (mean(N_estimates) - mean(true_N)) / mean(true_N) * 100,
    bias_pct_null = (mean(N_estimates_null) - mean(true_N)) / mean(true_N) * 100,
    cv_pct = sd(N_estimates) / mean(N_estimates) * 100,
    mean_p_hat = mean(p_hat_vec),
    p_bias_pct = (mean(p_hat_vec) - p_true) / p_true * 100,
    pct_converged = mean(converged) * 100,
    coverage_95 = mean(ci_covers[converged]) * 100,
    stringsAsFactors = FALSE
  )
}

# --- Run in parallel ---
n_cores <- 6
cat(sprintf("Starting parallel run on %d cores...\n", n_cores))
cat(sprintf("Start time: %s\n\n", Sys.time()))
t_start <- Sys.time()

cl <- makeCluster(n_cores)
clusterExport(cl, c("param_grid", "scenarios", "suit_values", "M",
                     "K_total", "K_sampled", "J", "n_sims",
                     "gamma0", "gamma1", "prog_dir"))
clusterEvalQ(cl, { suppressPackageStartupMessages(library(stats)) })

results_list <- parLapply(cl, 1:nrow(param_grid), run_combo,
                          param_grid = param_grid,
                          scenarios = scenarios,
                          suit_values = suit_values,
                          M = M, K_total = K_total, K_sampled = K_sampled,
                          J = J, n_sims = n_sims,
                          gamma0 = gamma0, gamma1 = gamma1,
                          prog_dir = prog_dir)

stopCluster(cl)

results <- bind_rows(results_list)

t_total <- difftime(Sys.time(), t_start, units = "hours")

# Results
cat("\n================================================================\n")
cat("RESULTS (p estimated, 1000 simulations, with coverage)\n")
cat("================================================================\n\n")

cat("0. BIAS COMPARISON: NULL vs TWO-LEVEL MODEL\n")
cat("--------------------------------------------\n")
bias_compare <- results %>%
  group_by(scenario, p_true) %>%
  summarise(
    null_bias     = round(mean(bias_pct_null), 2),
    twolevel_bias = round(mean(bias_pct), 2),
    .groups = "drop"
  )
print(bias_compare)

cat("\n1. BIAS BY SCENARIO AND DETECTION\n")
cat("----------------------------------\n")
bias_table <- results %>%
  group_by(scenario, p_true) %>%
  summarise(mean_bias = round(mean(bias_pct), 1), .groups = "drop") %>%
  pivot_wider(names_from = p_true, values_from = mean_bias, names_prefix = "p=")
print(bias_table)

cat("\n2. CV BY SAMPLE SIZE AND SCENARIO (p = 0.10)\n")
cat("---------------------------------------------\n")
cv_table <- results %>%
  filter(p_true == 0.10) %>%
  select(scenario, n, cv_pct) %>%
  mutate(cv_pct = round(cv_pct, 1)) %>%
  pivot_wider(names_from = scenario, values_from = cv_pct)
print(cv_table)

cat("\n3. MINIMUM SAMPLE SIZE FOR CV <= 15%\n")
cat("-------------------------------------\n")
min_n <- results %>%
  filter(cv_pct <= 15) %>%
  group_by(scenario, p_true) %>%
  summarise(min_n = min(n), cv_at_min = round(min(cv_pct), 1), .groups = "drop") %>%
  pivot_wider(names_from = p_true, values_from = c(min_n, cv_at_min), names_sep = "_p=")
print(min_n)

cat("\n4. DETECTION PROBABILITY ESTIMATION\n")
cat("-------------------------------------\n")
p_table <- results %>%
  group_by(p_true) %>%
  summarise(
    mean_p_hat = round(mean(mean_p_hat), 4),
    p_bias_pct = round(mean(p_bias_pct), 1),
    .groups = "drop"
  )
print(p_table)

cat("\n5. 95% CI COVERAGE (nominal = 95%)\n")
cat("------------------------------------\n")
cov_table <- results %>%
  group_by(scenario, p_true) %>%
  summarise(mean_coverage = round(mean(coverage_95), 1), .groups = "drop") %>%
  pivot_wider(names_from = p_true, values_from = mean_coverage, names_prefix = "p=")
print(cov_table)

cat("\n   Coverage by sample size (p = 0.10):\n")
cov_n_table <- results %>%
  filter(p_true == 0.10) %>%
  select(scenario, n, coverage_95) %>%
  mutate(coverage_95 = round(coverage_95, 1)) %>%
  pivot_wider(names_from = scenario, values_from = coverage_95)
print(cov_n_table)

cat("\n6. CONVERGENCE RATES\n")
cat("---------------------\n")
conv_table <- results %>%
  group_by(scenario, p_true) %>%
  summarise(mean_conv = round(mean(pct_converged), 1), .groups = "drop") %>%
  pivot_wider(names_from = p_true, values_from = mean_conv, names_prefix = "p=")
print(conv_table)

# Save
write_csv(results, "C:/Users/chrim/Dropbox/Current projects/PHEA/PHEA Surveys/PHEA H3/power_analysis_two_level_EstP.csv")
cat("\nResults saved to: power_analysis_two_level_EstP.csv\n")

# Figure — CV
fig_data <- results %>%
  mutate(p_label = paste0("p = ", p_true))

p1 <- ggplot(fig_data, aes(x = n, y = cv_pct, color = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "red") +
  facet_wrap(~p_label) +
  scale_color_brewer(name = "Scenario", palette = "Dark2") +
  scale_x_continuous(breaks = seq(150, 400, 50)) +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  labs(
    x = "Number of Res 6 cells surveyed",
    y = "Coefficient of Variation (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("C:/Users/chrim/Dropbox/Current projects/PHEA/PHEA Surveys/PHEA H3/Figure_CV_Two_Level_EstP.png",
       p1, width = 10, height = 6, dpi = 300)
cat("Figure saved: Figure_CV_Two_Level_EstP.png\n")

# Figure — Coverage
p2 <- ggplot(fig_data, aes(x = n, y = coverage_95, color = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
  facet_wrap(~p_label) +
  scale_color_brewer(name = "Scenario", palette = "Dark2") +
  scale_x_continuous(breaks = seq(150, 400, 50)) +
  scale_y_continuous(limits = c(50, 100), breaks = seq(50, 100, 10)) +
  labs(
    x = "Number of Res 6 cells surveyed",
    y = "95% CI Coverage (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("C:/Users/chrim/Dropbox/Current projects/PHEA/PHEA Surveys/PHEA H3/Figure_Coverage_Two_Level_EstP.png",
       p2, width = 10, height = 6, dpi = 300)
cat("Figure saved: Figure_Coverage_Two_Level_EstP.png\n")

cat(sprintf("\nTotal runtime: %.1f hours\n", as.numeric(t_total)))
cat("\n================================================================\n")
cat("COMPLETE\n")
cat("================================================================\n")
