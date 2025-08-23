#load packages ----------------------------------------------------------
library(tidyverse)
library(brms)
library(sn)
library(lme4)
library(arm)
--------------------------------------------------------------------------------
#model misspecification----
--------------------------------------------------------------------------------
##no theta in sim ----

run_no_theta_sim <- function(n_items, n_participants = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sigma_beta <- 1
  mean_beta <- 0
  thresholds <- c(0.8, 1.8)
  beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)
  
  # Create data
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(
      participant = 1:n_participants,
      item = 1:n_items
    )
  }
  
  sim_data <- sim_data |> 
    rowwise() |>
    mutate(
      theta = 0,
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_inc  = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_id   = 1 - pnorm(thresholds[2] - eta),
      y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
    ) |>
    ungroup() |>
    mutate(
      y = factor(y, levels = 1:3, labels = c("excl", "inc", "id"), ordered = TRUE)
    )
  
  # Fit model
  start_time <- Sys.time()
  brm_model <- brm(
    formula = y ~ (1 | participant) + (1 | item),
    data = sim_data,
    family = cumulative(probit),
    chains = 4,
    iter = 2000,
    silent = TRUE,
    refresh = 0
  )
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Extract summaries
  ranefs <- ranef(brm_model)
  
  theta_df <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
    mutate(participant = as.integer(participant)) |>
    rename(theta_estimate = Estimate, theta_lower = Q2.5, theta_upper = Q97.5)
  
  beta_df <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
    mutate(item = as.integer(item)) |>
    rename(beta_estimate = Estimate, beta_lower = Q2.5, beta_upper = Q97.5)
  
  sim_data <- sim_data |>
    left_join(theta_df, by = "participant") |>
    left_join(beta_df, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
}


# Run simulations
item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_no_theta_sim <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(.l = list(n_items = n_items, rep = rep), function(n_items, rep) {
    run_no_theta_sim(n_items = n_items, seed = 1000 + rep)
  }))

# --- Per-run results ---
results_no_theta_sim_per_run <- results_no_theta_sim |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    simulation = "No theta",
    model = "Theta"
  )

# --- Summary (average over reps) ---
results_no_theta_sim_summary <- results_no_theta_sim_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    simulation = first(simulation),
    model = first(model),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_no_theta_sim_per_run, "sim_results_no_theta_sim_per_run.csv")
write_csv(results_no_theta_sim_summary, "sim_results_no_theta_sim_summary.csv")
--------------------------------------------------------------------------------
##no theta in model----

run_no_theta_model <- function(n_items, n_participants = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sigma_theta <- 1
  sigma_beta <- 1 
  mean_beta <- 0
  mean_theta <- 0 
  thresholds <- c(0.8, 1.8)
  
  theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
  beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)
  
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(
      participant = 1:n_participants,
      item = 1:n_items
    )
  }
  
  sim_data <- sim_data |> 
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_id = 1 - pnorm(thresholds[2] - eta),
      y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
    ) |> 
    ungroup() |>
    mutate(
      y = factor(y, levels = 1:3, labels = c("excl", "inc", "id"), ordered = TRUE)
    )
  
  start_time <- Sys.time()
  brm_model <- brm(
    formula = y ~ (1 | item),
    data = sim_data,
    family = cumulative(probit),
    chains = 4,
    iter = 2000,
    silent = TRUE,
    refresh = 0
  )
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ranefs <- ranef(brm_model)
  
  item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
    mutate(item = as.integer(item)) |>
    rename(beta_estimate = Estimate, beta_lower = Q2.5, beta_upper = Q97.5)
  
  sim_data <- sim_data |>
    left_join(item_summary, by = "item") |>
    mutate(
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  list(
    time_to_fit = run_time,
    theta_coverage = NA,
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = NA,
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_no_theta_model <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(.l = list(n_items = n_items, rep = rep), function(n_items, rep) {
    run_no_theta_model(n_items = n_items)
  }))

# Per-run results
results_no_theta_model_per_run <- results_no_theta_model |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    simulation = "Theta",
    model = "No theta"
  )

# Summary (average over reps)
results_no_theta_model_summary <- results_no_theta_model_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    simulation = first(simulation),
    model = first(model),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = NA_real_,
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = NA_real_,
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_no_theta_model_per_run, "sim_results_no_theta_model_per_run.csv")
write_csv(results_no_theta_model_summary, "sim_results_no_theta_model_summary.csv")
--------------------------------------------------------------------------------
##skewed theta----

run_skewed_theta <- function(n_items, n_participants = 100, alpha_theta = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sigma_theta <- 1
  sigma_beta <- 1 
  mean_beta <- 0
  mean_theta <- 0 
  thresholds <- c(.8, 1.8)
  
  # Convert skew-normal params to mean/sd
  delta_theta <- alpha_theta / sqrt(1 + alpha_theta^2)
  omega_theta <- sqrt((sigma_theta^2) / (1 - (2 * delta_theta^2) / pi))
  xi_theta <- -(omega_theta * delta_theta * sqrt(2 / pi))
  
  theta <- rsn(n_participants, xi = xi_theta, omega = omega_theta, alpha = alpha_theta)
  beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)
  
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(
      participant = 1:n_participants,
      item = 1:n_items
    )
  }
  
  sim_data <- sim_data |> 
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_id = 1 - pnorm(thresholds[2] - eta),
      y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
    ) |>
    ungroup() |>
    mutate(
      y = factor(y, levels = 1:3, labels = c("excl", "inc", "id"), ordered = TRUE)
    )
  
  start_time <- Sys.time()
  brm_model <- brm(
    formula = y ~ (1 | participant) + (1 | item),
    data = sim_data,
    family = cumulative(probit),
    chains = 4,
    iter = 2000,
    silent = TRUE,
    refresh = 0
  )
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ranefs <- ranef(brm_model)
  
  participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
    mutate(participant = as.integer(participant)) |>
    rename(theta_estimate = Estimate, theta_lower = Q2.5, theta_upper = Q97.5)
  
  item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
    mutate(item = as.integer(item)) |>
    rename(beta_estimate = Estimate, beta_lower = Q2.5, beta_upper = Q97.5)
  
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_skewed_theta <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(.l = list(n_items = n_items, rep = rep), function(n_items, rep) {
    run_skewed_theta(n_items = n_items)
  }))

# Per-run results
results_skewed_theta_per_run <- results_skewed_theta |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    simulation = "Skewed theta",
    model = "Normal theta"
  )

# Summary (average over reps)
results_skewed_theta_summary <- results_skewed_theta_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    simulation = first(simulation),
    model = first(model),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_skewed_theta_per_run, "sim_results_skewed_theta_per_run.csv")
write_csv(results_skewed_theta_summary, "sim_results_skewed_theta_summary.csv")
--------------------------------------------------------------------------------
##skewed beta----

run_skewed_beta <- function(n_items, n_participants = 100, alpha_beta = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sigma_theta <- 1
  sigma_beta <- 1 
  mean_beta <- 0
  mean_theta <- 0 
  thresholds <- c(0.8, 1.8)
  
  # Convert skew-normal parameters for beta
  delta_beta <- alpha_beta / sqrt(1 + alpha_beta^2)
  omega_beta <- sqrt((sigma_beta^2) / (1 - (2 * delta_beta^2) / pi))
  xi_beta <- -(omega_beta * delta_beta * sqrt(2 / pi))
  
  theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
  beta <- rsn(n_items, xi = xi_beta, omega = omega_beta, alpha = alpha_beta)
  
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(
      participant = 1:n_participants,
      item = 1:n_items
    )
  }
  
  sim_data <- sim_data |>
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_id = 1 - pnorm(thresholds[2] - eta),
      y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
    ) |>
    ungroup() |>
    mutate(
      y = factor(y, levels = 1:3, labels = c("excl", "inc", "id"), ordered = TRUE)
    )
  
  start_time <- Sys.time()
  brm_model <- brm(
    formula = y ~ (1 | participant) + (1 | item),
    data = sim_data,
    family = cumulative(probit),
    chains = 4,
    iter = 2000,
    silent = TRUE,
    refresh = 0
  )
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ranefs <- ranef(brm_model)
  
  participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
    mutate(participant = as.integer(participant)) |>
    rename(theta_estimate = Estimate, theta_lower = Q2.5, theta_upper = Q97.5)
  
  item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
    mutate(item = as.integer(item)) |>
    rename(beta_estimate = Estimate, beta_lower = Q2.5, beta_upper = Q97.5)
  
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_skewed_beta <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(.l = list(n_items = n_items, rep = rep), function(n_items, rep) {
    run_skewed_beta(n_items = n_items)
  }))

# Per-run results
results_skewed_beta_per_run <- results_skewed_beta |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    simulation = "Skewed beta",
    model = "Normal beta"
  )

# Summary (average over reps)
results_skewed_beta_summary <- results_skewed_beta_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    simulation = first(simulation),
    model = first(model),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_skewed_beta_per_run, "sim_results_skewed_beta_per_run.csv")
write_csv(results_skewed_beta_summary, "sim_results_skewed_beta_summary.csv")
--------------------------------------------------------------------------------
##both skewed----

run_both_skewed <- function(n_items, n_participants = 100, alpha_theta = 5, alpha_beta = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sigma_theta <- 1
  sigma_beta <- 1 
  mean_beta <- 0
  mean_theta <- 0 
  thresholds <- c(0.8, 1.8)
  
  # Skew-normal parameters
  delta_theta <- alpha_theta / sqrt(1 + alpha_theta^2)
  omega_theta <- sqrt((sigma_theta^2) / (1 - (2 * delta_theta^2) / pi))
  xi_theta <- -(omega_theta * delta_theta * sqrt(2 / pi))
  
  delta_beta <- alpha_beta / sqrt(1 + alpha_beta^2)
  omega_beta <- sqrt((sigma_beta^2) / (1 - (2 * delta_beta^2) / pi))
  xi_beta <- -(omega_beta * delta_beta * sqrt(2 / pi))
  
  theta <- rsn(n_participants, xi = xi_theta, omega = omega_theta, alpha = alpha_theta)
  beta <- rsn(n_items, xi = xi_beta, omega = omega_beta, alpha = alpha_beta)
  
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(participant = 1:n_participants, item = 1:n_items)
  }
  
  sim_data <- sim_data |>
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_id = 1 - pnorm(thresholds[2] - eta),
      y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
    ) |>
    ungroup() |>
    mutate(y = factor(y, levels = 1:3, labels = c("excl", "inc", "id"), ordered = TRUE))
  
  start_time <- Sys.time()
  brm_model <- brm(
    formula = y ~ (1 | participant) + (1 | item),
    data = sim_data,
    family = cumulative(probit),
    chains = 4,
    iter = 2000,
    silent = TRUE,
    refresh = 0
  )
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ranefs <- ranef(brm_model)
  
  participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
    mutate(participant = as.integer(participant)) |>
    rename(theta_estimate = Estimate, theta_lower = Q2.5, theta_upper = Q97.5)
  
  item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
    mutate(item = as.integer(item)) |>
    rename(beta_estimate = Estimate, beta_lower = Q2.5, beta_upper = Q97.5)
  
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_both_skewed <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(.l = list(n_items = n_items, rep = rep), function(n_items, rep) {
    run_both_skewed(n_items = n_items)
  }))

# Per-run results
results_both_skewed_per_run <- results_both_skewed |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    simulation = "Both skewed",
    model = "Both normal"
  )

# Summary (average over reps)
results_both_skewed_summary <- results_both_skewed_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    simulation = first(simulation),
    model = first(model),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_both_skewed_per_run, "sim_results_both_skewed_per_run.csv")
write_csv(results_both_skewed_summary, "sim_results_both_skewed_summary.csv")
--------------------------------------------------------------------------------
#mixture ----

run_mixture <- function(
    n_items, n_participants = 100,
    alpha_theta = 5, alpha_beta = 5,
    p_skew_theta = 1/3, p_skew_beta = 1/3,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  sigma_theta <- 1
  sigma_beta <- 1 
  mean_theta <- 0
  mean_beta <- 0
  thresholds <- c(0.8, 1.8)
  
  # --- theta ---
  n_skew_theta <- round(n_participants * p_skew_theta)
  n_norm_theta <- n_participants - n_skew_theta
  
  # Skew-normal parameters for theta
  delta_theta <- alpha_theta / sqrt(1 + alpha_theta^2)
  omega_theta <- sqrt((sigma_theta^2) / (1 - (2 * delta_theta^2) / pi))
  xi_theta <- -(omega_theta * delta_theta * sqrt(2 / pi))
  
  theta <- c(
    rsn(n_skew_theta, xi = xi_theta, omega = omega_theta, alpha = alpha_theta),
    rnorm(n_norm_theta, mean = mean_theta, sd = sigma_theta)
  )
  theta <- sample(theta)
  
  # --- beta ---
  n_skew_beta <- round(n_items * p_skew_beta)
  n_norm_beta <- n_items - n_skew_beta
  
  delta_beta <- alpha_beta / sqrt(1 + alpha_beta^2)
  omega_beta <- sqrt((sigma_beta^2) / (1 - (2 * delta_beta^2) / pi))
  xi_beta <- -(omega_beta * delta_beta * sqrt(2 / pi))
  
  beta <- c(
    rsn(n_skew_beta, xi = xi_beta, omega = omega_beta, alpha = alpha_beta),
    rnorm(n_norm_beta, mean = mean_beta, sd = sigma_beta)
  )
  beta <- sample(beta)
  
  # --- data ---
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(participant = 1:n_participants, item = 1:n_items)
  }
  
  sim_data <- sim_data |>
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_id = 1 - pnorm(thresholds[2] - eta),
      y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
    ) |>
    ungroup() |>
    mutate(y = factor(y, levels = 1:3, labels = c("excl", "inc", "id"), ordered = TRUE))
  
  # --- fit model ---
  start_time <- Sys.time()
  brm_model <- brm(
    formula = y ~ (1 | participant) + (1 | item),
    data = sim_data,
    family = cumulative(probit),
    chains = 4,
    iter = 2000,
    silent = TRUE,
    refresh = 0
  )
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ranefs <- ranef(brm_model)
  
  participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
    mutate(participant = as.integer(participant)) |>
    rename(theta_estimate = Estimate, theta_lower = Q2.5, theta_upper = Q97.5)
  
  item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
    mutate(item = as.integer(item)) |>
    rename(beta_estimate = Estimate, beta_lower = Q2.5, beta_upper = Q97.5)
  
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

# Run simulations
results_mixture <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(.l = list(n_items = n_items, rep = rep), function(n_items, rep) {
    run_mixture(n_items = n_items)
  }))

# Per-run results
results_mixture_per_run <- results_mixture |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    simulation = "Mixture",
    model = "Both normal"
  )

# Summary (average over reps)
results_mixture_summary <- results_mixture_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    simulation = first(simulation),
    model = first(model),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_mixture_per_run, "sim_results_mixture_per_run.csv")
write_csv(results_mixture_summary, "sim_results_mixture_summary.csv")
--------------------------------------------------------------------------------
#combine prior misspec csvs----

#categories per run
prior_misspec_per_run_repeated <- sim_results_no_theta_sim_per_run |> 
  bind_rows(sim_results_no_theta_model_per_run) |> 
  bind_rows(sim_results_skewed_theta_per_run) |> 
  bind_rows(sim_results_skewed_beta_per_run) |> 
  bind_rows(sim_results_both_skewed_per_run) |> 
  bind_rows(sim_results_mixture_per_run)

# update the values
prior_misspec_per_run_repeated <- prior_misspec_per_run_repeated %>%
  mutate(
    simulation = recode(simulation,
                        "Partial skew" = "Mixture"))


prior_misspec_per_run_repeated <- prior_misspec_per_run_repeated %>%
  mutate(
    condition = paste(simulation, model, sep = " - ")
  )

prior_misspec_per_run_repeated <- prior_misspec_per_run_repeated %>%
  mutate(
    condition = case_when(
      simulation == "No theta" & model == "Theta" ~ "No theta sim",
      simulation == "Theta" & model == "No theta" ~ "No theta model",
      simulation == "Skewed theta" & model == "Normal theta" ~ "Skewed theta",
      simulation == "Skewed beta" & model == "Normal beta" ~ "Skewed beta",
      simulation == "Both skewed" & model == "Both normal" ~ "Both skewed",
      simulation == "Mixture" & model == "Mixed normal/skew" ~ "Mixture"
    )
  )

prior_misspec_per_run_repeated <- prior_misspec_per_run_repeated %>%
  mutate(
    theta_correlation = coalesce(theta_correlation, theta_corr),
    beta_correlation  = coalesce(beta_correlation, beta_corr)
  ) %>%
<<<<<<< Updated upstream
  select(-theta_corr, -beta_corr)   # drop the extras


write_csv(prior_misspec_per_run_repeated, "prior_misspec_per_run_repeated.csv")


#categories summary
prior_misspec_summary_repeated <- sim_results_no_theta_sim_summary |> 
  bind_rows(sim_results_no_theta_model_summary) |> 
  bind_rows(sim_results_skewed_theta_summary) |> 
  bind_rows(sim_results_skewed_beta_summary) |> 
  bind_rows(sim_results_both_skewed_summary) |> 
  bind_rows(sim_results_mixture_summary)


# update the values
prior_misspec_summary_repeated <- prior_misspec_summary_repeated %>%
  mutate(
    simulation = recode(simulation,
                        "Partial skew" = "Mixture"),
    model = recode(model,
                   "Mixed normal/skew" = "Both normal")
  )


prior_misspec_summary_repeated <- prior_misspec_summary_repeated %>%
  mutate(
    condition = paste(simulation, model, sep = " - ")
  )

prior_misspec_summary_repeated <- prior_misspec_summary_repeated %>%
  mutate(
    condition = case_when(
      simulation == "No theta" & model == "Theta" ~ "No theta sim",
      simulation == "Theta" & model == "No theta" ~ "No theta model",
      simulation == "Skewed theta" & model == "Normal theta" ~ "Skewed theta",
      simulation == "Skewed beta" & model == "Normal beta" ~ "Skewed beta",
      simulation == "Both skewed" & model == "Both normal" ~ "Both skewed",
      simulation == "Mixture" & model == "Both normal" ~ "Mixture"
    )
  )

write_csv(prior_misspec_summary_repeated, "prior_misspec_summary_repeated.csv")
=======
  dplyr::select(-beta_corr, -theta_corr)

write_csv(prior_misspec_per_run_repeated, "prior_misspec_per_run_repeated.csv")
>>>>>>> Stashed changes
--------------------------------------------------------------------------------
#3 vs 5 category----
--------------------------------------------------------------------------------
##3 cat ----

run_3cat_sim <- function(n_items, n_participants = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sigma_theta <- 1
  sigma_beta <- 1
  mean_theta <- 0
  mean_beta <- 0
  thresholds <- c(0.8, 1.8)
  
  theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
  beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)
  
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(participant = 1:n_participants, item = 1:n_items)
  }
  
  sim_data <- sim_data |>
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_id = 1 - pnorm(thresholds[2] - eta),
      y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
    ) |>
    ungroup() |>
    mutate(y = factor(y, levels = 1:3, labels = c("excl", "inc", "id"), ordered = TRUE))
  
  start_time <- Sys.time()
  brm_model <- brm(
    formula = y ~ (1 | participant) + (1 | item),
    data = sim_data,
    family = cumulative(probit),
    chains = 4,
    iter = 2000,
    silent = TRUE,
    refresh = 0
  )
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ranefs <- ranef(brm_model)
  
  participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
    mutate(participant = as.integer(participant)) |>
    rename(theta_estimate = Estimate, theta_lower = Q2.5, theta_upper = Q97.5)
  
  item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
    mutate(item = as.integer(item)) |>
    rename(beta_estimate = Estimate, beta_lower = Q2.5, beta_upper = Q97.5)
  
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_3cat <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(.l = list(n_items = n_items, rep = rep), function(n_items, rep) {
    run_3cat_sim(n_items = n_items)
  }))

# Per-run results
results_3cat_per_run <- results_3cat |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    category_system = "3-category"
  )

# Summary (averaged over reps)
results_3cat_summary <- results_3cat_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    category_system = first(category_system),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_3cat_per_run, "sim_results_3cat_per_run.csv")
write_csv(results_3cat_summary, "sim_results_3cat_summary.csv")
--------------------------------------------------------------------------------
##5 cat (a) ----

run_5cat_a_sim <- function(n_items, n_participants = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sigma_theta <- 1
  sigma_beta <- 1
  mean_theta <- 0
  mean_beta <- 0
  thresholds <- c(0.8, 1.1, 1.5, 1.8)
  
  theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
  beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)
  
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(participant = 1:n_participants, item = 1:n_items)
  }
  
  sim_data <- sim_data |>
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_lean_excl = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[3] - eta) - pnorm(thresholds[2] - eta),
      p_lean_id = pnorm(thresholds[4] - eta) - pnorm(thresholds[3] - eta),
      p_id = 1 - pnorm(thresholds[4] - eta),
      y = which.max(rmultinom(1, 1, c(p_excl, p_lean_excl, p_inc, p_lean_id, p_id)))
    ) |>
    ungroup() |>
    mutate(y = factor(y, levels = 1:5, labels = c("excl", "lean_excl", "inc", "lean_id", "id"), ordered = TRUE))
  
  start_time <- Sys.time()
  brm_model <- brm(
    formula = y ~ (1 | participant) + (1 | item),
    data = sim_data,
    family = cumulative(probit),
    chains = 4,
    iter = 2000,
    silent = TRUE,
    refresh = 0
  )
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ranefs <- ranef(brm_model)
  
  participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
    mutate(participant = as.integer(participant)) |>
    rename(theta_estimate = Estimate, theta_lower = Q2.5, theta_upper = Q97.5)
  
  item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
    mutate(item = as.integer(item)) |>
    rename(beta_estimate = Estimate, beta_lower = Q2.5, beta_upper = Q97.5)
  
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_5cat_a <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(.l = list(n_items = n_items, rep = rep), function(n_items, rep) {
    run_5cat_a_sim(n_items = n_items)
  }))

# Per-run results
results_5cat_a_per_run <- results_5cat_a |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    category_system = "5-category A"
  )

# Summary (averaged over reps)
results_5cat_a_summary <- results_5cat_a_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    category_system = first(category_system),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_5cat_a_per_run, "sim_results_5cat_a_per_run.csv")
write_csv(results_5cat_a_summary, "sim_results_5cat_a_summary.csv")
--------------------------------------------------------------------------------
##5 cat (b) ----

run_5cat_b_sim <- function(n_items, n_participants = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sigma_theta <- 1
  sigma_beta <- 1
  mean_theta <- 0
  mean_beta <- 0
  thresholds <- c(0.6, 0.9, 1.5, 2.0)
  
  theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
  beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)
  
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(participant = 1:n_participants, item = 1:n_items)
  }
  
  sim_data <- sim_data |>
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_lean_excl = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[3] - eta) - pnorm(thresholds[2] - eta),
      p_lean_id = pnorm(thresholds[4] - eta) - pnorm(thresholds[3] - eta),
      p_id = 1 - pnorm(thresholds[4] - eta),
      y = which.max(rmultinom(1, 1, c(p_excl, p_lean_excl, p_inc, p_lean_id, p_id)))
    ) |>
    ungroup() |>
    mutate(y = factor(y, levels = 1:5, labels = c("excl", "lean_excl", "inc", "lean_id", "id"), ordered = TRUE))
  
  start_time <- Sys.time()
  brm_model <- brm(
    formula = y ~ (1 | participant) + (1 | item),
    data = sim_data,
    family = cumulative(probit),
    chains = 4,
    iter = 2000,
    silent = TRUE,
    refresh = 0
  )
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ranefs <- ranef(brm_model)
  
  participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
    mutate(participant = as.integer(participant)) |>
    rename(theta_estimate = Estimate, theta_lower = Q2.5, theta_upper = Q97.5)
  
  item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
    mutate(item = as.integer(item)) |>
    rename(beta_estimate = Estimate, beta_lower = Q2.5, beta_upper = Q97.5)
  
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_5cat_b <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(.l = list(n_items = n_items, rep = rep), function(n_items, rep) {
    run_5cat_b_sim(n_items = n_items)
  }))

# Per-run detailed results
results_5cat_b_per_run <- results_5cat_b |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    category_system = "5-category B"
  )

# Summary (averaged over reps)
results_5cat_b_summary <- results_5cat_b_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    category_system = first(category_system),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_5cat_b_per_run, "sim_results_5cat_b_per_run.csv")
write_csv(results_5cat_b_summary, "sim_results_5cat_b_summary.csv")
--------------------------------------------------------------------------------
#combine cat csvs----

categories_per_run_repeated <- sim_results_3cat_per_run |> 
  bind_rows(sim_results_5cat_a_per_run) |> 
  bind_rows(sim_results_5cat_b_per_run)

write_csv(categories_per_run_repeated, "categories_per_run_repeated.csv")

categories_summary_repeated <- sim_results_3cat_summary |> 
  bind_rows(sim_results_5cat_a_summary) |> 
  bind_rows(sim_results_5cat_b_summary)

write_csv(categories_summary_repeated, "categories_summary_repeated.csv")
--------------------------------------------------------------------------------
#Non-bayes ----
--------------------------------------------------------------------------------
##3 cat frequentist ----

run_frequentist_3 <- function(n_items, n_participants = 100, seed = NULL, return_data = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  # Parameters
  sigma_theta <- 1
  sigma_beta <- 1
  mean_theta <- 0
  mean_beta <- 0
  thresholds <- c(0.8, 1.8)
  
  # Latent traits
  theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
  beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)
  
  # Build dataset
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(participant = 1:n_participants, item = 1:n_items)
  }
  
  # Simulate ordinal responses
  sim_data <- sim_data |>
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_id = 1 - pnorm(thresholds[2] - eta),
      y_num = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
    ) |>
    ungroup()
  
  # Fit linear mixed model treating numeric categories as outcome
  start_time <- Sys.time()
  lmer_model <- lmer(y_num ~ 1 + (1 | participant) + (1 | item), data = sim_data)
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Extract random effects and SEs
  re <- ranef(lmer_model)
  se_re <- se.ranef(lmer_model)
  
  participant_summary <- tibble(
    participant = as.integer(rownames(re$participant)),
    theta_estimate = as.numeric(re$participant[,1]),
    theta_lower = theta_estimate - 1.96 * se_re$participant[,1],
    theta_upper = theta_estimate + 1.96 * se_re$participant[,1]
  )
  
  item_summary <- tibble(
    item = as.integer(rownames(re$item)),
    beta_estimate = as.numeric(re$item[,1]),
    beta_lower = beta_estimate - 1.96 * se_re$item[,1],
    beta_upper = beta_estimate + 1.96 * se_re$item[,1]
  )
  
  # Merge with true values and compute coverage
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  # Summary metrics
  result_list <- list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
  
  if (return_data) result_list$sim_data <- sim_data
  return(result_list)
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_frequentist_3 <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(list(n_items = n_items, rep = rep),
                       function(n_items, rep) run_frequentist_3(n_items)))

# Per-run results
results_frequentist_3_per_run <- results_frequentist_3 |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    category_system = "3-category-freq"
  )

# Summary over reps
results_frequentist_3_summary <- results_frequentist_3_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    category_system = first(category_system),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_frequentist_3_per_run, "sim_results_frequentist_3_per_run.csv")
write_csv(results_frequentist_3_summary, "sim_results_frequentist_3_summary.csv")
--------------------------------------------------------------------------------
##5 cat a frequentist ----

run_frequentist_5a <- function(n_items, n_participants = 100, seed = NULL, return_data = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  # Parameters
  sigma_theta <- 1
  sigma_beta <- 1
  mean_theta <- 0
  mean_beta <- 0
  thresholds <- c(0.8, 1.1, 1.5, 1.8)
  
  # Latent traits
  theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
  beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)
  
  # Build dataset
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(participant = 1:n_participants, item = 1:n_items)
  }
  
  # Simulate ordinal responses
  sim_data <- sim_data |>
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_lean_excl = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[3] - eta) - pnorm(thresholds[2] - eta),
      p_lean_id = pnorm(thresholds[4] - eta) - pnorm(thresholds[3] - eta),
      p_id = 1 - pnorm(thresholds[4] - eta),
      y_num = which.max(rmultinom(1, 1, c(p_excl, p_lean_excl, p_inc, p_lean_id, p_id)))
    ) |>
    ungroup()
  
  # Fit linear mixed model treating numeric categories as outcome
  start_time <- Sys.time()
  lmer_model <- lmer(y_num ~ 1 + (1 | participant) + (1 | item), data = sim_data)
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Extract random effects and SEs
  re <- ranef(lmer_model)
  se_re <- se.ranef(lmer_model)
  
  participant_summary <- tibble(
    participant = as.integer(rownames(re$participant)),
    theta_estimate = as.numeric(re$participant[,1]),
    theta_lower = theta_estimate - 1.96 * se_re$participant[,1],
    theta_upper = theta_estimate + 1.96 * se_re$participant[,1]
  )
  
  item_summary <- tibble(
    item = as.integer(rownames(re$item)),
    beta_estimate = as.numeric(re$item[,1]),
    beta_lower = beta_estimate - 1.96 * se_re$item[,1],
    beta_upper = beta_estimate + 1.96 * se_re$item[,1]
  )
  
  # Merge with true values and compute coverage
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  # Summary metrics
  result_list <- list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
  
  if (return_data) result_list$sim_data <- sim_data
  return(result_list)
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_frequentist_5a <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(list(n_items = n_items, rep = rep),
                       function(n_items, rep) run_frequentist_5a(n_items)))

# Per-run results
results_frequentist_5a_per_run <- results_frequentist_5a |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    category_system = "5-category-a-freq"
  )

# Summary over reps
results_frequentist_5a_summary <- results_frequentist_5a_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    category_system = first(category_system),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_frequentist_5a_per_run, "sim_results_frequentist_5a_per_run.csv")
write_csv(results_frequentist_5a_summary, "sim_results_frequentist_5a_summary.csv")
--------------------------------------------------------------------------------
##5 cat b frequentist ----

run_frequentist_5b <- function(n_items, n_participants = 100, seed = NULL, return_data = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  # Parameters
  sigma_theta <- 1
  sigma_beta <- 1
  mean_theta <- 0
  mean_beta <- 0
  thresholds <- c(0.6, 0.9, 1.5, 2.0)
  
  # Latent traits
  theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
  beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)
  
  # Build dataset
  if (n_items > 100) {
    sim_data <- tibble(
      participant = rep(1:n_participants, each = 100),
      item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
    )
  } else {
    sim_data <- expand_grid(participant = 1:n_participants, item = 1:n_items)
  }
  
  # Simulate ordinal responses
  sim_data <- sim_data |>
    rowwise() |>
    mutate(
      theta = theta[participant],
      beta = beta[item],
      eta = theta + beta,
      p_excl = pnorm(thresholds[1] - eta),
      p_lean_excl = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
      p_inc = pnorm(thresholds[3] - eta) - pnorm(thresholds[2] - eta),
      p_lean_id = pnorm(thresholds[4] - eta) - pnorm(thresholds[3] - eta),
      p_id = 1 - pnorm(thresholds[4] - eta),
      y_num = which.max(rmultinom(1, 1, c(p_excl, p_lean_excl, p_inc, p_lean_id, p_id)))
    ) |>
    ungroup()
  
  # Fit linear mixed model treating numeric categories as outcome
  start_time <- Sys.time()
  lmer_model <- lmer(y_num ~ 1 + (1 | participant) + (1 | item), data = sim_data)
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Extract random effects and SEs
  re <- ranef(lmer_model)
  se_re <- se.ranef(lmer_model)
  
  participant_summary <- tibble(
    participant = as.integer(rownames(re$participant)),
    theta_estimate = as.numeric(re$participant[,1]),
    theta_lower = theta_estimate - 1.96 * se_re$participant[,1],
    theta_upper = theta_estimate + 1.96 * se_re$participant[,1]
  )
  
  item_summary <- tibble(
    item = as.integer(rownames(re$item)),
    beta_estimate = as.numeric(re$item[,1]),
    beta_lower = beta_estimate - 1.96 * se_re$item[,1],
    beta_upper = beta_estimate + 1.96 * se_re$item[,1]
  )
  
  # Merge with true values and compute coverage
  sim_data <- sim_data |>
    left_join(participant_summary, by = "participant") |>
    left_join(item_summary, by = "item") |>
    mutate(
      theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
      beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper)
    )
  
  # Summary metrics
  result_list <- list(
    time_to_fit = run_time,
    theta_coverage = mean(sim_data$theta_covered, na.rm = TRUE),
    beta_coverage = mean(sim_data$beta_covered, na.rm = TRUE),
    theta_correlation = cor(sim_data$theta, sim_data$theta_estimate, use = "complete.obs"),
    beta_correlation = cor(sim_data$beta, sim_data$beta_estimate, use = "complete.obs")
  )
  
  if (return_data) result_list$sim_data <- sim_data
  return(result_list)
}

item_sizes <- c(20, 50, 300, 600)
n_reps <- 10

results_frequentist_5b <- expand_grid(n_items = item_sizes, rep = 1:n_reps) |>
  mutate(result = pmap(list(n_items = n_items, rep = rep),
                       function(n_items, rep) run_frequentist_5b(n_items)))

# Per-run results
results_frequentist_5b_per_run <- results_frequentist_5b |>
  unnest_wider(result) |>
  mutate(
    n_participants = 100,
    category_system = "5-category-b-freq"
  )

# Summary over reps
results_frequentist_5b_summary <- results_frequentist_5b_per_run |>
  group_by(n_items) |>
  summarize(
    n_participants = first(n_participants),
    category_system = first(category_system),
    time_to_fit = mean(time_to_fit, na.rm = TRUE),
    theta_coverage = mean(theta_coverage, na.rm = TRUE),
    beta_coverage = mean(beta_coverage, na.rm = TRUE),
    theta_correlation = mean(theta_correlation, na.rm = TRUE),
    beta_correlation = mean(beta_correlation, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(results_frequentist_5b_per_run, "sim_results_frequentist_5b_per_run.csv")
write_csv(results_frequentist_5b_summary, "sim_results_frequentist_5b_summary.csv")
--------------------------------------------------------------------------------
#combine freq csvs----

frequentist_per_run_repeated <- sim_results_frequentist_3_per_run |> 
  bind_rows(sim_results_frequentist_5a_per_run) |> 
  bind_rows(sim_results_frequentist_5b_per_run)

write_csv(frequentist_per_run_repeated, "frequentist_per_run_repeated.csv")

frequentist_summary_repeated <- sim_results_frequentist_3_summary |> 
  bind_rows(sim_results_frequentist_5a_summary) |> 
  bind_rows(sim_results_frequentist_5b_summary)

write_csv(frequentist_summary_repeated, "frequentist_summary_repeated.csv")

--------------------------------------------------------------------------------
#scatterplot data ----
  
n_participants <- 100
n_items <- 50
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(0.8, 1.8)

set.seed(0825)

theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
beta <- rnorm(n_items, mean = mean_beta, sd = sigma_beta)

if (n_items == 600) {
  sim_data <- tibble(
    participant = rep(1:n_participants, each = 100),
    item = unlist(lapply(1:n_participants, function(p) sample(1:n_items, 100)))
  )
} else {
  sim_data <- expand_grid(
    participant = 1:n_participants,
    item = 1:n_items
  )
}

sim_data <- sim_data |> 
  rowwise() |>
  mutate(
    theta = theta[participant],
    beta = beta[item],
    eta = theta + beta,
    p_excl = pnorm(thresholds[1] - eta),
    p_inc = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
    p_id = 1 - pnorm(thresholds[2] - eta),
    y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
  ) |>
  ungroup() |>
  mutate(
    match_y_maxprob = y == max.col(cbind(p_excl, p_inc, p_id), ties.method = "first")
  )

sim_data$y <- factor(sim_data$y,
                     levels = 1:3,
                     labels = c("excl", "inc", "id"),
                     ordered = TRUE)

start_time <- Sys.time()

brm_model <- brm(
  formula = y ~ (1 | participant) + (1 | item),
  data = sim_data,
  family = cumulative(probit),
  chains = 4,
  iter = 2000,
  seed = 0825
)

end_time <- Sys.time()
run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

ranefs <- ranef(brm_model)

theta_estimates <- ranefs$participant[, , "Intercept"] |> 
  as_tibble(rownames = "participant") |> 
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  )

beta_estimates <- ranefs$item[, , "Intercept"] |> 
  as_tibble(rownames = "item") |> 
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )


sim_data_with_estimates <- sim_data |>
  mutate(
    participant = as.character(participant),
    item = as.character(item)
  ) |>
  left_join(theta_estimates, by = "participant") |>
  left_join(beta_estimates, by = "item") |>
  mutate(
    n_participants = n_participants,
    n_items = n_items,
    threshold_1 = thresholds[1],
    threshold_2 = thresholds[2],
    theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
    beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper),
    corr_theta = cor(theta, theta_estimate),
    corr_beta = cor(beta, beta_estimate),
    run_time = run_time
  )


# Save results
filename <- paste0("true_vs_est_scatter.csv")
write_csv(sim_data_with_estimates, filename)