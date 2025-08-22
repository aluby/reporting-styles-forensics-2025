#load packages ----------------------------------------------------------
library(tidyverse)
library(bayesplot)
library(rstanarm)
library(rstan)
library(brms)
library(ggplot2)
library(dplyr)
library(readr)
library(sn)
library(beepr)
library(broom)
library(MASS)
--------------------------------------------------------------------------------
#prior misspecification ----
--------------------------------------------------------------------------------
##No theta in simulation ----
n_participants <- 100
n_items <- 600
sigma_beta <- 1 
mean_beta <- 0
thresholds <- c(0.8, 1.8)

set.seed(0825)

#you shouldn't need to change anything below this for this section
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
    theta = 0,
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

beep(2)

ranefs <- ranef(brm_model)

# Get estimates for participants
participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
  mutate(participant = as.integer(participant)) |>
  select(participant, Estimate, Q2.5, Q97.5) |>
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  ) 

participant_summary <- participant_summary |>
  mutate(
    theta_significant = theta_lower > 0 | theta_upper < 0
  )

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  select(item, Estimate, Q2.5, Q97.5) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
  mutate(
    n_participants = n_participants,
    n_items = n_items,
    threshold_1 = thresholds[1],
    threshold_2 = thresholds[2],
    theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
    beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper),
    corr_theta = cor(theta, theta_estimate),
    corr_beta = cor(beta, beta_estimate),
    run_time = run_time,
    num_significant_thetas = sum(theta_significant)
  )

#save to csv
filename <- paste0("no_theta_sim_", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)
--------------------------------------------------------------------------------
##No theta in model ----
n_participants <- 100
n_items <- 600
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(0.8, 1.8)

set.seed(0825)

#you shouldn't need to change anything below this for this section
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
    p_id = 1 - pnorm(thresholds[2] - eta)
    ,
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
  formula = y ~ (1 | item),
  data = sim_data,
  family = cumulative(probit),
  chains = 4,
  iter = 2000,
  seed = 0825
)

end_time <- Sys.time()
run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

beep(2)

ranefs <- ranef(brm_model)

participant_summary <- tibble(
  participant = 1:n_participants,
  theta_estimate = NA,
  theta_lower = NA,
  theta_upper = NA,
  theta_covered = NA
)

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  select(item, Estimate, Q2.5, Q97.5) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
  mutate(
    n_participants = n_participants,
    n_items = n_items,
    threshold_1 = thresholds[1],
    threshold_2 = thresholds[2],
    theta_covered = NA,
    beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper),
    corr_theta = NA,
    corr_beta = cor(beta, beta_estimate),
    run_time = run_time
  )

#save to csv
filename <- paste0("no_theta_model_", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)

mean(no_theta_model_600_100$theta_covered)
mean(no_theta_model_600_100$beta_covered)
mean(no_theta_model_600_100$corr_theta)
mean(no_theta_model_600_100$corr_beta)
--------------------------------------------------------------------------------
##Skewed theta ----
n_participants <- 100
n_items <- 600
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(.8, 1.8)
alpha_theta <- 5

set.seed(0825)
  
#you shouldn't need to change anything below this for this section
delta_theta <- alpha_theta/sqrt(1+alpha_theta^2)
omega_theta <- sqrt((sigma_theta^2)/(1-(2*delta_theta^2)/pi))
xi_theta <- -(omega_theta * delta_theta * sqrt(2/pi))

theta <- rsn(n_participants, xi = xi_theta, omega = omega_theta, alpha = alpha_theta)
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
    p_id = 1 - pnorm(thresholds[2] - eta)
    ,
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

beep(2)

ranefs <- ranef(brm_model)

# Get estimates for participants
participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
  mutate(participant = as.integer(participant)) |>
  select(participant, Estimate, Q2.5, Q97.5) |>
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  )

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  select(item, Estimate, Q2.5, Q97.5) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
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

#save to csv
filename <- paste0("skewed_theta_", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)
--------------------------------------------------------------------------------
##Skewed beta ----
n_participants <- 100
n_items <- 600
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(.8, 1.8)
alpha_beta <- 5

set.seed(0825)
  
#you shouldn't need to change anything below this for this section
delta_beta <- alpha_beta/sqrt(1+alpha_beta^2)
omega_beta <- sqrt((sigma_beta^2)/(1-(2*delta_beta^2)/pi))
xi_beta <- -(omega_beta * delta_beta * sqrt(2/pi))

theta <- rnorm(n_participants, mean = mean_theta, sd = sigma_theta)
beta <- rsn(n_items, xi = xi_beta, omega = omega_beta, alpha = alpha_beta)

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
    p_id = 1 - pnorm(thresholds[2] - eta)
    ,
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

beep(2)

ranefs <- ranef(brm_model)

# Get estimates for participants
participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
  mutate(participant = as.integer(participant)) |>
  select(participant, Estimate, Q2.5, Q97.5) |>
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  )

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  select(item, Estimate, Q2.5, Q97.5) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
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

#save to csv
filename <- paste0("skewed_beta_", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)
--------------------------------------------------------------------------------
##both skewed ----

n_participants <- 100
n_items <- 600
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(.8, 1.8)
alpha_theta <- 5
alpha_beta <- 5

set.seed(0825)

#you shouldn't need to change anything below this for this section
delta_theta <- alpha_theta/sqrt(1+alpha_theta^2)
omega_theta <- sqrt((sigma_theta^2)/(1-(2*delta_theta^2)/pi))
xi_theta <- -(omega_theta * delta_theta * sqrt(2/pi))
delta_beta <- alpha_beta/sqrt(1+alpha_beta^2)
omega_beta <- sqrt((sigma_beta^2)/(1-(2*delta_beta^2)/pi))
xi_beta <- -(omega_beta * delta_beta * sqrt(2/pi))

theta <- rsn(n_participants, xi = xi_theta, omega = omega_theta, alpha = alpha_theta)
beta <- rsn(n_items, xi = xi_beta, omega = omega_beta, alpha = alpha_beta)

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
    p_id = 1 - pnorm(thresholds[2] - eta)
    ,
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

beep(2)

ranefs <- ranef(brm_model)

# Get estimates for participants
participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
  mutate(participant = as.integer(participant)) |>
  select(participant, Estimate, Q2.5, Q97.5) |>
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  )

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  select(item, Estimate, Q2.5, Q97.5) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
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

#save to csv
filename <- paste0("both_skewed_", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)
--------------------------------------------------------------------------------
#3 vs 5 categories ----
--------------------------------------------------------------------------------
##3 categories ----

n_participants <- 100
n_items <- 600
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(0.8, 1.8)

set.seed(0825)

#you shouldn't need to change anything below this for this section
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
    p_id = 1 - pnorm(thresholds[2] - eta)
    ,
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

beep(2)

ranefs <- ranef(brm_model)

# Get estimates for participants
participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
  mutate(participant = as.integer(participant)) |>
  select(participant, Estimate, Q2.5, Q97.5) |>
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  )

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  select(item, Estimate, Q2.5, Q97.5) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
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

#save to csv
filename <- paste0("3cat_", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)
--------------------------------------------------------------------------------
##5_a categories ----

n_participants <- 100
n_items <- 600
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(0.8, 1.1, 1.5, 1.8)

set.seed(0825)

#you shouldn't need to change anything below this for this section
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
    p_lean_excl = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
    p_inc = pnorm(thresholds[3] - eta) - pnorm(thresholds[2] - eta),
    p_lean_id = pnorm(thresholds[4] - eta) - pnorm(thresholds[3] - eta),
    p_id = 1 - pnorm(thresholds[4] - eta),
    y = which.max(rmultinom(1, 1, c(p_excl, p_lean_excl, p_inc, p_lean_id, p_id)))
  ) |>
  ungroup() |>
  mutate(
    match_y_maxprob = y == max.col(cbind(p_excl, p_lean_excl, p_inc, p_lean_id, p_id), ties.method = "first")
  )

sim_data$y <- factor(sim_data$y,
                     levels = 1:5,
                     labels = c("excl", "lean_excl", "inc", "lean_id", "id"),
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

beep(2)

ranefs <- ranef(brm_model)

# Get estimates for participants
participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
  mutate(participant = as.integer(participant)) |>
  select(participant, Estimate, Q2.5, Q97.5) |>
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  )

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  select(item, Estimate, Q2.5, Q97.5) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
  mutate(
    n_participants = n_participants,
    n_items = n_items,
    threshold_1 = thresholds[1],
    threshold_2 = thresholds[2],
    threshold_3 = thresholds[3],
    threshold_4 = thresholds[4],
    theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
    beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper),
    corr_theta = cor(theta, theta_estimate),
    corr_beta = cor(beta, beta_estimate),
    run_time = run_time
  )

#save to csv
filename <- paste0("5cat_a_", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)
--------------------------------------------------------------------------------
##5_b categories ----

n_participants <- 100
n_items <- 600
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(0.6, 0.9, 1.5, 2.0)

set.seed(0825)

#you shouldn't need to change anything below this for this section
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
    p_lean_excl = pnorm(thresholds[2] - eta) - pnorm(thresholds[1] - eta),
    p_inc = pnorm(thresholds[3] - eta) - pnorm(thresholds[2] - eta),
    p_lean_id = pnorm(thresholds[4] - eta) - pnorm(thresholds[3] - eta),
    p_id = 1 - pnorm(thresholds[4] - eta),
    y = which.max(rmultinom(1, 1, c(p_excl, p_lean_excl, p_inc, p_lean_id, p_id)))
  ) |>
  ungroup() |>
  mutate(
    match_y_maxprob = y == max.col(cbind(p_excl, p_lean_excl, p_inc, p_lean_id, p_id), ties.method = "first")
  )

sim_data$y <- factor(sim_data$y,
                     levels = 1:5,
                     labels = c("excl", "lean_excl", "inc", "lean_id", "id"),
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

beep(2)

ranefs <- ranef(brm_model)

# Get estimates for participants
participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
  mutate(participant = as.integer(participant)) |>
  select(participant, Estimate, Q2.5, Q97.5) |>
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  )

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  select(item, Estimate, Q2.5, Q97.5) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
  mutate(
    n_participants = n_participants,
    n_items = n_items,
    threshold_1 = thresholds[1],
    threshold_2 = thresholds[2],
    threshold_3 = thresholds[3],
    threshold_4 = thresholds[4],
    theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
    beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper),
    corr_theta = cor(theta, theta_estimate),
    corr_beta = cor(beta, beta_estimate),
    run_time = run_time
  )

#save to csv
filename <- paste0("5cat_b_", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)

--------------------------------------------------------------------------------
#Bayes vs non-Bayes ----
--------------------------------------------------------------------------------
##Bayes ----

n_participants <- 100
n_items <- 600
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(0.8, 1.8)

set.seed(0825)

#you shouldn't need to change anything below this for this section
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
    p_id = 1 - pnorm(thresholds[2] - eta)
    ,
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

beep(2)

ranefs <- ranef(brm_model)

# Get estimates for participants
participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
  mutate(participant = as.integer(participant)) |>
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  )

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
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

#save to csv
filename <- paste0("bayes_sim", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)

--------------------------------------------------------------------------------
##non-Bayes -----

n_participants <- 100
n_items <- 600
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(.3, 1.3)

# set seed so that results are reproducible

set.seed(6127128) 

# simulate parameters for participants and items

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

sim_data <- expand_grid(participant = 1:n_participants,
                                   item = 1:n_items) |>
  rowwise() |>
  mutate(
    theta = theta[participant],
    beta = beta[item],
    p_excl = pnorm(theta + beta + thresholds[1]),
    p_inc = pnorm(theta + beta + thresholds[2]) - pnorm(theta + beta + thresholds[1]),
    p_id = 1 - pnorm(theta + beta + thresholds[2])
  ) |>
  mutate(
    y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
  )
view(sim_data)

sim_data$y <- factor(sim_data$y, levels = 1:3, ordered = TRUE)

# frequentest test 
frequentist_test <- polr(y ~ factor(item) + factor(participant), data = sim_data)

# pulling out item and participant estimates from the model

frequentist_effects <- tidy(frequentist_test)

item_effects <- frequentist_effects %>%
  rename(beta_estimate = estimate) %>% 
  rename(item = term)

item_effects <- item_effects %>% 
  mutate(item = as.integer(gsub("[^0-9]", "", as.character(item)))) %>% 
  mutate(item = as.integer(as.character(item))) 

participant_effects <- frequentist_effects %>%
  rename(theta_estimate = estimate) %>% 
  rename(participant = term)

participant_effects <- participant_effects %>% 
  mutate(participant = as.integer(gsub("[^0-9]", "", as.character(participant)))) %>% 
  mutate(participant = as.integer(as.character(participant))) 

sim_data_with_estimates <- sim_data |>
  left_join(participant_effects, by = "participant") |>
  left_join(item_effects, by = "item") |>
  mutate(
    n_participants = n_participants,
    n_items = n_items,
    threshold_1 = thresholds[1],
    threshold_2 = thresholds[2],
    theta_covered = as.integer(theta >= theta_lower & theta <= theta_upper),
    beta_covered = as.integer(beta >= beta_lower & beta <= beta_upper),
    corr_theta = cor(theta, theta_estimate),
    corr_beta = cor(beta, beta_estimate)
  )

#save to csv
filename <- paste0("frequentest_sim", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)
--------------------------------------------------------------------------------
#brms vs raw STAN ----
n_participants <- 100
n_items <- 20
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(0.8, 1.8)

set.seed(0825)

#you shouldn't need to change anything below this for this section
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
    p_id = 1 - pnorm(thresholds[2] - eta)
    ,
    y = which.max(rmultinom(1, 1, c(p_excl, p_inc, p_id)))
  ) |>
  ungroup() |>
  mutate(
    match_y_maxprob = y == max.col(cbind(p_excl, p_inc, p_id), ties.method = "first")
  )

sim_3cat <- sim_data |> 
  mutate(decision3 = y)

sim_3cat$decision3 <- as.integer(sim_3cat$decision3)
sim_3cat$item <- as.factor(sim_3cat$item)
sim_3cat$participant <- as.factor(sim_3cat$participant)

x <- model.matrix(~ item + 0, data = sim_3cat)
new_order <- order(as.numeric(gsub("item", "", colnames(x))))
X <- x[, new_order]

z <- model.matrix(~ participant + 0, data = sim_3cat)
new_order <- order(as.numeric(gsub("participant", "", colnames(z))))
Z <- z[, new_order]

# Data list to pass into stan
grm_probit_stan_datalist_sim_3cat <- list(
  K = nlevels(as.factor(sim_3cat$decision3)),
  P = length(unique(sim_3cat$participant)), # 100 examiners
  D = length(unique(sim_3cat$item)), # 600 image pairs
  N = length(sim_3cat$decision3),
  y = sim_3cat$decision3,
  # item = as.integer(factor(sim_3cat$item)),
  # par = as.integer(factor(sim_3cat$participant)),
  x = X,
  z = Z
  #c = c(1.5, 4.5)
)

# Fit the grm to Busey data
grm_probit_normprior_3cat_WithTheta2_stan_sim_3cat <-
  stan(file = "~/Documents/GitHub/forensic-stats-research/notebooks/vivian-notebook/Raw Stan Code/grm_normprior_3cat_WithTheta2.stan", 
       data = grm_probit_stan_datalist_sim_3cat, 
       chains = 4, 
       iter = 2000,
       seed = 84735)

#########

# Check correlation between our beta estimates and simulated beta 
grm_probit_normprior_stansummary_mu_b_df_sim_3cat <- 
  rstan::summary(grm_probit_normprior_3cat_WithTheta2_stan_sim_3cat, pars = "beta")$summary |>
  as.data.frame() |>
  rownames_to_column(var = "item")

grm_probit_normprior_stansummary_mu_b_df_sim_3cat$item <- as.numeric(gsub(".*\\[|\\]", "", grm_probit_normprior_stansummary_mu_b_df_sim_3cat$item))

grm_probit_normprior_stansummary_mu_b_df_sim_3cat |>
  left_join(sim_3cat, by = "item") |>
  select(mean, `50%`, beta, item) |>
  mutate(label = ifelse(`50%` < quantile(`50%`, 0.10) | `50%` > quantile(`50%`, 0.90),
                        as.character(item), NA)) |>
  ggplot(aes(x = `50%`, y = beta)) +
  geom_point(color = "blue") +
  geom_text_repel(aes(label = label), max.overlaps = 10) +
  # geom_vline(xintercept = c(1.5, 4.5), color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  # scale_x_continuous(breaks = c(pretty(grm_probit_normprior_stansummary_mu_b_df_sim_3cat$`50%`), 1.5, 4.5), labels = c(pretty(grm_probit_normprior_stansummary_mu_b_df_sim_3cat$`50%`), "c1 (1.5)", "c2 (4.5)")) +
  labs(x = "Posterior beta estimates from our grm model with a probit link (w/ theta's)", y = "Simulated (true) beta (w/o theta's)")



beta_summary <- rstan::summary(
  grm_probit_normprior_3cat_WithTheta2_stan_sim_3cat,
  pars = "beta"
)$summary |>
  as.data.frame() |>
  rownames_to_column(var = "item")

# Extract item number from "beta[##]"
beta_summary$item <- as.numeric(gsub(".*\\[|\\]", "", beta_summary$item))

# Merge with simulated betas
beta_with_truth <- beta_summary |>
  left_join(sim_3cat, by = "item")

# Compute coverage (95% CI contains true value)
beta_with_truth <- beta_with_truth |>
  mutate(
    beta_covered = beta >= `2.5%` & beta <= `97.5%`
  )

beta_coverage_rate <- mean(beta_with_truth$beta_covered)
cat("Beta coverage rate:", beta_coverage_rate, "\n")

if ("theta[1]" %in% rownames(rstan::summary(
  grm_probit_normprior_3cat_WithTheta2_stan_sim_3cat
)$summary)) {
  
  theta_summary <- rstan::summary(
    grm_probit_normprior_3cat_WithTheta2_stan_sim_3cat,
    pars = "theta"
  )$summary |>
    as.data.frame() |>
    rownames_to_column(var = "participant")
  
  theta_summary$participant <- as.numeric(gsub(".*\\[|\\]", "", theta_summary$participant))
  
  # Merge with simulated thetas
  theta_with_truth <- theta_summary |>
    left_join(sim_3cat_theta, by = "participant")  # sim_3cat_theta should have columns participant + theta
  
  theta_with_truth <- theta_with_truth |>
    mutate(
      theta_covered = theta >= `2.5%` & theta <= `97.5%`
    )
  
  theta_coverage_rate <- mean(theta_with_truth$theta_covered)
  cat("Theta coverage rate:", theta_coverage_rate, "\n")
}

--------------------------------------------------------------------------------
##brms ----

--------------------------------------------------------------------------------
##raw STAN ----


stancode(y ~ (1 | participant) + (1 | item),
         data = sim_data, family = cumulative(probit))


--------------------------------------------------------------------------------
#logit test ----
n_participants <- 100
n_items <- 50
sigma_theta <- 1
sigma_beta <- 1 
mean_beta <- 0
mean_theta <- 0 
thresholds <- c(0.8, 1.8)

set.seed(0825)

#you shouldn't need to change anything below this for this section
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
    p_id = 1 - pnorm(thresholds[2] - eta)
    ,
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
  family = cumulative(logit),
  chains = 4,
  iter = 2000,
  seed = 0825
)

end_time <- Sys.time()
run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

beep(2)

ranefs <- ranef(brm_model)

# Get estimates for participants
participant_summary <- as_tibble(ranefs$participant[, , "Intercept"], rownames = "participant") |>
  mutate(participant = as.integer(participant)) |>
  select(participant, Estimate, Q2.5, Q97.5) |>
  rename(
    theta_estimate = Estimate,
    theta_lower = Q2.5,
    theta_upper = Q97.5
  )

# get estimates for items
item_summary <- as_tibble(ranefs$item[, , "Intercept"], rownames = "item") |>
  mutate(item = as.integer(item)) |>
  select(item, Estimate, Q2.5, Q97.5) |>
  rename(
    beta_estimate = Estimate,
    beta_lower = Q2.5,
    beta_upper = Q97.5
  )

# Join summaries with sim_data
sim_data_with_estimates <- sim_data |>
  left_join(participant_summary, by = "participant") |>
  left_join(item_summary, by = "item") |>
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

#save to csv
filename <- paste0("logit_test_", n_items, "_", n_participants, ".csv")
write_csv(sim_data_with_estimates, filename)

mean(logit_test_50_100$theta_covered)
mean(logit_test_50_100$beta_covered)
mean(logit_test_50_100$corr_theta)
mean(logit_test_50_100$corr_beta)


library(shinystan)
launch_shinystan(brm_model)
