# Copyright Darjus Hosszejni, 2021
# See README and LICENSE for further information

message("Start")

options(warn = 1)

library("tseries")
library("sample.student.asis")

data("NelPlo")

version <- "5"
task_id <- as.integer(Sys.getenv("SGE_TASK_ID", -1))
message("Task ID: ", task_id)
data_grid <- expand.grid(data_name = colnames(NelPlo),
                         strategy = c("a", "s", "sa"),
                         stringsAsFactors = FALSE)
data_grid <- cbind(ID = seq_len(NROW(data_grid)), data_grid)
data_row <- task_id
parameters <- as.list(data_grid[data_row, ])

message("Parameters:")
print(as.data.frame(parameters))

attach(parameters)

mu_delta <- 0
sigma2_delta <- 0.05^2
sigma2_gamma <- 10^2
lambda_nu <- 0.333
mu_a <- rep_len(0, 4)
sigma2_a <- 0.731 * 0.342^(1:4)

draw_multivariate_normal <- function (weighted_observations, weighted_covariates) {
  post_mean <- solve(crossprod(weighted_covariates),
                     crossprod(weighted_covariates, weighted_observations))
  post_stdev <- t(chol(solve(crossprod(weighted_covariates))))
  drop(post_stdev %*% rnorm(NCOL(weighted_covariates)) + post_mean)
}

draw_delta_gamma <- function (y, rho, a, sigma2, tau) {
  n <- length(y) - 5L
  mu_gamma <- y[5L]
  w <- numeric(n + 2L)
  z1 <- numeric(length(w))
  z2 <- numeric(length(w))
  weights <- numeric(length(w))

  # priors
  w[n + 1L] <- mu_gamma
  z1[n + 1L] <- 1
  z2[n + 1L] <- 0
  weights[n + 1L] <- sigma2_gamma
  w[n + 2L] <- mu_delta
  z1[n + 2L] <- 0
  z2[n + 2L] <- 1
  weights[n + 2L] <- sigma2_delta

  # likelihood
  w[1:n] <- y[5 + 1:n] - (rho + a[1]) * y[4 + 1:n] + (a[1] - a[2]) * y[3 + 1:n] +
    (a[2] - a[3]) * y[2 + 1:n] + (a[3] - a[4]) * y[1 + 1:n] + a[4] * y[1:n]
  z1[1:n] <- 1 - rho
  z2[1:n] <- rho - sum(a) + (1 - rho) * (1:n)
  weights[1:n] <- sigma2 * tau

  # draw (gamma, delta)
  sqrt_precision <- 1 / sqrt(weights)
  weighted_covariates <- cbind(z1, z2) * sqrt_precision
  draw_multivariate_normal(weighted_observations = w * sqrt_precision,
                           weighted_covariates = weighted_covariates)
}

draw_a <- function (y, gamma, rho, delta, sigma2, tau) {
  n <- length(y) - 5L
  w <- numeric(n + 4L)
  z <- matrix(nrow = n + 4L, ncol = 4L)
  weights <- numeric(n + 4L)

  # priors
  w[n + 1:4] <- 0
  z[n + 1:4, ] <- diag(x = rep_len(1, 4))
  weights[n + 1:4] <- sigma2_a

  # likelihood
  w[1:n] <- y[5 + 1:n] - gamma * (1 - rho) - delta * rho -
    delta * (1 - rho) * (1:n) - rho * y[4 + 1:n]
  z[1:n, ] <- cbind(y[4 + 1:n] - y[3 + 1:n],
                    y[3 + 1:n] - y[2 + 1:n],
                    y[2 + 1:n] - y[1 + 1:n],
                    y[1 + 1:n] - y[0 + 1:n]) - delta
  weights[1:n] <- sigma2 * tau

  sqrt_precision <- 1 / sqrt(weights)
  weighted_covariates <- z * sqrt_precision
  draw_multivariate_normal(weighted_observations = w * sqrt_precision,
                           weighted_covariates = weighted_covariates)
}

draw_sigma2 <- function (y, gamma, rho, delta, a, tau) {
  n <- length(y) - 5L
  epsilon <- y[5 + 1:n] - (rho + a[1]) * y[4 + 1:n] + (a[1] - a[2]) * y[3 + 1:n] +
    (a[2] - a[3]) * y[2 + 1:n] + (a[3] - a[4]) * y[1 + 1:n] + a[4] * y[1:n] -
    gamma * (1 - rho) - delta * (1 - rho) * (1:n) -
    delta * (rho - sum(a))
  coefficient <- sum(epsilon^2 / tau)
  coefficient / rchisq(1, n)
}

draw_rho <- function (y, gamma, delta, a, sigma2, tau) {
  n <- length(y) - 5L
  s <- 4

  w <- y[5 + 1:n] - a[1] * y[4 + 1:n] + (a[1] - a[2]) * y[3 + 1:n] +
    (a[2] - a[3]) * y[2 + 1:n] + (a[3] - a[4]) * y[1 + 1:n] + a[4] * y[1:n] -
    gamma - delta * (1:n) + delta * sum(a)
  z <- y[4 + 1:n] - gamma + delta - delta * (1:n)

  lambda2 <- sigma2 / sum(z^2 / tau)
  rho_hat <- sum(w * z / tau) * lambda2 / sigma2

  # notation from page 26-27 of Geweke (1992): sigma2 := lambda2 and mu := rho_hat
  # https://www.minneapolisfed.org/research/discussion-papers/priors-for-macroeconomic-time-series-and-their-application
  v <- 0.5 * (rho_hat + sqrt(rho_hat^2 + 4 * lambda2 * s))
  v <- if (s * lambda2 <= v - rho_hat) {
    v
  } else {
    1
  }
  acceptance_rate_fix_part <- v^(-s) * exp(0.5 * (v - rho_hat)^2 / lambda2)
  accepted <- FALSE
  while (!accepted) {
    x <- -1
    while (x < 0 || x >= 1) {
      x <- rnorm(1, v, sqrt(lambda2))
    }
    acceptance_rate <- acceptance_rate_fix_part *
      exp(0.5 * ((x - v)^2 - (x - rho_hat)^2) / lambda2) * x^s
    accepted <- runif(1) <= acceptance_rate
  }
  x
}

sample_us_macro <- function (y, draws = 10000L, burnin = 1000L,
                             nu_strategy = "sa") {
  n <- length(y) - 5L

  # initialize tau and nu
  tau <- rep_len(1, n)
  nu <- 1 / lambda_nu

  # initialize gamma and delta using least squares for y_t = gamma + delta * t + u_t
  time_trend <- seq_along(y)
  init_gd <- lm(y ~ time_trend)$coefficients
  gamma <- c(init_gd["(Intercept)"])
  delta <- c(init_gd["time_trend"])

  # initialize rho given gamma and delta
  init_r <- ar(y - gamma + delta * (1 - time_trend), order.max = 1L, aic = FALSE, demean = FALSE)
  rho <- c(init_r$ar)

  # initialize 'a' given everything else above
  w <- y[5 + 1:n] - gamma * (1 - rho) - delta * rho -
    delta * (1 - rho) * (1:n) - rho * y[4 + 1:n]
  z <- cbind(y[4 + 1:n] - y[3 + 1:n],
             y[3 + 1:n] - y[2 + 1:n],
             y[2 + 1:n] - y[1 + 1:n],
             y[1 + 1:n] - y[0 + 1:n]) - delta
  init_a <- lm(w ~ 0 + z)
  a <- c(init_a$coefficients)

  # initialize sigma2
  sigma2 <- var(w - as.vector(z %*% a))

  # storage
  store_tau_indices <- ceiling(seq(1, n, length.out = 6))
  store_rho <- numeric(draws)
  store_a <- matrix(nrow = 4, ncol = draws)
  store_gamma <- numeric(draws)
  store_delta <- numeric(draws)
  store_tau <- matrix(nrow = length(store_tau_indices), ncol = draws)
  store_nu <- numeric(draws)
  store_sigma2 <- numeric(draws)
  store_underflow <- numeric(draws)

  # samplers
  parameterization <- strsplit(nu_strategy, "")[[1]]
  do_ancillarity <- "a" %in% parameterization
  do_sufficiency <- "s" %in% parameterization
  target_acceptance <- 0.44
  draw_nu_ancillary <-
    construct_nu_ancillarity_sampler(prior_rate = lambda_nu,
                                   lower_bound = 0,
                                   on_underflow = "reject",
                                   target_acceptance = target_acceptance)
  #draw_nu_sufficient <-
  #  construct_nu_sufficiency_sampler(prior_rate = lambda_nu,
  #                                   lower_bound = 0,
  #                                   target_acceptance = target_acceptance)
  compute_student_t_y <- function (y, gamma, delta, rho, a, sigma2) {
    n <- length(y) - 5L
    w <- y[5 + 1:n] - gamma * (1 - rho) - delta * rho -
      delta * (1 - rho) * (1:n) - rho * y[4 + 1:n]
    z <- cbind(y[4 + 1:n] - y[3 + 1:n],
               y[3 + 1:n] - y[2 + 1:n],
               y[2 + 1:n] - y[1 + 1:n],
               y[1 + 1:n] - y[0 + 1:n]) - delta
    (w - as.vector(z %*% a)) / sqrt(sigma2)
  }

  # MCMC loop
  numwidth <- ceiling(log10(burnin + draws + 1))
  start_time <- Sys.time()
  for (m in seq_len(burnin + draws)) {
    if (m %% 200L == 0L) {
      now <- Sys.time()
      message("\r",
              "Passes: ", formatC(m, width = numwidth), " / ", burnin + draws,
              ", Elapsed: ", formatC(difftime(now, start_time, units = "secs"), digits = 1, format = "f", width = 5), "s",
              ", Remaining: ", formatC((burnin + draws - m) / m * difftime(now, start_time, units = "secs"), digits = 1, format = "f", width = 5), "s",
              appendLF = FALSE)
    }

    # sample parameters
    gd <- draw_delta_gamma(y = y, rho = rho, a = a, sigma2 = sigma2, tau = tau)
    gamma <- gd[1]
    delta <- gd[2]
    a <- draw_a(y = y, gamma = gamma, rho = rho, delta = delta, sigma2 = sigma2, tau = tau)
    rho <- draw_rho(y = y, gamma = gamma, delta = delta, a = a, sigma2 = sigma2, tau = tau)
    student_t_y <- compute_student_t_y(y = y, gamma = gamma, delta = delta, a = a, sigma2 = sigma2, rho = rho)
    tau <- update_tau_vector(y = student_t_y, nu = nu)
    if (do_sufficiency) {
      nu <- update_nu_sufficiency(tau = tau, prior_rate = lambda_nu, lower_bound = 0)
    }
    if (do_ancillarity) {
      u <- pinvgamma(tau = tau, nu = nu)
      nu <- draw_nu_ancillary(nu = nu, y = student_t_y, u = u)
      if (!is.finite(nu[1])) {
        nu <- nu[2]
        store_underflow[m] <- TRUE
      }
      tau <- qinvgamma(u = u, nu = nu)
    }
    sigma2 <- draw_sigma2(y = y, gamma = gamma, rho = rho, delta = delta, a = a, tau = tau)

    # store draws
    if (m > burnin) {
      index <- m - burnin
      store_rho[index] <- rho
      store_a[, index] <- a
      store_gamma[index] <- gamma
      store_delta[index] <- delta
      store_tau[, index] <- tau[store_tau_indices]
      store_nu[index] <- nu
      store_sigma2[index] <- sigma2
    }
  }
  message("\r                                                                                          ", appendLF = FALSE)
  message("\rElapsed: ", round(difftime(Sys.time(), start_time, units = "secs"), 1), "s")

  list(rho = store_rho,
       a = store_a,
       gamma = store_gamma,
       delta = store_delta,
       tau = store_tau,
       nu = store_nu,
       sigma2 = store_sigma2,
       underflow = store_underflow,
       y = y)
}

# get random seed
mcmc_seed <- 999476L + as.integer(exp(4) * task_id)
set.seed(mcmc_seed)
message("Random seed generated: ", mcmc_seed)

input <- NelPlo[, data_name]

message("Start sampler")
result <- switch(strategy,
                 s = sample_us_macro(na.omit(input), 10000, 1000, nu_strategy = "s"),
                 a = sample_us_macro(na.omit(input), 10000, 1000, nu_strategy = "a"),
                 sa = sample_us_macro(na.omit(input), 10000, 1000, nu_strategy = "sa"),
                 stop("Unknown strategy"))

to_save <- list(posterior = result,
                setup = parameters,
                mcmc_seed = mcmc_seed,
                data = input)

new_folder <- paste0("results-application/results-app", version, "/", task_id)
dir.create(new_folder, recursive = TRUE, showWarnings = FALSE)
setwd(new_folder)
saveRDS(to_save, "result.RDS")

coda::effectiveSize(result$nu)
coda::effectiveSize(result$sigma2)
quantile(result$nu, probs = c(1, 5, 10, 25, 50, 75, 90, 95, 99) / 100)
quantile(result$sigma2, probs = c(1, 5, 10, 25, 50, 75, 90, 95, 99) / 100)

