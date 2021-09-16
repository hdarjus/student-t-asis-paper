# Copyright Darjus Hosszejni, 2021
# See README and LICENSE for further information

message("Start")

options(warn = 1)

library(sample.student.asis)
library(coda)  # effectiveSize and geweke.diag
library(posterior)  # rhat

version <- "3"
# task_id <- 18731 for testing asis
# 9938 for testing ancillary
# 1466 for testing sufficient
task_id <- as.integer(Sys.getenv("SGE_TASK_ID", -1))
message("Task ID: ", task_id)

n_repetitions <- seq_len(5L)
n_values <- c(1,3,10,30,100,300,1000,3000,10000)#,30000,100000)
nu_values <- c(1,1.5,2,2.5,3,4,5,10,20,50,100)
lambda_values <- c(0.05, 0.1, 0.2, 0.5, 1)
parameterizations <- c("s", "a", "sa")
data_grid <- expand.grid(repetition = n_repetitions,
                         n = n_values,
                         nu = nu_values,
                         lambda = lambda_values,
                         parameterization = parameterizations,
                         stringsAsFactors = FALSE)
data_grid <- cbind(ID = seq_len(NROW(data_grid)), data_grid)
# shuffle tasks for better load balancing (longer tasks are clustered)
set.seed(1)
data_row <- sample(NROW(data_grid))[task_id]
parameters <- as.list(data_grid[data_row, ])

message("Parameters:")
print(as.data.frame(parameters))

attach(parameters)
if (parameterization != "sa") {
  quit(status = 0, save = "no")
}

# generate data
data_seed <- repetition + floor(100000 * sqrt(pi))
set.seed(data_seed)
data <- simulate_ancillarity(n, nu)$y
message("Head of data:")
print(head(data))

# set "hyper-values"
n_chains <- 4L
n_sim <- 10000L
burnin <- 1000L
store_latent_index <- as.integer(c(1, median(c(1, n)), n))

# get random seed
mcmc_seed <- as.integer(exp(3) * task_id)
set.seed(mcmc_seed)
message("Random seed generated: ", mcmc_seed)

# run sampler
results <-
  lapply(seq_len(n_chains), function (chain) {
           message("Chain ", chain)
           # initilize nu
           init_nu <- c(0.5, 2, 10, 100)[chain]
           runtime <- system.time({
             result <- dfsample(data, draws = n_sim, burnin = burnin, init_nu = init_nu,
                                strategy = parameterization, on_underflow = "reject",
                                store_latent_indices = store_latent_index,
                                lower_bound = 0,
                                target_acceptance = 0.44)
           })
           if (parameterization == "sa") {
             result$nu <- result$nu[seq(2, length(result$nu), by = 2)]
           }

           n_underflow_rejection <- sum(result$underflow)
           #rne_nu <- effectiveSize(result$nu) / length(result$nu)
           #rne_u <- apply(result$u, 1, effectiveSize) / NCOL(result$u)
           #rne_tau <- apply(result$tau, 1, effectiveSize) / NCOL(result$tau)
           acceptance_rate <- (function (x) mean(head(x,-1) != tail(x,-1)))(result$nu)

           geweke_diagnostic <- c(nu = geweke.diag(result$nu)$z,
                                  tryCatch(sapply(apply(result$tau, 1, geweke.diag), function (x) x$z), error = function (e) rep_len(NA, 3)),
                                  tryCatch(sapply(apply(result$u, 1, geweke.diag), function (x) x$z), error = function (e) rep_len(NA, 3)))
           list(runtime = runtime,
                result = result,
                n_underflow_rejection = n_underflow_rejection,
                acceptance_rate = acceptance_rate,
                geweke_diagnostic = geweke_diagnostic)
           }
  )

nu <- rvar(sapply(results, function (x) x$result$nu), with_chains = TRUE)
tau <- array(do.call(cbind, lapply(results, function (x) t(x$result$tau))),
             c(rev(dim(results[[1]]$result$tau)), n_chains))
tau <- aperm(tau, c(1, 3, 2))
tau <- rvar(tau, with_chains = TRUE)
u <- array(do.call(cbind, lapply(results, function (x) t(x$result$u))),
           c(rev(dim(results[[1]]$result$u)), n_chains))
u <- aperm(u, c(1, 3, 2))
u <- rvar(u, with_chains = TRUE)

post <-
  switch(parameterization,
         a = draws_rvars(nu = nu, u = u),
         s = draws_rvars(nu = nu, tau = tau),
         sa = draws_rvars(nu = nu, u = u, tau = tau),
         stop("Unknown parameterization"))

runtime <- sapply(results, function (x) x$runtime)
acceptance_rate <- sapply(results, function (x) x$acceptance_rate)
geweke_diagnostic <- sapply(results, function (x) x$geweke_diagnostic)
underflow_rejection <- sapply(results, function (x) x$n_underflow_rejection)

to_save <- list(data = data,
                posterior = post,
                runtime = runtime,
                setup = parameters,
                data_seed = data_seed,
                mcmc_seed = mcmc_seed,
                ess_bulk = sapply(post, ess_bulk),
                ess_tail = sapply(post, ess_tail),
                rhat = sapply(post, rhat),
                #rne = c(nu = rne_nu, u = rne_u, tau = rne_tau),
                geweke_diagnostic = geweke_diagnostic,
                underflow_rejection = underflow_rejection,
                acceptance_rate = acceptance_rate)

new_folder <- paste0("results-simulation/results", version, "/", task_id)
dir.create(new_folder, recursive = TRUE, showWarnings = FALSE)
setwd(new_folder)
saveRDS(to_save, "result.RDS")

message("Runtime:")
print(runtime)

