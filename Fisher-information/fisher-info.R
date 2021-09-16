# Copyright Darjus Hosszejni, 2021
# See README and LICENSE for further information

message("Start")

options(warn = 1)

version <- "2"
task_id <- as.integer(Sys.getenv("SGE_TASK_ID", -1))
message("Task ID: ", task_id)
data_grid <- expand.grid(y = seq(0, 10, length.out = 101L),
                         nu = exp(seq(log(0.1), log(100), length.out = 100L)))
data_grid <- cbind(ID = seq_len(NROW(data_grid)), data_grid)
data_row <- task_id
parameters <- as.list(data_grid[data_row, ])

message("Parameters:")
print(as.data.frame(parameters))

attach(parameters)

# get random seed
rng <- file("/dev/urandom","rb", raw = TRUE) # open connection
mcmc_seed <- readBin(rng,what="integer",n=1) # read some 8-byte integers
close(rng) # close the connection
mcmc_seed <- abs(as.integer(Sys.time()) %% 1000003L + mcmc_seed %% 1000003L) %% 1000003L
message("Random seed generated: ", mcmc_seed)

library(sample.student.asis)
library(numDeriv)

dd_logdnu_aa <- function (nu, u, y) {
  hessian(function (x) log_likelihood_ancillarity(y, u, x), nu,
          method = "Richardson",  # citation in help(grad)
          method.args = list(r = 6))  # recommended by help
}

evaluate_efficiency_aa <- function (nu, y, len) {
  u <- update_u_vector_gibbs(rep_len(y, len), nu)  # iid sampling
  store_values <- sapply(u, function (x, nu, y) -dd_logdnu_aa(nu, x, y), nu = nu, y = y)
  list(quantiles = quantile(store_values, probs = (0:100) / 100),
       moments = c(sd = sd(store_values), mean = mean(store_values)))
}

message("Start sampler")
result <- evaluate_efficiency_aa(nu, y, 10000L)

message("Save results and clean up")
to_save <- list(y = y,
                nu = nu,
                quantiles = result$quantiles,
                moments = result$moments,
                mcmc_seed = mcmc_seed)

new_folder <- paste0("results-fisher/results-fisher", version, "/", task_id)
dir.create(new_folder, recursive = TRUE, showWarnings = FALSE)
setwd(new_folder)
saveRDS(to_save, "result.RDS")

