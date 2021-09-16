# Copyright Darjus Hosszejni, 2021
# See README and LICENSE for further information

library(coda)
library(posterior)
library(tidyverse)
library(magrittr)

folder_top <- "results-simulation"
version <- "3"
main_filename <- file.path(folder_top, paste0("summary-table", version, ".RDS"))
if (!file.exists(main_filename)) {
  message("File not found")
  tab <- list()

  parameterizations <- c(a = "AA", s = "SA", sa = "ASIS")
  folder <- file.path(folder_top, paste0("results", version))
  for (subfolder in list.dirs(folder, recursive = FALSE)) {
    message(subfolder)
    filepath <- file.path(subfolder, "result.RDS")
    contents <- readRDS(filepath)
    chains <- cbind(draws_of(contents$posterior$nu)[1:10000L],
                    draws_of(contents$posterior$nu)[10001L:20000L],
                    draws_of(contents$posterior$nu)[20001L:30000L],
                    draws_of(contents$posterior$nu)[30001L:40000L])
    rne <- as.vector(apply(chains, 2, effectiveSize)) / 10000
    quantiles <- apply(chains, 2, quantile, probs = c(10, 50, 90) / 100)
    quantile_joint <- quantile(chains, probs = c(10, 50, 90) / 100)
    moments <- apply(chains, 2, function (x) c(mean = mean(x), sd = sd(x)))
    moments_joint <- c(mean = mean(chains), sd = sd(chains))
    tab <-
      c(tab, list(tibble(Path = filepath,
                         ID = contents$setup$ID,
                         Length = contents$setup$n,
                         Nu = contents$setup$nu,
                         DataID = contents$setup$repetition,
                         Lambda = contents$setup$lambda,
                         Parameterization =
                           as.vector(parameterizations[contents$setup$parameterization]),
                         DataSeed = contents$data_seed,
                         McmcSeed = contents$mcmc_seed,
                         Q10 = as.vector(quantiles["10%", ]),
                         Q50 = as.vector(quantiles["50%", ]),
                         Q90 = as.vector(quantiles["90%", ]),
                         "Joint Q10" = as.vector(quantile_joint["10%"]),
                         "Joint Q50" = as.vector(quantile_joint["50%"]),
                         "Joint Q90" = as.vector(quantile_joint["90%"]),
                         Mean = as.vector(moments["mean", ]),
                         Sd = as.vector(moments["sd", ]),
                         "Joint Mean" = as.vector(moments_joint["mean"]),
                         "Joint Sd" = as.vector(moments_joint["sd"]),
                         Chain = 1:4,
                         RNE = rne,
                         Runtime = as.vector(contents$runtime["user.self", ]),
                         Rhat = as.vector(contents$rhat$nu),
                         GewekeDiagnostic = as.vector(contents$geweke_diagnostic[1, ]),
                         UnderflowRejection = contents$underflow_rejection)))
  }

  tab <- do.call(bind_rows, tab)
  tab <- tab %>%
    mutate(Parameterization = factor(Parameterization, levels = c("AA", "SA", "ASIS")))

  saveRDS(tab, main_filename)
} else {
  message("File found")
  tab <- readRDS(main_filename)
}

# 1. Show examples
n <- 1000L
nus <- c(1, 100)
lambda <- 0.2

pdf("example1.pdf", width = 7, height = 4)
tab %>%
  filter(Length == n, Lambda == lambda, Nu == nus[1],
         Chain == 2L, DataID == 1L) %>%
  select(Path, Parameterization) %>%
  mutate(Draws = map(Path, function (p) tibble(Step = 5001L:10000L, Value = draws_of(readRDS(p)$posterior$nu)[15001L:20000L]))) %>%
  unnest(Draws) %>%
  ggplot(aes(x = Step, y = Value)) +
  geom_line() +
  facet_wrap(Parameterization ~ ., nrow = 3, scales = "free_y") +
  ylab(expression(nu)) +
  #ggtitle(expression(paste("Example Traceplot: ", nu==1, ", ", n==1000, ", ", lambda==0.2))) +
  theme_bw()
dev.off()

pdf("example2.pdf", width = 7, height = 4)
tab %>%
  filter(Length == n, Lambda == lambda, Nu == nus[2],
         Chain == 2L, DataID == 1L) %>%
  select(Path, Parameterization) %>%
  mutate(Draws = map(Path, function (p) tibble(Step = 5001L:10000L, Value = draws_of(readRDS(p)$posterior$nu)[15001L:20000L]))) %>%
  unnest(Draws) %>%
  ggplot(aes(x = Step, y = Value)) +
  geom_line() +
  facet_wrap(Parameterization ~ ., nrow = 3) +
  ylab(expression(nu)) +
  #ggtitle(expression(paste("Example Traceplot: ", nu==100, ", ", n==1000, ", ", lambda==0.2))) +
  theme_bw()
dev.off()

# 2. Show consistency
n <- 100L
pdf("consistency.pdf", width = 7, height = 5)
tab %>%
  filter(Length == n, Nu %in% c(1, 2, 3, 5, 10), DataID < 5, Rhat < 1.1) %>%
  mutate(DataID = as.factor(as.character(DataID)),
         RNE = 100 * RNE) %>%
  ggplot(aes(x = Parameterization, y = RNE)) +
  geom_jitter(aes(color = Lambda, shape = DataID)) +
  facet_wrap(~ Nu, nrow = 1, labeller = label_bquote(nu==.(Nu))) +
  xlab("Algorithm") +
  ylab("RNE (%)") +
  #ggtitle(expression(paste("RNE of AA, SA, and ASIS for ", n==100))) +
  theme_bw()
dev.off()

# 3. Show length-dependence
nus <- c(1.5, 10, 100)
pdf("length-dependence.pdf", width = 7, height = 4)
tab %>%
  filter(Nu %in% nus, Rhat < 1.1) %>%
  mutate(RNE = 100 * RNE) %>%
  group_by(Parameterization, Length, Nu) %>%
  summarize("90%" = quantile(RNE, probs = 0.9),
            "50%" = median(RNE),
            "10%" = quantile(RNE, probs = 0.1)) %>%
  ungroup() %>%
  pivot_longer(c("10%", "50%", "90%"), names_to = "Quantile", values_to = "RNE") %>%
  mutate(Quantile = factor(Quantile, levels = c("90%", "50%", "10%"))) %>%
  ggplot(aes(x = Length, color = Parameterization)) +
  geom_line(aes(y = RNE, linetype = Quantile, size = Quantile)) +
  facet_wrap(~ Nu, ncol = 1, labeller = label_bquote(nu==.(Nu))) +
  xlab(expression(n)) +
  ylab("RNE (%)") +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  scale_linetype_manual(values = c("dashed", "solid", "dashed")) +
  scale_size_manual(values = c(0.7, 1, 0.7)) +
  #ggtitle("RNE in Dependence on the Length of the Data") +
  theme_bw()
dev.off()

# 4. Table of values
library(kableExtra)
library(knitr)
options(knitr.kable.NA = "-")

ns <- c(10, 100, 1000, 10000)
tab %>%
  filter(Rhat < 1.1, Length %in% ns) %>%
  group_by(Parameterization, Length, Nu) %>%
  summarize(RNE = 100 * mean(RNE)) %>%
  ungroup() %>%
  pivot_wider(names_from = Nu, values_from = RNE) %>%
  arrange(Parameterization, Length) %>%
  rename(Algorithm = Parameterization, n = Length) %>%
  kable("latex", digits = 1, caption = "TODO", position = "!t",
        booktabs = TRUE) %>%
  collapse_rows(columns = 1) %>%
  writeClipboard()

# 5. Parameter estimates
tab %>%
  filter(DataID == 1, Lambda == 0.2, Length == 1000, Nu %in% c(1, 2, 5, 20, 100)) %>%
  select(DataID, Parameterization, Q10, Q90, "Joint Q10", "Joint Q90", Nu, Chain, Rhat) %>%
  mutate(CredibleInt = paste0("$(", formatC(Q10, 3, 4, "fg"), ", ", formatC(Q90, 3, 4, "fg"), ")",
                              ifelse(is.na(Rhat) | Rhat >= 1.1, "^\\ast$", "$"))) %>%
  select(Chain, Nu, Parameterization, CredibleInt) %>%
  mutate(Initialization = c(0.5, 2, 10, 100)[Chain]) %>%
  select(-Chain) %>%
  pivot_wider(names_from = Initialization, values_from = CredibleInt,
              names_sort = TRUE) %>%
  arrange(Nu, Parameterization) %>%
  rename("True Nu" = Nu,
         "Initial $\\nu=0.5$" = "0.5") %>%
  kable("latex", caption = "Estimated 10th and 90th percentiles for $\\nu$.
	      The same data of length $n=1000$ is used inside each group of $\\nu_\\text{true}$;
        further, $\\lambda=0.2$ for each entry.
        Sets of chains with $\\hat R\\ge1.1$ are marked.",
        position = "!t",
        booktabs = TRUE, escape = FALSE) %>%
  collapse_rows(columns = 1) %>%
  writeClipboard()

tab %>%
  filter(DataID == 1, Lambda == 0.2, Length == 1000, Nu %in% c(1, 2, 5, 20, 100)) %>%
  select(Chain, Nu, Parameterization, Rhat) %>%
  mutate(Initialization = c(0.5, 2, 10, 100)[Chain]) %>%
  select(-Chain) %>%
  pivot_wider(names_from = Initialization, values_from = Rhat,
              names_sort = TRUE) %>%
  arrange(Nu, Parameterization) %>%
  rename("True Nu" = Nu) %>%
  View
