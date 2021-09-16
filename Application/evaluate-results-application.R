# Copyright Darjus Hosszejni, 2021
# See README and LICENSE for further information

library(coda)
library(posterior)
library(tidyverse)
library(magrittr)

folder_top <- "results-application"
version <- "5"
main_filename <- file.path(folder_top, paste0("summary-table", version, ".RDS"))
if (!file.exists(main_filename)) {
  message("File not found")
  tab <- list()

  parameterizations <- c(a = "AA", s = "SA", sa = "ASIS")
  folder <- file.path(folder_top, paste0("results-app", version))
  parameter_names <- c("nu", "rho", "gamma", "delta", "a[1]", "a[2]", "a[3]",
                       "a[4]")
  variable_names <- c("cpi" = "CPI", "gnp.real" = "Real GNP",
                      "stock.prices" = "Stock Prices",
                      "gnp.capita" = "Real per Capita GNP",
                      "real.wages" = "Real Wages", "unemp" = "Unemployment Rate",
                      "ip" = "Industrial Production",
                      "gnp.nom" = "Nominal GNP", "vel" = "Velocity",
                      "emp" = "Employment", "int.rate" = "Interest Rate",
                      "nom.wages" = "Wages", "gnp.def" = "GNP Deflator",
                      "money.stock" = "Money Stock")
  for (subfolder in list.dirs(folder, recursive = FALSE)) {
    message(subfolder)
    filepath <- file.path(subfolder, "result.RDS")
    contents <- readRDS(filepath)
    posterior <- list(contents$posterior$nu,
                      contents$posterior$rho,
                      contents$posterior$gamma,
                      contents$posterior$delta,
                      contents$posterior$a[1, ],
                      contents$posterior$a[2, ],
                      contents$posterior$a[3, ],
                      contents$posterior$a[4, ])
    names(posterior) <- parameter_names
    tab <-
      c(tab, list(tibble(Path = filepath,
                         ID = contents$setup$ID,
                         Data = contents$setup$data_name,
                         Parameterization =
                           as.vector(parameterizations[contents$setup$strategy]),
                         McmcSeed = contents$mcmc_seed,
                         RNE = as.vector(sapply(posterior, effectiveSize)) / 10000,
                         Mean = as.vector(sapply(posterior, mean)),
                         Sd = as.vector(sapply(posterior, sd)),
                         Median = as.vector(sapply(posterior, median)),
                         Low = as.vector(sapply(posterior, quantile, probs = 0.10)),
                         High = as.vector(sapply(posterior, quantile, probs = 0.90)),
                         Underflow = sum(contents$posterior$underflow),
                         Parameter = parameter_names)))
  }

  tab <- do.call(bind_rows, tab)
  tab <- tab %>%
    mutate(Parameterization = factor(Parameterization, levels = c("AA", "SA", "ASIS")),
           Parameter = factor(Parameter, levels = parameter_names),
           Data = factor(as.vector(variable_names[Data]), levels = sort(variable_names)))

  saveRDS(tab, main_filename)
} else {
  message("File found")
  tab <- readRDS(main_filename)
}

library(kableExtra)
library(knitr)
options(knitr.kable.NA = "-")

tab %>%
  select(Data, Parameter, Parameterization, Mean) %>%
  filter(Parameter %in% c("nu", "delta", "rho")) %>%
  pivot_wider(names_from = c(Parameter, Parameterization), values_from = Mean,
              names_sort = TRUE) %>%
  arrange(Data)

tab %>%
  select(Data, Parameter, Parameterization, Sd) %>%
  filter(Parameter %in% c("nu", "delta", "rho")) %>%
  pivot_wider(names_from = c(Parameter, Parameterization), values_from = Sd,
              names_sort = TRUE) %>%
  arrange(Data)

tab %>%
  select(Data, Parameter, Parameterization, Sd, Mean, Median) %>%
  filter(Parameter == "nu") %>%
  select(-Parameter) %>%
  pivot_longer(c(Sd, Mean, Median), names_to = "Quantity", values_to = "Value") %>%
  pivot_wider(names_from = c(Quantity, Parameterization), values_from = Value,
              names_sort = TRUE) %>%
  arrange(Data)

# Table of estimates
tab %>%
  select(Data, Parameter, Parameterization, Low, High, Median) %>%
  mutate(CredibleInt = paste0("$(", formatC(Low, 3, 4, "fg"), ", ", formatC(High, 3, 4, "fg"), ")$"),
         Median = paste0("$", formatC(Median, 3, 4, "fg"), "$")) %>%
  filter(Parameter == "nu") %>%
  select(Data, Parameterization, CredibleInt, Median) %>%
  pivot_wider(names_from = Parameterization, values_from = c(Median, CredibleInt),
              names_sort = TRUE) %>%
  arrange(Data) %>%
  kable("latex", digits = 1, caption = "TODO", position = "!t",
        escape = FALSE, booktabs = TRUE)

tab %>%
  filter(Underflow > 0) %>%
  select(Parameterization, Underflow)

# Table of efficiencies
tab %>%
  filter(Parameter == "nu") %>%
  select(Data, Parameterization, RNE, Median) %>%
  group_by(Data) %>%
  mutate(Median = mean(Median)) %>%
  ungroup() %>%
  mutate(Median = paste0("$\\sim", formatC(Median, 3, 4, "fg"), "$")) %>%
  mutate(RNE = 100 * RNE) %>%
  unite("Data_Median", Data, Median) %>%
  pivot_wider(id_cols = Data_Median, names_from = Parameterization, values_from = RNE,
              names_sort = TRUE) %>%
  separate(Data_Median, c("Data", "Median"), sep = "_") %>%
  mutate("AA / SA" = AA / SA) %>%
  arrange(Data) %>%
  select(Data, AA, SA, ASIS, "AA / SA", Median) %>%
  kable("latex", digits = 1, caption = "TODO", position = "!t",
        escape = FALSE, booktabs = TRUE) %>%
  writeClipboard()

# Find differences
par1 <- readRDS("results-application/results-app3/18/result.RDS")$posterior
par2 <- readRDS("results-application/results-app3/32/result.RDS")$posterior
par3 <- readRDS("results-application/results-app3/4/result.RDS")$posterior
param <- "nu"
bind_rows(tibble(Par = "SS", Value = head(par1[[param]], 50000L), Step = 1:50000),
          tibble(Par = "ASIS", Value = par2[[param]], Step = 1:50000),
          tibble(Par = "AA", Value = par3[[param]], Step = 1:50000)) %>%
  ggplot(aes(x = Step, y = Value)) +
  geom_line() +
  facet_wrap(~ Par, ncol = 1)

param <- "tau"
bind_rows(tibble(Par = "SS", Value = head(par1[[param]][1, ], 50000L), Step = 1:50000),
          tibble(Par = "ASIS", Value = par2[[param]][1, ], Step = 1:50000),
          tibble(Par = "AA", Value = par3[[param]][1, ], Step = 1:50000)) %>%
  ggplot(aes(x = Step, y = Value)) +
  geom_line() +
  facet_wrap(~ Par, ncol = 1)
