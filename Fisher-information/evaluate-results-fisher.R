# Copyright Darjus Hosszejni, 2021
# See README and LICENSE for further information

library(tidyverse)

folder_top <- "results-fisher"
main_filename <- file.path(folder_top, "summary-table.RDS")
if (!file.exists(main_filename)) {
  message("File not found")
  tab <- list()

  version <- "2"
  folder <- file.path(folder_top, paste0("results-fisher", version))
  for (subfolder in list.dirs(folder, recursive = FALSE)) {
    message(subfolder)
    filepath <- file.path(subfolder, "result.RDS")
    contents <- readRDS(filepath)
    tab <-
      c(tab, list(tibble(Path = filepath,
                         Nu = contents$nu,
                         Y = contents$y,
                         Mean = contents$moments["mean"],
                         Sd = contents$moments["sd"],
                         Low = contents$quantiles["5%"],
                         Median = contents$quantiles["50%"],
                         High = contents$quantiles["95%"],
                         McmcSeed = contents$mcmc_seed)))
  }

  tab <- do.call(bind_rows, tab)
  tab <- tab %>%
    mutate(Mean = as.vector(Mean),
           Sd = as.vector(Sd),
           Low = as.vector(Low),
           Median = as.vector(Median),
           High = as.vector(High),
           SA_Efficiency = -0.5 * (1 / Nu - digamma(Nu / 2)))

  saveRDS(tab, main_filename)
} else {
  message("File found")
  tab <- readRDS(main_filename)
}

# tab %>%
#   select(Y, Nu, Mean) %>%
#   ggplot(aes(x = Nu, y = Y, z = Mean)) +
#   geom_contour_filled() +
#   scale_x_continuous(trans = "log10") +
#   theme_bw()

pdf("efficiency-contour.pdf", width = 7, height = 5)
tab %>%
  select(Y, Nu, Mean, SA_Efficiency) %>%
  mutate(Value = Mean - SA_Efficiency) %>%
  ggplot(aes(x = Nu, y = Y, z = Value)) +
  geom_contour_filled() +
  scale_x_continuous(trans = "log10") +
  ylab(expression(y)) +
  xlab(expression(nu)) +
  theme_bw()
dev.off()

nu <- sort(unique(tab$Nu))[55]  # ~ 4
pdf("efficiency-line-nu4.pdf", width = 7, height = 7)
tab %>%
  select(Y, Nu, Low, Median, High, SA_Efficiency) %>%
  filter(Nu == nu) %>%
  rename("5%" = Low, "50%" = Median, "95%" = High) %>%
  pivot_longer(c("5%", "50%", "95%"), names_to = "Quantile", values_to = "Value") %>%
  mutate(Quantile = factor(Quantile, levels = c("95%", "50%", "5%"))) %>%
  ggplot(aes(x = Y)) +
  geom_line(aes(y = Value, linetype = Quantile)) +
  scale_linetype_manual(values = c("dashed", "solid", "dashed")) +
  scale_size_manual(values = c(0.7, 1, 0.7)) +
  geom_line(aes(y = SA_Efficiency), color = "red") +
  ggtitle("Efficiency of AA (black) and of SA (red)") +
  theme_bw()
dev.off()

pdf("appendix-efficiency-line-y0.pdf", width = 7, height = 4)
tab %>%
  select(Y, Nu, Low, Median, High) %>%
  filter(Y == 0) %>%
  pivot_longer(c(Low, Median, High), names_to = "Quantile", values_to = "Value") %>%
  ggplot(aes(x = Nu, y = Value)) +
  geom_line(aes(linetype = Quantile)) +
  scale_x_continuous(trans = "log10") +
  theme_bw()
dev.off()

pdf("appendix-efficiency-line-y0.pdf", width = 7, height = 4)
tab %>%
  select(Y, Nu, Low, Median, High) %>%
  filter(Y == 9) %>%
  pivot_longer(c(Low, Median, High), names_to = "Quantile", values_to = "Value") %>%
  ggplot(aes(x = Nu, y = Value)) +
  geom_line(aes(linetype = Quantile)) +
  scale_x_continuous(trans = "log10") +
  theme_bw()
dev.off()

# The contour plot is robust
tab %>%
  select(Y, Nu, Low, Median, High, SA_Efficiency) %>%
  mutate(Value = Low - SA_Efficiency) %>%
  ggplot(aes(x = Nu, y = Y, z = Value)) +
  geom_contour_filled() +
  scale_x_continuous(trans = "log10") +
  ggtitle("Difference in the Efficiencies of AA and SA") +
  theme_bw()

# The contour plot is robust
tab %>%
  select(Y, Nu, Low, Median, High, SA_Efficiency) %>%
  mutate(Value = High - SA_Efficiency) %>%
  ggplot(aes(x = Nu, y = Y, z = Value)) +
  geom_contour_filled() +
  scale_x_continuous(trans = "log10") +
  ggtitle("Difference in the Efficiencies of AA and SA") +
  theme_bw()
