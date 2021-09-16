# Copyright Darjus Hosszejni, 2021
# See README and LICENSE for further information

# ggplot2
library(tidyverse)

pdf("contour-y-small-top.pdf", height=4, width=7)
# Contour plot of p(nu, tau | y = 0)
y <- 0
# Mesh
tibble(Tau = rep(exp(seq(from = log(0.01), to = log(4), length.out = 101)), times = 101),
       Nu = rep(exp(seq(from = log(1), to = log(30), length.out = 101)), each = 101)) %>%
  mutate(Posterior = dnorm(y, 0, sqrt(Tau), log = TRUE) +
           dgamma(1 / Tau, .5 * Nu, .5 * Nu, log = TRUE) - 2 * log(Tau) +
           dexp(Nu, 0.1, log = TRUE)) %>%
  mutate(Posterior = exp(Posterior - max(Posterior))) %>%
  mutate(Posterior = Posterior / sum(Posterior)) %>%
  ggplot(aes(x = Tau, y = Nu, z = Posterior)) +
  geom_contour_filled() +
  xlab(expression(tau)) +
  ylab(expression(nu)) +
  theme_bw()
dev.off()

pdf("contour-y-small-bottom.pdf", height=4, width=7)
# Contour plot of p(nu, tau | y = 0)
y <- 0
# Mesh
tibble(U = rep(seq(from = 0.01, to = 0.99, length.out = 101), times = 101),
       Nu = rep(exp(seq(from = log(1), to = log(30), length.out = 101)), each = 101)) %>%
  mutate(Posterior = dnorm(y, 0, sqrt(1 / qgamma(1 - U, .5 * Nu, .5 * Nu)), log = TRUE) +
             dexp(Nu, 0.1, log = TRUE)) %>%
  mutate(Posterior = exp(Posterior - max(Posterior))) %>%
  mutate(Posterior = Posterior / sum(Posterior)) %>%
  ggplot(aes(x = U, y = Nu, z = Posterior)) +
  geom_contour_filled() +
  xlab(expression(u)) +
  ylab(expression(nu)) +
  theme_bw()
dev.off()

pdf("contour-y-large-top.pdf", height=4, width=7)
# Contour plot of p(nu, tau | y = 4)
y <- 4
# Mesh
tibble(Tau = rep(exp(seq(from = log(0.01), to = log(20), length.out = 101)), times = 101),
       Nu = rep(exp(seq(from = log(1), to = log(4), length.out = 101)), each = 101)) %>%
  mutate(Posterior = dnorm(y, 0, sqrt(Tau), log = TRUE) +
           dgamma(1 / Tau, .5 * Nu, .5 * Nu, log = TRUE) - 2 * log(Tau) +
           dexp(Nu, 0.1, log = TRUE)) %>%
  mutate(Posterior = exp(Posterior - max(Posterior))) %>%
  mutate(Posterior = Posterior / sum(Posterior)) %>%
  ggplot(aes(x = Tau, y = Nu, z = Posterior)) +
  geom_contour_filled() +
  xlab(expression(tau)) +
  ylab(expression(nu)) +
  theme_bw()
dev.off()

pdf("contour-y-large-bottom.pdf", height=4, width=7)
# Contour plot of p(nu, tau | y = 4)
y <- 4
# Mesh
tibble(U = rep(seq(from = 0.01, to = 0.99, length.out = 101), times = 101),
       Nu = rep(exp(seq(from = log(1), to = log(4), length.out = 101)), each = 101)) %>%
  mutate(Posterior = dnorm(y, 0, sqrt(1 / qgamma(1 - U, .5 * Nu, .5 * Nu)), log = TRUE) +
             dexp(Nu, 0.1, log = TRUE)) %>%
  mutate(Posterior = exp(Posterior - max(Posterior))) %>%
  mutate(Posterior = Posterior / sum(Posterior)) %>%
  ggplot(aes(x = U, y = Nu, z = Posterior)) +
  geom_contour_filled() +
  xlab(expression(u)) +
  ylab(expression(nu)) +
  theme_bw()
dev.off()
