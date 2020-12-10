## GLM-derived fish mortality estimates and comparisons:
## a proposed error structure selection approach for catch-curve analyses
## Mainguy & Moral

## This script generates the half-normal plots with a simulated envelope
## for each of the six models explored on the paper (Poisson, Quasi-Poisson,
## Negative binomial type 2, Mean-parameterized Conway-Maxwell-Poisson, Generalized Poisson and
## Poisson-normal) using the following datasets:
## 1. Walleye full dataset
## 2. Walleye 2012 data only
## 3. Walleye 2017 data only
## 4. Arctic charr full dataset
## 5. Arctic charr Tasiujaq only
## 6. Arctic charr Salluit only

## loading packages
library(hnp)
library(glmmTMB)
library(tidyverse)

## reading datasets
walleye <- read.csv("walleyePP.csv", stringsAsFactors = TRUE, header = TRUE) %>%
  mutate(YEAR = as.factor(YEAR))
walleye <- walleye %>%
  mutate(obs = factor(1:nrow(walleye)))

charr <- read.csv("charrPP.csv", stringsAsFactors = TRUE, header = TRUE) %>%
  mutate(SITE = as.factor(SITE))
charr <- charr %>%
  mutate(obs = factor(1:nrow(charr)))

walleye12 <- walleye %>%
  filter(YEAR == "2012")
walleye12 <- walleye12 %>%
  mutate(obs = factor(1:nrow(walleye12)))
walleye17 <- walleye %>%
  filter(YEAR == "2017")
walleye17 <- walleye17 %>%
  mutate(obs = factor(1:nrow(walleye17)))

charr_tas <- charr %>%
  filter(SITE == "TASIUJAQ")
charr_tas <- charr_tas %>%
  mutate(obs = factor(1:nrow(charr_tas)))
charr_sal <- charr %>%
  filter(SITE == "SALLUIT")
charr_sal <- charr_sal %>%
  mutate(obs = factor(1:nrow(charr_sal)))

####################################################################
## Model fitting
####################################################################

## Walleye full models
walleye_pois <- glm(N ~ AGE * YEAR, family = poisson, data = walleye)
walleye_qpois <- glm(N ~ AGE * YEAR, family = quasipoisson, data = walleye)
walleye_gp <- glmmTMB(N ~ AGE * YEAR, family = genpois, data = walleye)
walleye_nb1 <- glmmTMB(N ~ AGE * YEAR, family = nbinom1, data = walleye)
walleye_nb2 <- glm.nb(N ~ AGE * YEAR, data = walleye)
walleye_comp <- glmmTMB(N ~ AGE * YEAR, family = compois, data = walleye)
walleye_pn <- glmmTMB(N ~ AGE * YEAR + (1 | obs), family = poisson, data = walleye)

## Walleye 2012 models
walleye12_pois <- glm(N ~ AGE, family = poisson, data = walleye12)
walleye12_qpois <- glm(N ~ AGE, family = quasipoisson, data = walleye12)
walleye12_gp <- glmmTMB(N ~ AGE, family = genpois, data = walleye12)
walleye12_nb1 <- glmmTMB(N ~ AGE, family = nbinom1, data = walleye12)
walleye12_nb2 <- glm.nb(N ~ AGE, data = walleye12)
walleye12_comp <- glmmTMB(N ~ AGE, family = compois, data = walleye12)
walleye12_pn <- glmmTMB(N ~ AGE + (1 | obs), family = poisson, data = walleye12)

## Walleye 2017 models
walleye17_pois <- glm(N ~ AGE, family = poisson, data = walleye17)
walleye17_qpois <- glm(N ~ AGE, family = quasipoisson, data = walleye17)
walleye17_gp <- glmmTMB(N ~ AGE, family = genpois, data = walleye17)
walleye17_nb1 <- glmmTMB(N ~ AGE, family = nbinom1, data = walleye17)
walleye17_nb2 <- glm.nb(N ~ AGE, data = walleye17)
walleye17_comp <- glmmTMB(N ~ AGE, family = compois, data = walleye17)
walleye17_pn <- glmmTMB(N ~ AGE + (1 | obs), family = poisson, data = walleye17)

## Charr full models
charr_pois <- glm(N ~ AGE * SITE, family = poisson, data = charr)
charr_qpois <- glm(N ~ AGE * SITE, family = quasipoisson, data = charr)
charr_gp <- glmmTMB(N ~ AGE * SITE, family = genpois, data = charr)
#charr_nb1 <- glmmTMB(N ~ AGE * SITE, family = nbinom1, data = charr)
charr_nb2 <- glm.nb(N ~ AGE * SITE, data = charr)
charr_comp <- glmmTMB(N ~ AGE * SITE, family = compois, data = charr)
charr_pn <- glmmTMB(N ~ AGE * SITE + (1 | obs), family = poisson, data = charr)

## Charr Tasiujaq models
charr_tas_pois <- glm(N ~ AGE, family = poisson, data = charr_tas)
charr_tas_qpois <- glm(N ~ AGE, family = quasipoisson, data = charr_tas)
charr_tas_gp <- glmmTMB(N ~ AGE, family = genpois, data = charr_tas)
#charr_tas_nb1 <- glmmTMB(N ~ AGE, family = nbinom1, data = charr_tas)
charr_tas_nb2 <- glm.nb(N ~ AGE, data = charr_tas)
#charr_tas_comp <- glmmTMB(N ~ AGE, family = compois, data = charr_tas)
charr_tas_pn <- glmmTMB(N ~ AGE + (1 | obs), family = poisson, data = charr_tas)

## Charr Salluit models
charr_sal_pois <- glm(N ~ AGE, family = poisson, data = charr_sal)
charr_sal_qpois <- glm(N ~ AGE, family = quasipoisson, data = charr_sal)
charr_sal_gp <- glmmTMB(N ~ AGE, family = genpois, data = charr_sal)
charr_sal_nb1 <- glmmTMB(N ~ AGE, family = nbinom1, data = charr_sal)
charr_sal_nb2 <- glm.nb(N ~ AGE, data = charr_sal)
charr_sal_comp <- glmmTMB(N ~ AGE, family = compois, data = charr_sal)
charr_sal_pn <- glmmTMB(N ~ AGE + (1 | obs), family = poisson, data = charr_sal)

#########################################################################
## Computing and plotting the half-normal plots with a simulated envelope
#########################################################################

## Walleye full models - half-normal plots
dfun <- function(obj) {
  residuals(obj, type = "pearson")
}
sfun <- function(n, obj) {
  simulate(obj)[[1]]
}
ffun_gp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR, family = genpois, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye_gp)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR, family = genpois, data = walleye), silent = TRUE)
  }
  return(fit)
}
ffun_nb1 <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR, family = nbinom1, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye_gp)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR, family = nbinom1, data = walleye), silent = TRUE)
  }
  return(fit)
}
ffun_comp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR, family = compois, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye_comp)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR, family = compois, data = walleye), silent = TRUE)
  }
  return(fit)
}
ffun_pn <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR + (1 | obs), family = poisson, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye_pn)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR + (1 | obs), family = poisson, data = walleye), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
hnp_walleye_pois <- hnp(walleye_pois, plot = FALSE, resid.type = "pearson")
hnp_walleye_qpois <- hnp(walleye_qpois, plot = FALSE, resid.type = "pearson")
hnp_walleye_gp <- hnp(walleye_gp, plot = FALSE, newclass = TRUE,
                      diagfun = dfun, simfun = sfun, fitfun = ffun_gp)
hnp_walleye_nb1 <- hnp(walleye_nb1, plot = FALSE, newclass = TRUE,
                       diagfun = dfun, simfun = sfun, fitfun = ffun_nb1)
hnp_walleye_nb2 <- hnp(walleye_nb2, plot = FALSE, resid.type = "pearson")
hnp_walleye_comp <- hnp(walleye_comp, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_comp)
hnp_walleye_pn <- hnp(walleye_pn, plot = FALSE, newclass = TRUE,
                      diagfun = dfun, simfun = sfun, fitfun = ffun_pn)

hnp_walleye <- data_frame(residuals = c(hnp_walleye_pois$residuals,
                                        hnp_walleye_qpois$residuals,
                                        hnp_walleye_gp$residuals,
                                        hnp_walleye_nb1$residuals,
                                        hnp_walleye_nb2$residuals,
                                        hnp_walleye_comp$residuals,
                                        hnp_walleye_pn$residuals),
                          lower = c(hnp_walleye_pois$lower,
                                    hnp_walleye_qpois$lower,
                                    hnp_walleye_gp$lower,
                                    hnp_walleye_nb1$lower,
                                    hnp_walleye_nb2$lower,
                                    hnp_walleye_comp$lower,
                                    hnp_walleye_pn$lower),
                          median = c(hnp_walleye_pois$median,
                                     hnp_walleye_qpois$median,
                                     hnp_walleye_gp$median,
                                     hnp_walleye_nb1$median,
                                     hnp_walleye_nb2$median,
                                     hnp_walleye_comp$median,
                                     hnp_walleye_pn$median),
                          upper = c(hnp_walleye_pois$upper,
                                    hnp_walleye_qpois$upper,
                                    hnp_walleye_gp$upper,
                                    hnp_walleye_nb1$upper,
                                    hnp_walleye_nb2$upper,
                                    hnp_walleye_comp$upper,
                                    hnp_walleye_pn$upper),
                          x = c(hnp_walleye_pois$x,
                                hnp_walleye_qpois$x,
                                hnp_walleye_gp$x,
                                hnp_walleye_nb1$x,
                                hnp_walleye_nb2$x,
                                hnp_walleye_comp$x,
                                hnp_walleye_pn$x),
                          model = factor(rep(c("Poisson",
                                               "Quasi-Poisson",
                                               "Generalized Poisson",
                                               "Negative binomial type 1",
                                               "Negative binomial type 2",
                                               "Mean-parameterized Conway-Maxwell-Poisson",
                                               "Poisson-normal"), each = 99),
                                         levels = c("Poisson",
                                                    "Quasi-Poisson",
                                                    "Negative binomial type 1",
                                                    "Negative binomial type 2",
                                                    "Mean-parameterized Conway-Maxwell-Poisson",
                                                    "Generalized Poisson",
                                                    "Poisson-normal")))

plot_walleye <- hnp_walleye %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Walleye 2012 models - half-normal plots
dfun <- function(obj) {
  residuals(obj, type = "pearson")
}
sfun <- function(n, obj) {
  simulate(obj)[[1]]
}
ffun_nb1 <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = nbinom1, data = walleye12), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye12_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = nbinom1, data = walleye12), silent = TRUE)
  }
  return(fit)
}
ffun_gp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = genpois, data = walleye12), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye12_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = genpois, data = walleye12), silent = TRUE)
  }
  return(fit)
}
ffun_comp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = compois, data = walleye12), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye12_comp)
    fit <- try(glmmTMB(response2 ~ AGE, family = compois, data = walleye12), silent = TRUE)
  }
  return(fit)
}
ffun_pn <- function(response) {
  fit <- try(glmmTMB(response ~ AGE + (1 | obs), family = poisson, data = walleye12), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye12_pn)
    fit <- try(glmmTMB(response2 ~ AGE + (1 | obs), family = poisson, data = walleye12), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
hnp_walleye12_pois <- hnp(walleye12_pois, plot = FALSE, resid.type = "pearson")
hnp_walleye12_qpois <- hnp(walleye12_qpois, plot = FALSE, resid.type = "pearson")
hnp_walleye12_gp <- hnp(walleye12_gp, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_gp)
hnp_walleye12_nb1 <- hnp(walleye12_nb1, plot = FALSE, newclass = TRUE,
                         diagfun = dfun, simfun = sfun, fitfun = ffun_nb1)
hnp_walleye12_nb2 <- hnp(walleye12_nb2, plot = FALSE, resid.type = "pearson")
hnp_walleye12_comp <- hnp(walleye12_comp, plot = FALSE, newclass = TRUE,
                          diagfun = dfun, simfun = sfun, fitfun = ffun_comp)
hnp_walleye12_pn <- hnp(walleye12_pn, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_pn)

hnp_walleye12 <- data_frame(residuals = c(hnp_walleye12_pois$residuals,
                                          hnp_walleye12_qpois$residuals,
                                          hnp_walleye12_gp$residuals,
                                          hnp_walleye12_nb1$residuals,
                                          hnp_walleye12_nb2$residuals,
                                          hnp_walleye12_comp$residuals,
                                          hnp_walleye12_pn$residuals),
                            lower = c(hnp_walleye12_pois$lower,
                                      hnp_walleye12_qpois$lower,
                                      hnp_walleye12_gp$lower,
                                      hnp_walleye12_nb1$lower,
                                      hnp_walleye12_nb2$lower,
                                      hnp_walleye12_comp$lower,
                                      hnp_walleye12_pn$lower),
                            median = c(hnp_walleye12_pois$median,
                                       hnp_walleye12_qpois$median,
                                       hnp_walleye12_gp$median,
                                       hnp_walleye12_nb1$median,
                                       hnp_walleye12_nb2$median,
                                       hnp_walleye12_comp$median,
                                       hnp_walleye12_pn$median),
                            upper = c(hnp_walleye12_pois$upper,
                                      hnp_walleye12_qpois$upper,
                                      hnp_walleye12_gp$upper,
                                      hnp_walleye12_nb1$upper,
                                      hnp_walleye12_nb2$upper,
                                      hnp_walleye12_comp$upper,
                                      hnp_walleye12_pn$upper),
                            x = c(hnp_walleye12_pois$x,
                                  hnp_walleye12_qpois$x,
                                  hnp_walleye12_gp$x,
                                  hnp_walleye12_nb1$x,
                                  hnp_walleye12_nb2$x,
                                  hnp_walleye12_comp$x,
                                  hnp_walleye12_pn$x),
                            model = factor(rep(c("Poisson",
                                                 "Quasi-Poisson",
                                                 "Generalized Poisson",
                                                 "Negative binomial type 1",
                                                 "Negative binomial type 2",
                                                 "Mean-parameterized Conway-Maxwell-Poisson",
                                                 "Poisson-normal"), each = 44),
                                           levels = c("Poisson",
                                                      "Quasi-Poisson",
                                                      "Negative binomial type 1",
                                                      "Negative binomial type 2",
                                                      "Mean-parameterized Conway-Maxwell-Poisson",
                                                      "Generalized Poisson",
                                                      "Poisson-normal")))

plot_walleye12 <- hnp_walleye12 %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Walleye 2017 models - half-normal plots
dfun <- function(obj) {
  residuals(obj, type = "pearson")
}
sfun <- function(n, obj) {
  simulate(obj)[[1]]
}
ffun_nb1 <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = nbinom1, data = walleye17), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye17_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = nbinom1, data = walleye17), silent = TRUE)
  }
  return(fit)
}
ffun_gp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = genpois, data = walleye17), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye17_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = genpois, data = walleye17), silent = TRUE)
  }
  return(fit)
}
ffun_comp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = compois, data = walleye17), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye17_comp)
    fit <- try(glmmTMB(response2 ~ AGE, family = compois, data = walleye17), silent = TRUE)
  }
  return(fit)
}
ffun_pn <- function(response) {
  fit <- try(glmmTMB(response ~ AGE + (1 | obs), family = poisson, data = walleye17), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, walleye17_pn)
    fit <- try(glmmTMB(response2 ~ AGE + (1 | obs), family = poisson, data = walleye17), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
hnp_walleye17_pois <- hnp(walleye17_pois, plot = FALSE, resid.type = "pearson")
hnp_walleye17_qpois <- hnp(walleye17_qpois, plot = FALSE, resid.type = "pearson")
hnp_walleye17_gp <- hnp(walleye17_gp, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_gp)
hnp_walleye17_nb1 <- hnp(walleye17_nb1, plot = FALSE, newclass = TRUE,
                         diagfun = dfun, simfun = sfun, fitfun = ffun_nb1)
hnp_walleye17_nb2 <- hnp(walleye17_nb2, plot = FALSE, resid.type = "pearson")
hnp_walleye17_comp <- hnp(walleye17_comp, plot = FALSE, newclass = TRUE,
                          diagfun = dfun, simfun = sfun, fitfun = ffun_comp)
hnp_walleye17_pn <- hnp(walleye17_pn, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_pn)

hnp_walleye17 <- data_frame(residuals = c(hnp_walleye17_pois$residuals,
                                          hnp_walleye17_qpois$residuals,
                                          hnp_walleye17_gp$residuals,
                                          hnp_walleye17_nb1$residuals,
                                          hnp_walleye17_nb2$residuals,
                                          hnp_walleye17_comp$residuals,
                                          hnp_walleye17_pn$residuals),
                            lower = c(hnp_walleye17_pois$lower,
                                      hnp_walleye17_qpois$lower,
                                      hnp_walleye17_gp$lower,
                                      hnp_walleye17_nb1$lower,
                                      hnp_walleye17_nb2$lower,
                                      hnp_walleye17_comp$lower,
                                      hnp_walleye17_pn$lower),
                            median = c(hnp_walleye17_pois$median,
                                       hnp_walleye17_qpois$median,
                                       hnp_walleye17_gp$median,
                                       hnp_walleye17_nb1$median,
                                       hnp_walleye17_nb2$median,
                                       hnp_walleye17_comp$median,
                                       hnp_walleye17_pn$median),
                            upper = c(hnp_walleye17_pois$upper,
                                      hnp_walleye17_qpois$upper,
                                      hnp_walleye17_gp$upper,
                                      hnp_walleye17_nb1$upper,
                                      hnp_walleye17_nb2$upper,
                                      hnp_walleye17_comp$upper,
                                      hnp_walleye17_pn$upper),
                            x = c(hnp_walleye17_pois$x,
                                  hnp_walleye17_qpois$x,
                                  hnp_walleye17_gp$x,
                                  hnp_walleye17_nb1$x,
                                  hnp_walleye17_nb2$x,
                                  hnp_walleye17_comp$x,
                                  hnp_walleye17_pn$x),
                            model = factor(rep(c("Poisson",
                                                 "Quasi-Poisson",
                                                 "Generalized Poisson",
                                                 "Negative binomial type 1",
                                                 "Negative binomial type 2",
                                                 "Mean-parameterized Conway-Maxwell-Poisson",
                                                 "Poisson-normal"), each = 55),
                                           levels = c("Poisson",
                                                      "Quasi-Poisson",
                                                      "Negative binomial type 1",
                                                      "Negative binomial type 2",
                                                      "Mean-parameterized Conway-Maxwell-Poisson",
                                                      "Generalized Poisson",
                                                      "Poisson-normal")))

plot_walleye17 <- hnp_walleye17 %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Charr full models - half-normal plots
dfun <- function(obj) {
  residuals(obj, type = "pearson")
}
sfun <- function(n, obj) {
  simulate(obj)[[1]]
}
ffun_nb1 <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * SITE, family = nbinom1, data = charr), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_gp)
    fit <- try(glmmTMB(response2 ~ AGE * SITE, family = nbinom1, data = charr), silent = TRUE)
  }
  return(fit)
}
ffun_gp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * SITE, family = genpois, data = charr), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_gp)
    fit <- try(glmmTMB(response2 ~ AGE * SITE, family = genpois, data = charr), silent = TRUE)
  }
  return(fit)
}
ffun_comp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * SITE, family = compois, data = charr), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_comp)
    fit <- try(glmmTMB(response2 ~ AGE * SITE, family = compois, data = charr), silent = TRUE)
  }
  return(fit)
}
ffun_pn <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * SITE + (1 | obs), family = poisson, data = charr), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_pn)
    fit <- try(glmmTMB(response2 ~ AGE * SITE + (1 | obs), family = poisson, data = charr), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
hnp_charr_pois <- hnp(charr_pois, plot = FALSE, resid.type = "pearson")
hnp_charr_qpois <- hnp(charr_qpois, plot = FALSE, resid.type = "pearson")
hnp_charr_gp <- hnp(charr_gp, plot = FALSE, newclass = TRUE,
                    diagfun = dfun, simfun = sfun, fitfun = ffun_gp)
#hnp_charr_nb1 <- hnp(charr_nb1, plot = FALSE, newclass = TRUE,
#                     diagfun = dfun, simfun = sfun, fitfun = ffun_nb1)
hnp_charr_nb2 <- hnp(charr_nb2, plot = FALSE, resid.type = "pearson")
hnp_charr_comp <- hnp(charr_comp, plot = FALSE, newclass = TRUE,
                      diagfun = dfun, simfun = sfun, fitfun = ffun_comp)
hnp_charr_pn <- hnp(charr_pn, plot = FALSE, newclass = TRUE,
                    diagfun = dfun, simfun = sfun, fitfun = ffun_pn)

hnp_charr <- data_frame(residuals = c(hnp_charr_pois$residuals,
                                      hnp_charr_qpois$residuals,
                                      hnp_charr_gp$residuals,
                                      #hnp_charr_nb1$residuals,
                                      hnp_charr_nb2$residuals,
                                      hnp_charr_comp$residuals,
                                      hnp_charr_pn$residuals),
                        lower = c(hnp_charr_pois$lower,
                                  hnp_charr_qpois$lower,
                                  hnp_charr_gp$lower,
                                  #hnp_charr_nb1$lower,
                                  hnp_charr_nb2$lower,
                                  hnp_charr_comp$lower,
                                  hnp_charr_pn$lower),
                        median = c(hnp_charr_pois$median,
                                   hnp_charr_qpois$median,
                                   hnp_charr_gp$median,
                                   #hnp_charr_nb1$median,
                                   hnp_charr_nb2$median,
                                   hnp_charr_comp$median,
                                   hnp_charr_pn$median),
                        upper = c(hnp_charr_pois$upper,
                                  hnp_charr_qpois$upper,
                                  hnp_charr_gp$upper,
                                  #hnp_charr_nb1$upper,
                                  hnp_charr_nb2$upper,
                                  hnp_charr_comp$upper,
                                  hnp_charr_pn$upper),
                        x = c(hnp_charr_pois$x,
                              hnp_charr_qpois$x,
                              hnp_charr_gp$x,
                              #hnp_charr_nb1$x,
                              hnp_charr_nb2$x,
                              hnp_charr_comp$x,
                              hnp_charr_pn$x),
                        model = factor(rep(c("Poisson",
                                             "Quasi-Poisson",
                                             "Generalized Poisson",
                                             #"Negative binomial type 1",
                                             "Negative binomial type 2",
                                             "Mean-parameterized Conway-Maxwell-Poisson",
                                             "Poisson-normal"), each = 84),
                                       levels = c("Poisson",
                                                  "Quasi-Poisson",
                                                  #"Negative binomial type 1",
                                                  "Negative binomial type 2",
                                                  "Mean-parameterized Conway-Maxwell-Poisson",
                                                  "Generalized Poisson",
                                                  "Poisson-normal")))

plot_charr <- hnp_charr %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Charr Tasiujaq models half-normal plots
dfun <- function(obj) {
  residuals(obj, type = "pearson")
}
sfun <- function(n, obj) {
  simulate(obj)[[1]]
}
ffun_nb1 <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = nbinom1, data = charr_tas), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_tas_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = nbinom1, data = charr_tas), silent = TRUE)
  }
  return(fit)
}
ffun_gp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = genpois, data = charr_tas), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_tas_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = genpois, data = charr_tas), silent = TRUE)
  }
  return(fit)
}
ffun_pn <- function(response) {
  fit <- try(glmmTMB(response ~ AGE + (1 | obs), family = poisson, data = charr_tas), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_tas_pn)
    fit <- try(glmmTMB(response2 ~ AGE + (1 | obs), family = poisson, data = charr_tas), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
hnp_charr_tas_pois <- hnp(charr_tas_pois, plot = FALSE, resid.type = "pearson")
hnp_charr_tas_qpois <- hnp(charr_tas_qpois, plot = FALSE, resid.type = "pearson")
hnp_charr_tas_gp <- hnp(charr_tas_gp, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_gp)
#hnp_charr_tas_nb1 <- hnp(charr_tas_nb1, plot = FALSE, newclass = TRUE,
#                         diagfun = dfun, simfun = sfun, fitfun = ffun_nb1)
hnp_charr_tas_nb2 <- hnp(charr_tas_nb2, plot = FALSE, resid.type = "pearson")
hnp_charr_tas_pn <- hnp(charr_tas_pn, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_pn)

hnp_charr_tas <- data_frame(residuals = c(hnp_charr_tas_pois$residuals,
                                          hnp_charr_tas_qpois$residuals,
                                          hnp_charr_tas_gp$residuals,
                                          #hnp_charr_tas_nb1$residuals,
                                          hnp_charr_tas_nb2$residuals,
                                          hnp_charr_tas_pn$residuals),
                            lower = c(hnp_charr_tas_pois$lower,
                                      hnp_charr_tas_qpois$lower,
                                      hnp_charr_tas_gp$lower,
                                      #hnp_charr_tas_nb1$lower,
                                      hnp_charr_tas_nb2$lower,
                                      hnp_charr_tas_pn$lower),
                            median = c(hnp_charr_tas_pois$median,
                                       hnp_charr_tas_qpois$median,
                                       hnp_charr_tas_gp$median,
                                       #hnp_charr_tas_nb1$median,
                                       hnp_charr_tas_nb2$median,
                                       hnp_charr_tas_pn$median),
                            upper = c(hnp_charr_tas_pois$upper,
                                      hnp_charr_tas_qpois$upper,
                                      hnp_charr_tas_gp$upper,
                                      #hnp_charr_tas_nb1$upper,
                                      hnp_charr_tas_nb2$upper,
                                      hnp_charr_tas_pn$upper),
                            x = c(hnp_charr_tas_pois$x,
                                  hnp_charr_tas_qpois$x,
                                  hnp_charr_tas_gp$x,
                                  #hnp_charr_tas_nb1$x,
                                  hnp_charr_tas_nb2$x,
                                  hnp_charr_tas_pn$x),
                            model = factor(rep(c("Poisson",
                                                 "Quasi-Poisson",
                                                 "Generalized Poisson",
                                                 #"Negative binomial type 1",
                                                 "Negative binomial type 2",
                                                 "Poisson-normal"), each = 43),
                                           levels = c("Poisson",
                                                      "Quasi-Poisson",
                                                      #"Negative binomial type 1",
                                                      "Negative binomial type 2",
                                                      "Generalized Poisson",
                                                      "Poisson-normal")))

plot_charr_tas <- hnp_charr_tas %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Charr Salluit models half-normal plots
dfun <- function(obj) {
  residuals(obj, type = "pearson")
}
sfun <- function(n, obj) {
  simulate(obj)[[1]]
}
ffun_nb1 <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = nbinom1, data = charr_sal), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_sal_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = nbinom1, data = charr_sal), silent = TRUE)
  }
  return(fit)
}
ffun_gp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = genpois, data = charr_sal), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_sal_gp)
    fit <- try(glmmTMB(response2 ~ AGE, family = genpois, data = charr_sal), silent = TRUE)
  }
  return(fit)
}
ffun_comp <- function(response) {
  fit <- try(glmmTMB(response ~ AGE, family = compois, data = charr_sal), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_sal_comp)
    fit <- try(glmmTMB(response2 ~ AGE, family = compois, data = charr_sal), silent = TRUE)
  }
  return(fit)
}
ffun_pn <- function(response) {
  fit <- try(glmmTMB(response ~ AGE + (1 | obs), family = poisson, data = charr_sal), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- sfun(1, charr_sal_pn)
    fit <- try(glmmTMB(response2 ~ AGE + (1 | obs), family = poisson, data = charr_sal), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
hnp_charr_sal_pois <- hnp(charr_sal_pois, plot = FALSE, resid.type = "pearson")
hnp_charr_sal_qpois <- hnp(charr_sal_qpois, plot = FALSE, resid.type = "pearson")
hnp_charr_sal_gp <- hnp(charr_sal_gp, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_gp)
hnp_charr_sal_nb1 <- hnp(charr_sal_nb1, plot = FALSE, newclass = TRUE,
                         diagfun = dfun, simfun = sfun, fitfun = ffun_nb1)
hnp_charr_sal_nb2 <- hnp(charr_sal_nb2, plot = FALSE, resid.type = "pearson")
hnp_charr_sal_comp <- hnp(charr_sal_comp, plot = FALSE, newclass = TRUE,
                          diagfun = dfun, simfun = sfun, fitfun = ffun_comp)
hnp_charr_sal_pn <- hnp(charr_sal_pn, plot = FALSE, newclass = TRUE,
                        diagfun = dfun, simfun = sfun, fitfun = ffun_pn)

hnp_charr_sal <- data_frame(residuals = c(hnp_charr_sal_pois$residuals,
                                          hnp_charr_sal_qpois$residuals,
                                          hnp_charr_sal_gp$residuals,
                                          hnp_charr_sal_nb1$residuals,
                                          hnp_charr_sal_nb2$residuals,
                                          hnp_charr_sal_comp$residuals,
                                          hnp_charr_sal_pn$residuals),
                            lower = c(hnp_charr_sal_pois$lower,
                                      hnp_charr_sal_qpois$lower,
                                      hnp_charr_sal_gp$lower,
                                      hnp_charr_sal_nb1$lower,
                                      hnp_charr_sal_nb2$lower,
                                      hnp_charr_sal_comp$lower,
                                      hnp_charr_sal_pn$lower),
                            median = c(hnp_charr_sal_pois$median,
                                       hnp_charr_sal_qpois$median,
                                       hnp_charr_sal_gp$median,
                                       hnp_charr_sal_nb1$median,
                                       hnp_charr_sal_nb2$median,
                                       hnp_charr_sal_comp$median,
                                       hnp_charr_sal_pn$median),
                            upper = c(hnp_charr_sal_pois$upper,
                                      hnp_charr_sal_qpois$upper,
                                      hnp_charr_sal_gp$upper,
                                      hnp_charr_sal_nb1$upper,
                                      hnp_charr_sal_nb2$upper,
                                      hnp_charr_sal_comp$upper,
                                      hnp_charr_sal_pn$upper),
                            x = c(hnp_charr_sal_pois$x,
                                  hnp_charr_sal_qpois$x,
                                  hnp_charr_sal_gp$x,
                                  hnp_charr_sal_nb1$x,
                                  hnp_charr_sal_nb2$x,
                                  hnp_charr_sal_comp$x,
                                  hnp_charr_sal_pn$x),
                            model = factor(rep(c("Poisson",
                                                 "Quasi-Poisson",
                                                 "Generalized Poisson",
                                                 "Negative binomial type 1",
                                                 "Negative binomial type 2",
                                                 "Mean-parameterized Conway-Maxwell-Poisson",
                                                 "Poisson-normal"), each = 41),
                                           levels = c("Poisson",
                                                      "Quasi-Poisson",
                                                      "Negative binomial type 1",
                                                      "Negative binomial type 2",
                                                      "Mean-parameterized Conway-Maxwell-Poisson",
                                                      "Generalized Poisson",
                                                      "Poisson-normal")))

plot_charr_sal <- hnp_charr_sal %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## saving all plots

png("walleye_full.tiff", res = 800, units = "in", w = 10, h = 8)
print(plot_walleye)
dev.off()

png("walleye_2012.tiff", res = 800, units = "in", w = 10, h = 8)
print(plot_walleye12)
dev.off()

png("walleye_2017.tiff", res = 800, units = "in", w = 10, h = 8)
print(plot_walleye17)
dev.off()

png("charr_full.tiff", res = 800, units = "in", w = 10, h = 8)
print(plot_charr)
dev.off()

png("charr_tas.tiff", res = 800, units = "in", w = 10, h = 6)
print(plot_charr_tas)
dev.off()

png("charr_sal.tiff", res = 800, units = "in", w = 10, h = 8)
print(plot_charr_sal)
dev.off()

figure1_data <- data_frame(residuals = c(hnp_walleye17_nb2$residuals,
                                         hnp_walleye12_qpois$residuals,
                                         hnp_walleye17_pois$residuals,
                                         hnp_charr_sal_gp$residuals),
                           lower = c(hnp_walleye17_nb2$lower,
                                     hnp_walleye12_qpois$lower,
                                     hnp_walleye17_pois$lower,
                                     hnp_charr_sal_gp$lower),
                           median = c(hnp_walleye17_nb2$median,
                                      hnp_walleye12_qpois$median,
                                      hnp_walleye17_pois$median,
                                      hnp_charr_sal_gp$median),
                           upper = c(hnp_walleye17_nb2$upper,
                                     hnp_walleye12_qpois$upper,
                                     hnp_walleye17_pois$upper,
                                     hnp_charr_sal_gp$upper),
                           x = c(hnp_walleye17_nb2$x,
                                 hnp_walleye12_qpois$x,
                                 hnp_walleye17_pois$x,
                                 hnp_charr_sal_gp$x),
                           model = factor(c(rep("Walleye 2017: NB2", 55),
                                            rep("Walleye 2012: QP", 44),
                                            rep("Walleye 2017: Poisson", 55),
                                            rep("Charr Salluit: GP", 41)),
                                          levels = c("Walleye 2017: NB2",
                                                     "Walleye 2012: QP",
                                                     "Walleye 2017: Poisson",
                                                     "Charr Salluit: GP")))

fig1 <- figure1_data %>%
  ggplot(aes(x = x, y = residuals)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free") +
  geom_point(cex = 1, pch = 16, alpha = .75) +
  geom_line(aes(y = median),
            lty = 2, lwd = .2) +
  geom_line(aes(y = upper),
            lty = 1, lwd = .2) +
  geom_line(aes(y = lower),
            lty = 1, lwd = .2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray", alpha = .2) +
  xlab("Half-normal scores") +
  ylab("Pearson residuals") +
  theme(axis.text.x = element_text(color = 1),
        axis.text.y = element_text(color = 1),
        panel.border = element_rect(colour = 1, fill = NA),
        panel.background = element_rect(colour = 1, size = 1),
        strip.background = element_rect(colour = 1, size = 1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face = "bold", size = 10))

png("figure1_hnp.tiff", res = 800, units = "in", w = 8, h = 6)
print(fig1)
dev.off()
