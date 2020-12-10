## GLM-derived fish mortality estimates and comparisons:
## a proposed error structure selection approach for catch-curve analyses
## Mainguy & Moral

## This script generates the summary statistics to reproduce Table 3 of the paper

## loading packages
library(hnp)
library(glmmTMB)
library(tidyverse)

## reading datasets
walleye <- read.csv("walleyePP.csv", stringsAsFactors = TRUE, header = TRUE) %>%
  mutate(YEAR = as.factor(YEAR))
charr <- read.csv("charrPP.csv", stringsAsFactors = TRUE, header = TRUE) %>%
  mutate(SITE = as.factor(SITE))

## Poisson fits
walleye_pois <- glm(N ~ AGE * YEAR, family = poisson, data = walleye)
charr_pois <- glm(N ~ AGE * SITE, family = poisson, data = charr)

## Poisson hnp
set.seed(2020)
walleye_hnp_pois <- charr_hnp_pois <- list()
for(i in 1:100) {
  walleye_hnp_pois[[i]] <- hnp(walleye_pois,
                               resid.type = "pearson",
                               how.many.out = TRUE,
                               plot.sim = FALSE)
  charr_hnp_pois[[i]] <- hnp(charr_pois,
                             resid.type = "pearson",
                             how.many.out = TRUE,
                             plot.sim = FALSE)
}

## Poisson hnp summary
walleye_pois_hnp_summary <- sapply(walleye_hnp_pois, function(x) x$out/x$total * 100)
charr_pois_hnp_summary <- sapply(charr_hnp_pois, function(x) x$out/x$total * 100)

## Quasi-Poisson fits
walleye_qpois <- glm(N ~ AGE * YEAR, family = quasipoisson, data = walleye)
charr_qpois <- glm(N ~ AGE * SITE, family = quasipoisson, data = charr)

## Quasi-Poisson hnp
set.seed(2020)
walleye_hnp_qpois <- charr_hnp_qpois <- list()
for(i in 1:100) {
  walleye_hnp_qpois[[i]] <- hnp(walleye_qpois,
                                resid.type = "pearson",
                                how.many.out = TRUE,
                                plot.sim = FALSE)
  charr_hnp_qpois[[i]] <- hnp(charr_qpois,
                              resid.type = "pearson",
                              how.many.out = TRUE,
                              plot.sim = FALSE)
}

## Quasi-Poisson hnp summary
walleye_qpois_hnp_summary <- sapply(walleye_hnp_qpois, function(x) x$out/x$total * 100)
charr_qpois_hnp_summary <- sapply(charr_hnp_qpois, function(x) x$out/x$total * 100)

## Generalized Poisson fit
walleye_gp <- glmmTMB(N ~ AGE * YEAR, family = genpois, data = walleye)
charr_gp <- glmmTMB(N ~ AGE * SITE, family = genpois, data = charr)

## Generalized Poisson hnp
d_fun <- function(obj) {
  residuals(obj, type = "pearson")
}
s_fun <- function(n, obj) {
  simulate(obj)[[1]] 
}
f_fun_w <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR, family = genpois, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- s_fun(1, walleye_gp)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR, family = genpois, data = walleye), silent = TRUE)
  }
  return(fit)
}
f_fun_c <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * SITE, family = genpois, data = charr), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- s_fun(1, charr_gp)
    fit <- try(glmmTMB(response2 ~ AGE * SITE, family = genpois, data = charr), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
walleye_hnp_gp <- charr_hnp_gp <- list()
for(i in 1:100) {
  walleye_hnp_gp[[i]] <- hnp(walleye_gp, newclass = TRUE,
                              diagfun = d_fun, simfun = s_fun, fitfun = f_fun_w,
                              how.many.out = TRUE,
                              plot.sim = FALSE)
  charr_hnp_gp[[i]] <- hnp(charr_gp, newclass = TRUE,
                            diagfun = d_fun, simfun = s_fun, fitfun = f_fun_c,
                            how.many.out = TRUE,
                            plot.sim = FALSE)
}

## Generalized Poisson hnp summary
walleye_gp_hnp_summary <- sapply(walleye_hnp_gp, function(x) x$out/x$total * 100)
charr_gp_hnp_summary <- sapply(charr_hnp_gp, function(x) x$out/x$total * 100)

## Negative binomial type-1 fit
walleye_nb1 <- glmmTMB(N ~ AGE * YEAR, family = nbinom1, data = walleye)
#charr_nb1 <- glmmTMB(N ~ AGE * SITE, family = nbinom1, data = charr)

## Negative binomial type-1 hnp
d_fun <- function(obj) {
  residuals(obj, type = "pearson")
}
s_fun <- function(n, obj) {
  simulate(obj)[[1]] 
}
f_fun_w <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR, family = nbinom1, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- s_fun(1, walleye_nb1)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR, family = nbinom1, data = walleye), silent = TRUE)
  }
  return(fit)
}
f_fun_c <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * SITE, family = nbinom1, data = charr), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- s_fun(1, charr_nb1)
    fit <- try(glmmTMB(response2 ~ AGE * SITE, family = nbinom1, data = charr), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
walleye_hnp_nb1 <- charr_hnp_nb1 <- list()
for(i in 1:100) {
  walleye_hnp_nb1[[i]] <- hnp(walleye_nb1, newclass = TRUE,
                               diagfun = d_fun, simfun = s_fun, fitfun = f_fun_w,
                               how.many.out = TRUE,
                               plot.sim = FALSE)
#  charr_hnp_nb1[[i]] <- hnp(charr_nb1, newclass = TRUE,
#                             diagfun = d_fun, simfun = s_fun, fitfun = f_fun_c,
#                             how.many.out = TRUE,
#                             plot.sim = FALSE)
}

## Negative binomial type-1 hnp summary
walleye_nb1_hnp_summary <- sapply(walleye_hnp_nb1, function(x) x$out/x$total * 100)
#charr_nb1_hnp_summary <- sapply(charr_hnp_nb1, function(x) x$out/x$total * 100)
charr_nb1_hnp_summary <- rep(NA, 100)

## Negative binomial type-2 fits
walleye_nb2 <- glm.nb(N ~ AGE * YEAR, data = walleye)
charr_nb2 <- glm.nb(N ~ AGE * SITE, data = charr)

## Negative binomial type-2 hnp
set.seed(2020)
walleye_hnp_nb2 <- charr_hnp_nb2 <- list()
for(i in 1:100) {
  walleye_hnp_nb2[[i]] <- hnp(walleye_nb2,
                             resid.type = "pearson",
                             how.many.out = TRUE,
                             plot.sim = FALSE)
  charr_hnp_nb2[[i]] <- hnp(charr_nb2,
                           resid.type = "pearson",
                           how.many.out = TRUE,
                           plot.sim = FALSE)
}

## Negative binomial type-2 hnp summary
walleye_nb2_hnp_summary <- sapply(walleye_hnp_nb2, function(x) x$out/x$total * 100)
charr_nb2_hnp_summary <- sapply(charr_hnp_nb2, function(x) x$out/x$total * 100)

## COM-Poisson fit
walleye_comp <- glmmTMB(N ~ AGE * YEAR, family = compois, data = walleye)
charr_comp <- glmmTMB(N ~ AGE * SITE, family = compois, data = charr)

## COM-Poisson hnp
d_fun <- function(obj) {
  residuals(obj, type = "pearson")
}
s_fun <- function(n, obj) {
  simulate(obj)[[1]] 
}
f_fun_w <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR, family = compois, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- s_fun(1, walleye_comp)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR, family = compois, data = walleye), silent = TRUE)
  }
  return(fit)
}
f_fun_c <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * SITE, family = compois, data = charr), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- s_fun(1, charr_comp)
    fit <- try(glmmTMB(response2 ~ AGE * SITE, family = compois, data = charr), silent = TRUE)
  }
  return(fit)
}

set.seed(2020)
walleye_hnp_comp <- charr_hnp_comp <- list()
for(i in 1:100) {
  walleye_hnp_comp[[i]] <- hnp(walleye_comp, newclass = TRUE,
                               diagfun = d_fun, simfun = s_fun, fitfun = f_fun_w,
                               how.many.out = TRUE,
                               plot.sim = FALSE)
  charr_hnp_comp[[i]] <- hnp(charr_comp, newclass = TRUE,
                             diagfun = d_fun, simfun = s_fun, fitfun = f_fun_c,
                             how.many.out = TRUE,
                             plot.sim = FALSE)
}

## COM-Poisson hnp summary
walleye_comp_hnp_summary <- sapply(walleye_hnp_comp, function(x) x$out/x$total * 100)
charr_comp_hnp_summary <- sapply(charr_hnp_comp, function(x) x$out/x$total * 100)

## Poisson-normal fits
walleye$obs <- factor(1:nrow(walleye))
charr$obs <- factor(1:nrow(charr))
walleye_pn <- glmmTMB(N ~ AGE * YEAR + (1 | obs), family = poisson, data = walleye)
charr_pn <- glmmTMB(N ~ AGE * SITE + (1 | obs), family = poisson, data = charr)

f_fun_pn_w <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * YEAR + (1 | obs), family = poisson, data = walleye), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- s_fun(1, walleye_comp)
    fit <- try(glmmTMB(response2 ~ AGE * YEAR + (1 | obs), family = poisson, data = walleye), silent = TRUE)
  }
  return(fit)
}
f_fun_pn_c <- function(response) {
  fit <- try(glmmTMB(response ~ AGE * SITE + (1 | obs), family = poisson, data = charr), silent = TRUE)
  while(class(fit) == "try-error") {
    response2 <- s_fun(1, charr_comp)
    fit <- try(glmmTMB(response2 ~ AGE * SITE + (1 | obs), family = poisson, data = charr), silent = TRUE)
  }
  return(fit)
}

## Poisson-normal hnp
set.seed(2020)
walleye_hnp_pn <- charr_hnp_pn <- list()
for(i in 1:100) {
  walleye_hnp_pn[[i]] <- hnp(walleye_pn,
                             newclass = TRUE,
                             diagfun = d_fun, simfun = s_fun, fitfun = f_fun_pn_w,
                             how.many.out = TRUE,
                             plot.sim = FALSE)
  charr_hnp_pn[[i]] <- hnp(charr_pn,
                           newclass = TRUE,
                           diagfun = d_fun, simfun = s_fun, fitfun = f_fun_pn_c,
                           how.many.out = TRUE,
                           plot.sim = FALSE)
}

## Poisson hnp summary
walleye_pn_hnp_summary <- sapply(walleye_hnp_pn, function(x) x$out/x$total * 100)
charr_pn_hnp_summary <- sapply(charr_hnp_pn, function(x) x$out/x$total * 100)

## All summaries
results <- data.frame(model = factor(rep(c("Poisson",
                                           "Quasi-Poisson",
                                           "Generalized Poisson",
                                           "Negative binomial type-1",
                                           "Negative binomial type-2",
                                           "COM Poisson",
                                           "Poisson-normal"), each = 100),
                                     levels = c("Poisson",
                                                "Quasi-Poisson",
                                                "Generalized Poisson",
                                                "Negative binomial type-1",
                                                "Negative binomial type-2",
                                                "COM Poisson",
                                                "Poisson-normal")),
                      species = factor(rep(c("Walleye","Arctic charr"), each = 700)),
                      perc_out = c(walleye_pois_hnp_summary,
                                   walleye_qpois_hnp_summary,
                                   walleye_gp_hnp_summary,
                                   walleye_nb1_hnp_summary,
                                   walleye_nb2_hnp_summary,
                                   walleye_comp_hnp_summary,
                                   walleye_pn_hnp_summary,
                                   charr_pois_hnp_summary,
                                   charr_qpois_hnp_summary,
                                   charr_gp_hnp_summary,
                                   charr_nb1_hnp_summary,
                                   charr_nb2_hnp_summary,
                                   charr_comp_hnp_summary,
                                   charr_pn_hnp_summary))

results_full <- results