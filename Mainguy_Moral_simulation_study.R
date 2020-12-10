## GLM-derived fish mortality estimates and comparisons:
## a proposed error structure selection approach for catch-curve analyses
## Mainguy & Moral

## This script reproduces the simulation studies in the paper

## loading packages
library(hnp)
library(glmmTMB)
library(tidyverse)
library(gamlss)
library(gridExtra)

## loading helper functions
source("helper_functions.R")

## reading walleye dataset to serve as basis for simulation study
walleye <- read.csv("walleyePP.csv", stringsAsFactors = TRUE, header = TRUE) %>%
  mutate(YEAR = as.factor(YEAR)) %>%
  filter(YEAR == "2012")

## setting up true values
## scenario 1: true model = equidispersion
## scenario 2: true model = overdispersion (linear)
## scenario 3: true model = overdispersion (quadratic)
## scenario 4: true model = underdispersion

true_beta <- round(coef(glm(N ~ AGE, family = poisson, data = walleye)), 1)
true_disp_nb <- 2
true_disp_dp <- .5

n_sim <- 20
beta <- true_beta
age <- walleye$AGE

## run studies for each scenario - takes approx. 10 hours in total
#results1 <- run_study(age = age, beta = beta, disp = 1, n_sim = n_sim, family = "poisson")
#results2 <- run_study(age = age, beta = beta, disp = 2, n_sim = n_sim, family = "nb1")
#results3 <- run_study(age = age, beta = beta, disp = 2, n_sim = n_sim, family = "nb2")
#results4 <- run_study(age = age, beta = beta, disp = .5, n_sim = n_sim, family = "dp")

## compile results
source("get_summaries.R")

aicc_all %>%
  as_tibble %>%
  mutate(model = rownames(aicc_all)) %>%
  pivot_longer(cols = 1:4,
               names_to = "scenario",
               values_to = "aicc") %>%
  ggplot(aes(x = model, y = aicc)) +
  theme_bw() +
  facet_wrap(~ scenario, scales = "free_y") +
  geom_point()

hnp_all %>%
  as_tibble %>%
  mutate(model = rownames(hnp_all)) %>%
  pivot_longer(cols = 1:4,
               names_to = "scenario",
               values_to = "hnp") %>%
  ggplot(aes(x = model, y = hnp)) +
  theme_bw() +
  facet_wrap(~ scenario, scales = "free_y") +
  geom_point()

rmse_slope_all %>%
  as_tibble %>%
  mutate(model = rownames(rmse_slope_all)) %>%
  pivot_longer(cols = 1:4,
               names_to = "scenario",
               values_to = "rmse_slope") %>%
  ggplot(aes(x = model, y = rmse_slope)) +
  theme_bw() +
  facet_wrap(~ scenario, scales = "free_y") +
  geom_point()

se_slope_all %>%
  as_tibble %>%
  mutate(model = rownames(se_slope_all)) %>%
  pivot_longer(cols = 1:4,
               names_to = "scenario",
               values_to = "se_slope") %>%
  ggplot(aes(x = model, y = se_slope)) +
  theme_bw() +
  facet_wrap(~ scenario, scales = "free_y") +
  geom_point()

###################### plots for paper

fish1 <- read.csv("equi_fish.csv", header = TRUE)[,-1]
fish2 <- read.csv("linear_fish.csv", header = TRUE)[,-1]
fish3 <- read.csv("quad_fish.csv", header = TRUE)[,-1]
fish4 <- read.csv("under_fish.csv", header = TRUE)[,-1]

fish1$scenario <- "equidispersion"
fish2$scenario <- "linear overdispersion"
fish3$scenario <- "quadratic overdispersion"
fish4$scenario <- "underdispersion"

fish <- rbind(fish1, fish2, fish3, fish4)

fish_est <- fish %>%
  filter(type == "estimate") %>%
  dplyr::select(- type) %>%
  pivot_longer(cols = 1:6,
               names_to = "model",
               values_to = "slope") %>%
  dplyr::select(model, slope, scenario)

all_slopes <- rbind(slope1, slope2, slope3, slope4)

model_est <- - all_slopes %>%
  t %>%
  as_tibble %>%
  pivot_longer(cols = 1:28,
               names_to = "model",
               values_to = "slope")
model_est$model <- as.factor(model_est$model)
levels(model_est$model) <- rep(c("CMP","GP","NB1","NB2","PN","Poisson","QP"), each = 4)
model_est$scenario <- gl(4, 7, labels = c("equidispersion",
                                          "linear overdispersion",
                                          "quadratic overdispersion",
                                          "underdispersion"))

fish_est$method <- "fishmethods"
model_est$method <- "GLM"
all_est <- rbind(fish_est, model_est)
all_est$model <- factor(all_est$model, levels = c("LR","WLR","Heincke","CR","CRCB","PM",
                                                  "Poisson","QP","NB1","NB2","CMP","GP","PN"))

fish_se <- fish %>%
  filter(type == "standard_error") %>%
  dplyr::select(- type) %>%
  pivot_longer(cols = 1:6,
               names_to = "model",
               values_to = "se") %>%
  dplyr::select(model, se, scenario)

all_ses <- rbind(se1_slo, se2_slo, se3_slo, se4_slo)
rownames(all_ses) <- substr(rownames(all_ses), 1, nchar(rownames(all_ses)) - 4)

model_se <- all_ses %>%
  t %>%
  as_tibble
colnames(model_se) <- c(colnames(model_se)[1:7], paste0(colnames(model_se)[1:7], "1"),
                        paste0(colnames(model_se)[1:7], "2"), paste0(colnames(model_se)[1:7], "3"))

model_se <- model_se %>%
  pivot_longer(cols = 1:28,
               names_to = "model",
               values_to = "se")
model_se$model <- as.factor(model_se$model)
levels(model_se$model) <- rep(c("CMP","GP","NB1","NB2","PN","Poisson","QP"), each = 4)
model_se$scenario <- gl(4, 7, labels = c("equidispersion",
                                         "linear overdispersion",
                                         "quadratic overdispersion",
                                         "underdispersion"))

fish_se$method <- "fishmethods"
model_se$method <- "GLM"
all_se <- rbind(fish_se, model_se)
all_se$model <- factor(all_se$model, levels = c("LR","WLR","Heincke","CR","CRCB","PM",
                                                "Poisson","QP","NB1","NB2","CMP","GP","PN"))

fake_data <- data.frame(scenario = rep(c("Equi-dispersion","Linear over-dispersion","Quadratic over-dispersion","Under-dispersion"), 2),
                        model = "LR",
                        method = "fishmethods",
                        slope = c(.2,.2,-.2,.2,.7,.7,1.8,.7))

all_est$scenario <- as.factor(all_est$scenario)
levels(all_est$scenario) <- c("Equi-dispersion","Linear over-dispersion",
                              "Quadratic over-dispersion","Under-dispersion")
all_se$scenario <- as.factor(all_se$scenario)
levels(all_se$scenario) <- c("Equi-dispersion","Linear over-dispersion",
                             "Quadratic over-dispersion","Under-dispersion")

weighted_se <- data.frame(scenario = c("Equi-dispersion","Linear over-dispersion","Quadratic over-dispersion","Under-dispersion"),
                          wse = c(0.025,0.042,0.103,0.019))
                                 
p1 <- all_est %>%
  ggplot(aes(x = model, y = slope, fill = method)) +
  theme_bw() +
  facet_wrap(~ scenario, scales = "free_y") +
  #geom_hline(yintercept = 0.5, col = "lightgray", lty = 1, lwd = 1) +
  geom_boxplot(outlier.size = .75, outlier.alpha = .5) +
  geom_hline(yintercept = 0.5, lty = "11", lwd = 1) +
  ylab(expression(italic(hat(Z)))) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, color = 1),
        axis.text.y = element_text(color = 1),
        panel.border = element_rect(colour = 1, fill = NA),
        panel.background = element_rect(colour = 1, size = 1),
        strip.background = element_rect(colour = 1, size = 1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face = "bold", size = 10)) +
  scale_fill_grey(start = .7, end = 1) +
  ggtitle("(a)") +
  geom_point(data = fake_data, col = NA)

p2 <- all_se %>%
  ggplot(aes(x = model, y = se, fill = method)) +
  theme_bw() +
  facet_wrap(~ scenario, scales = "free_y") +
  #geom_hline(data = weighted_se,
  #           aes(yintercept = wse),
  #           col = "lightgray",
  #           lty = 1, lwd = 1) +
  geom_boxplot(outlier.size = .75, outlier.alpha = .5) +
  ylab(expression(paste("SE(",italic(hat(Z)),")"))) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, color = 1),
        axis.text.y = element_text(color = 1),
        panel.border = element_rect(colour = 1, fill = NA),
        panel.background = element_rect(colour = 1, size = 1),
        strip.background = element_rect(colour = 1, size = 1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face = "bold", size = 10)) +
  scale_fill_grey(start = .7, end = 1) +
  ggtitle("(b)") +
  geom_hline(data = weighted_se,
             aes(yintercept = wse),
             lty = "11", lwd = 1)

png("figure1.tiff", res = 800, units = "in", w = 10, h = 12)
grid.arrange(p1, p2, ncol = 1)
dev.off()

png("figure1a.tiff", res = 800, units = "in", w = 10, h = 6)
p1 + ggtitle("")
dev.off()

png("figure1b.tiff", res = 800, units = "in", w = 10, h = 6)
p2 + ggtitle("")
dev.off()
