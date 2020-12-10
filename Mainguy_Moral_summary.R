## GLM-derived fish mortality estimates and comparisons:
## a proposed error structure selection approach for catch-curve analyses
## Mainguy & Moral

## This script summarises all results from the following files:
## Mainguy_Moral_full_model_simulations.R
## Mainguy_Moral_split_model_simulations.R

library(tidyverse)

load("Mainguy_Moral_sim_results.RData")

## Walleye and Atlantic charr -- full model
results_full %>%
  group_by(species, model) %>%
  summarise(mean = mean(perc_out),
            sd = sd(perc_out))

results_full %>%
  ggplot(aes(x = model, y = perc_out, col = model)) +
  theme_bw() +
  geom_boxplot() +
  #geom_jitter(aes(col = model), alpha = .4) +
  facet_wrap(~ species) +
  ylab("% points out of envelope") +
  xlab("Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Split between sites and years
results_split %>%
  group_by(species, split, model) %>%
  summarise(mean = mean(perc_out),
            sd = sd(perc_out))

results_split %>%
  ggplot(aes(x = model, y = perc_out, col = model)) +
  theme_bw() +
  geom_boxplot() +
  #geom_jitter(aes(col = model), alpha = .4) +
  facet_wrap(~ species + split) +
  ylab("% points out of envelope") +
  xlab("Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))