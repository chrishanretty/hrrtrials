library(tidyverse)
set.seed(843)

psw <- readRDS("working/munged_psw.rds")

alpha <- -0.33
beta_female <- -0.2
beta_owner <- 0.9
beta_ab <- 0.2
beta_over65 <- 1.2

area_pop <- 100000
psw <- psw %>%
    mutate(count = round(w8 * area_pop)) %>% 
    mutate(linpred = alpha +
               beta_female * isFemale +
               beta_owner * isOwner + 
               beta_ab * p.ab +
               beta_over65 * p.over65) %>%
    mutate(prob = exp(linpred) / (1 + exp(linpred)),
           successes = rbinom(n(),
                              size = count,
                              prob = prob))

### Create our individual data, sampled from the full data
ind <- psw %>%
    group_by(group, isFemale, isOwner, p.ab, p.over65) %>%
    slice(c(1, 1)) %>%
    mutate(y = c(1, 0),
           weight = c(unique(successes),
                      unique(count) - unique(successes))) %>%
    dplyr::select(group,
                  isFemale, isOwner, 
                  p.ab, p.over65,
                  y, weight) %>%
    ungroup() %>% 
    sample_n(1000)

### Note our aggregate results
agg <- psw %>%
    group_by(group) %>%
    summarize(successes = sum(successes),
              count = sum(count))

saveRDS(ind, file = "working/ind_data.rds")
saveRDS(agg, file = "working/agg_outcomes.rds")
