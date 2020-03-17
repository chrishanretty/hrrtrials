library(tidyverse)
set.seed(843)

psw <- readRDS("working/munged_psw.rds")

alpha <- -0.33
beta_age <- seq(-2, 2, length.out = 8)
### Make the fourth group zero
beta_age <- beta_age - beta_age[4]
beta_educ <- c(seq(0, 4, length.out = 4),
               -1, 0)
beta_ab <- 0.2
beta_owns <- 1.2

area_pop <- 100000
psw <- psw %>%
    mutate(count = round(w8 * area_pop)) %>% 
    mutate(linpred = alpha +
               beta_age[as.numeric(as.factor(psw$age))] +
               beta_educ[as.numeric(as.factor(psw$education))] + 
               beta_ab * p.ab +
               beta_owns * p.owns) %>%
    mutate(prob = exp(linpred) / (1 + exp(linpred)),
           successes = rbinom(n(),
                              size = count,
                              prob = prob))

### Create our individual data, sampled from the full data
ind <- psw %>%
    group_by(group, age, education, p.ab, p.owns) %>%
    slice(c(1, 1)) %>%
    mutate(y = c(1, 0),
           weight = c(unique(successes),
                      unique(count) - unique(successes))) %>%
    dplyr::select(group,
                  age, education,
                  p.ab, p.owns,
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
