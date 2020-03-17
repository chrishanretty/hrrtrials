library(tidyverse)
library(rstan)

psw <- readRDS("working/munged_psw.rds")
ind <- readRDS("working/ind_data.rds")
agg <- readRDS("working/agg_outcomes.rds")

### Add on summary statistics to the aggregate outcomes
agg <- merge(agg,
             psw %>%
             group_by(group) %>%
             summarize(isFemale = weighted.mean(isFemale, w8),
                       isOwner = weighted.mean(isOwner, w8),
                       p.ab = unique(p.ab),
                       p.over65 = unique(p.over65)),
             all = TRUE,
             sort = FALSE)

stan_data <- list()
stan_data$npsw <- nrow(psw)
stan_data$nobs <- nrow(ind)
stan_data$npredictors <- 2
stan_data$nareas <- length(unique(psw$group))
stan_data$startpos <- psw %>%
    mutate(pos = 1:n()) %>% 
    group_by(group) %>%
    summarize(pos = min(pos)) %>%
    pull(pos)

stan_data$endpos <- psw %>%
    mutate(pos = 1:n()) %>% 
    group_by(group) %>%
    summarize(pos = max(pos)) %>%
    pull(pos)

stan_data$area <- psw %>%
    mutate(area = factor(group)) %>%
    pull(area) %>%
    as.numeric()

stan_data$indy <- ind$y

agg <- psw %>%
    group_by(group) %>%
    summarize(successes = sum(successes),
              count = sum(count))

stan_data$aggy <- agg$successes
stan_data$aggcount <- agg$count

my_form <- formula(~ p.ab + p.over65)

stan_data$contx <- model.matrix(my_form,
                               data = ind)[,-1]

stan_data$pswx <- model.matrix(my_form,
                               data = psw)[,-1]

stan_data$w8 <- psw$w8

stan_data$ind_own <- ind$isOwner
stan_data$ind_gender <- ind$isFemale

stan_data$psw_own <- psw$isOwner
stan_data$psw_gender <- psw$isFemale

stan_data$ncontpredictors <- 2

stan_code <- '
data {
  int<lower=0> nobs;
  int<lower=0> nareas;
  int ncontpredictors;
  int npsw;

  int<lower=1,upper = nareas> area[npsw];
  matrix[nobs, ncontpredictors] contx;
  matrix[npsw, ncontpredictors] pswx;

  vector[nobs] ind_gender;
  vector[nobs] ind_own;
  vector[npsw] psw_gender;
  vector[npsw] psw_own;

  int<lower=0,upper=1> indy[nobs];
  int<lower=0> aggy[nareas];
  int<lower=0> aggcount[nareas];
  int<lower=1, upper = npsw> startpos[nareas];
  int<lower=1, upper = npsw> endpos[nareas];
  vector[npsw] w8;
}
parameters {
  real alpha;
  vector [ncontpredictors] beta;
  real beta_own;
  real beta_gender;
}
transformed parameters {
  vector[npsw] predp;
  vector[nareas] mu;

  // generate the (weighted) probability of turnout for each psw group
  // Here I use element wise multiplication (.*) for the rows in the frame
  predp = inv_logit(pswx * beta + alpha +
     beta_own * psw_own + beta_gender * psw_gender) .* w8;

  // get the sum by group
  for (a in 1:nareas) {
     mu[a] = sum(predp[startpos[a]:endpos[a]]);
  }
}
model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 1);
  beta_own ~ normal(0, 1); 
  beta_gender ~ normal(0, 1);

  indy ~ bernoulli_logit(contx * beta + alpha +
     beta_own * ind_own + beta_gender * ind_gender);
  aggy ~ binomial(aggcount, mu);
}

'

writeLines(stan_code, con = "m1.stan")

the_model <- stan_model(file = "m1.stan")

### This took about two and a half minutes on my laptop
the_samples <- sampling(the_model,
                        data = stan_data,
                        cores = 3,
                        chains = 3,
                        iter = 1000,
                        control = list(max_treedepth = 14),
                        pars = c("alpha", "beta",
                                 "beta_own", "beta_gender"))

### We can amend the data and see whether it matters for the speed of
### processing. It does, but only very marginally (reduction of 12%)
### The high posterior correlation between the gender and intercept
### didn't go away.

gender_mean <- mean(ind$isFemale)
stan_data$ind_gender <- stan_data$ind_gender - gender_mean

own_mean <- mean(ind$isFemale)
stan_data$ind_own <- stan_data$ind_own - own_mean

the_samples2 <- sampling(the_model,
                        data = stan_data,
                        cores = 3,
                        chains = 3,
                        iter = 1000,
                        control = list(max_treedepth = 14),
                        pars = c("alpha", "beta",
                                 "beta_own", "beta_gender"))
