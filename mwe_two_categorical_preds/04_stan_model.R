library(tidyverse)
library(rstan)

psw <- readRDS("working/munged_psw.rds")
ind <- readRDS("working/ind_data.rds")
agg <- readRDS("working/agg_outcomes.rds")

### Add on summary statistics to the aggregate outcomes
agg <- merge(agg,
             psw %>%
             group_by(group) %>%
             summarize(p.ab = unique(p.ab),
                       p.owns = unique(p.owns)),
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

stan_data$aggy <- agg$successes
stan_data$aggcount <- agg$count

my_form <- formula(~ p.ab + p.owns)

stan_data$contx <- model.matrix(my_form,
                               data = ind)[,-1]

stan_data$pswx <- model.matrix(my_form,
                               data = psw)[,-1]

stan_data$w8 <- psw$w8

stan_data$ind_age <- as.numeric(factor(ind$age))
stan_data$ind_educ <- as.numeric(factor(ind$education))

stan_data$psw_age <- as.numeric(factor(psw$age))
stan_data$psw_educ <- as.numeric(factor(psw$education))

stan_data$ncontpredictors <- 2
stan_data$age_cats <- length(unique(ind$age))
stan_data$educ_cats <- length(unique(ind$education))

stan_code <- '
data {
  int<lower=0> nobs;
  int<lower=0> nareas;
  int ncontpredictors;
  int npsw;

  int<lower=1,upper = nareas> area[npsw];
  matrix[nobs, ncontpredictors] contx;
  matrix[npsw, ncontpredictors] pswx;

  int<lower=0> age_cats;
  int<lower=0> educ_cats;
  int<lower=1, upper = age_cats> ind_age[nobs];
  int<lower=1, upper = educ_cats> ind_educ[nobs];
  int<lower=1, upper = age_cats> psw_age[npsw];
  int<lower=1, upper = educ_cats> psw_educ[npsw];

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
  vector [age_cats] beta_age;
  vector [educ_cats] beta_educ;
  real<lower=0> age_sigma;
  real<lower=0> educ_sigma;

}
transformed parameters {
  vector[npsw] predp;
  vector[nareas] mu;

  // generate the (weighted) probability of turnout for each psw group
  // Here I use element wise multiplication (.*) for the rows in the frame
  predp = inv_logit(pswx * beta + alpha +
     beta_educ[psw_educ] + beta_age[psw_age]) .* w8;

  // get the sum by group
  for (a in 1:nareas) {
     mu[a] = sum(predp[startpos[a]:endpos[a]]);
  }
}
model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 1);
  beta_educ ~ normal(0, educ_sigma); 
  beta_age ~ normal(0, age_sigma);
  educ_sigma ~ normal(0, 1);
  age_sigma ~ normal(0, 1);

  indy ~ bernoulli_logit(contx * beta + alpha +
     beta_educ[ind_educ] + beta_age[ind_age]);
  aggy ~ binomial(aggcount, mu);
}

'

writeLines(stan_code, con = "m1.stan")

the_model <- stan_model(file = "m1.stan")

### This took way too long
the_samples <- sampling(the_model,
                        data = stan_data,
                        cores = 3,
                        chains = 3,
                        iter = 1000,
                        control = list(max_treedepth = 14),
                        pars = c("alpha", "beta",
                                 "beta_educ", "beta_age"))

