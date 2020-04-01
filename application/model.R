### Load libraries

library(rstan)
library(tidyverse)

### Load the different sets of data
### individual data
### aggregate data
### post-strat frame
load("stan_data_workspace.RData")

stan_code <- '
data {
  int<lower=0> nobs;
  int<lower=0> nareas;
  int ncontpredictors;
  int npsw;
  int nregions;

  int<lower=1,upper = nareas> area[npsw];
  matrix[nobs, ncontpredictors] contx;
  matrix[npsw, ncontpredictors] pswx;

  int<lower=0> age_cats;
  int<lower=0> sex_cats;
  int<lower=0> ethn_cats;

  int<lower=1, upper = age_cats> ind_age[nobs];
  int<lower=1, upper = nregions> ind_region[nobs];
  int<lower=1, upper = nareas> ind_const[nobs];

  vector[nobs] ind_sex;
  vector[nobs] ind_ethn;

  int<lower=1, upper = age_cats> psw_age[npsw];
  int<lower=1, upper = nregions> psw_region[npsw];

  vector[npsw] psw_sex;
  vector[npsw] psw_ethn;


  int<lower=0,upper=1> indy[nobs];
  int<lower=0> aggy[nareas];
  int<lower=0> aggcount[nareas];
  int<lower=1, upper = npsw> startpos[nareas];
  int<lower=1, upper = npsw> endpos[nareas];
  vector<lower=0, upper = 1>[npsw] w8;
}
parameters {
  real alpha;
  vector [ncontpredictors] beta;
  vector [age_cats] beta_age;
  vector [nregions] beta_region;
  vector [nareas] beta_const;
  vector<lower=0, upper=1>[nareas] theta;
  real beta_sex;
  real beta_ethn;
  real<lower=0> age_sigma;
  real<lower=0> region_sigma;
  real<lower=0> const_sigma;
  real<lower=0> inv_kappa;
}
transformed parameters {
  vector[npsw] predp;
  vector<lower=0, upper=1>[nareas] mu;
  real<lower=0> kappa;
  vector [age_cats] beta_age_tr;
  vector [nregions] beta_region_tr;
  vector [nareas] beta_const_tr;

  beta_region_tr = region_sigma * beta_region;
  beta_age_tr = age_sigma * beta_age;
  beta_const_tr = const_sigma * beta_const;
  kappa = 1 / inv_kappa;

  // generate the (weighted) probability of turnout for each psw group
  // Here I use element wise multiplication (.*) for the rows in the frame
  predp = inv_logit(pswx * beta + alpha +
     beta_age_tr[psw_age] +
     beta_region_tr[psw_region] +
     beta_const_tr[area] +
     beta_ethn * psw_ethn +
     beta_sex * psw_sex) .* w8;

  // get the sum by group
  for (a in 1:nareas) {
     mu[a] = sum(predp[startpos[a]:endpos[a]]);
  }
}
model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 1);
  inv_kappa ~ normal(0, 0.015);  // kappa btwn 33 and Inf

  beta_age ~ normal(0, 1);
  beta_region ~ normal(0, 1);
  beta_const ~ normal(0, 1);
  beta_ethn ~ normal(0, 1);
  beta_sex ~ normal(0, 1);
  age_sigma ~ normal(0, 1);
  region_sigma ~ normal(0, 1);

  indy ~ bernoulli_logit(contx * beta + alpha +
     beta_age_tr[ind_age] +
     beta_region_tr[ind_region] +
     beta_const_tr[ind_const] +
     beta_sex * (ind_sex) +
     beta_ethn * (ind_ethn));

  theta ~ beta_proportion(mu,
      kappa);
  aggy ~ binomial(aggcount, theta);
}

'

writeLines(stan_code, con = "m4.stan")

the_model <- stan_model(file = "m4.stan")

stan_data$ind_const <- as.numeric(ind$PCON16CD)
stan_data$psw_const <- as.numeric(psw$PCON16CD)

## Generate inits
ind$age <- factor(ind$age)
ind$age <- relevel(ind$age, "35 to 44")
ind$y <- ind$voted == "voted"
init_mod <- glm(y ~ turn05 + Maj05 + age + sex + ethnicity,
                data = ind,
                family = binomial)

init_coefs  <- coef(init_mod)
my_inits <- list()
my_inits$alpha <- init_coefs[["(Intercept)"]]
my_inits$beta <- init_coefs[c("turn05", "Maj05")]
my_inits$beta_sex <- init_coefs[grep("^sex", names(init_coefs))]
my_inits$beta_ethn <- init_coefs[grep("^ethn", names(init_coefs))]
my_inits$beta_age <- c(0, init_coefs[grep("^age", names(init_coefs))])
my_inits$beta_const <- rep(0, 200)
my_inits$beta_region <- rep(0, 11)

my_inits$inv_kappa <- 1/33

my_inits <- list(my_inits,
     lapply(my_inits, function(x) x * 2),
     lapply(my_inits, function(x) x * 0.5))

the_samples <- sampling(the_model,
                        data = stan_data,
                        cores = 3,
                        chains = 3,
                        iter = 2000,
                        init = my_inits,
                        seed = 158,
                        control = list(max_treedepth = 14,
                                           adapt_delta = 0.95),
                        pars = c("alpha", "beta", "kappa",
                                 "beta_sex", 
                                 "age_sigma",
                                 "region_sigma",
                                 "const_sigma",
                                 "beta_ethn",
                                 "beta_const_tr",
                                 "mu",
                                 "theta",
                                 "beta_age_tr", "beta_region_tr"))


saveRDS(the_samples, file = "mcmc_samples.rds")

