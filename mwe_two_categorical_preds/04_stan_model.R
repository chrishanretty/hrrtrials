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
  vector[age_cats] beta_age_starred;
  vector[educ_cats] beta_educ_starred;

  // Set first element to zero to match ecoreg
  beta_educ_starred = beta_educ - beta_educ[1];
  // Set fourth element to zero to match ecoreg
  beta_age_starred = beta_age - beta_age[4];

  // generate the (weighted) probability of turnout for each psw group
  // Here I use element wise multiplication (.*) for the rows in the frame
  predp = inv_logit(pswx * beta + alpha +
     beta_educ_starred[psw_educ] + beta_age_starred[psw_age]) .* w8;

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
     beta_educ_starred[ind_educ] + beta_age_starred[ind_age]);
  aggy ~ binomial(aggcount, mu);
}

'

writeLines(stan_code, con = "m1.stan")

the_model <- stan_model(file = "m1.stan")

## Generate inits
init_mod <- glm(y ~ p.ab + p.owns + age + education,
                data = ind,
                family = binomial)
init_coefs  <- coef(init_mod)
my_inits <- list()
my_inits$alpha <- init_coefs[["(Intercept)"]]
my_inits$beta <- init_coefs[c("p.owns", "p.ab")]
my_inits$beta_educ <- c(0,
                        init_coefs[grep("^educ", names(init_coefs))])
my_inits$beta_age <- c(0,
                        init_coefs[grep("^age", names(init_coefs))])

### Estimate the model
mle_samples <- optimizing(the_model,
                          data = stan_data,
                          hessian = TRUE,
                          draws = 100,
                          init = my_inits,
                          as_vector = FALSE)

### Did we achieve convergence?
(mle_samples$return_code == 0)

### Get the coefficients out from the draw
coef_names <- c("alpha", "beta[1]", "beta[2]",
                paste0("beta_age_starred[", 1:8, "]"),
                paste0("beta_educ_starred[", 1:6, "]"))
                
### How do they look?
coef_hat <- apply(mle_samples$theta_tilde[,coef_names], 2, median)
coef_sd <- apply(mle_samples$theta_tilde[,coef_names], 2, sd)

coef.df <- data.frame(var = coef_names,
                      hat = coef_hat,
                      sd = coef_sd,
                      row.names = NULL)

### Compare to the known values
### 
beta_age <- seq(-2, 2, length.out = 8)
### Make the fourth group zero
beta_age <- beta_age - beta_age[4]
beta_educ <- c(seq(0, 4, length.out = 4),
               -1, 0)

coef.df$known <- c(-0.33,
                0.2, 1.2,
                beta_age,
                beta_educ)

### Plot this
ggplot(coef.df, aes(x = var,
                    y = hat,
                    ymin = hat - 1.96 * sd,
                    ymax = hat + 1.96 * sd)) +
    geom_pointrange() +
    geom_point(aes(x = var, y = known),
               colour = "red",
               size = 4,
               alpha = 0.2) + 
    coord_flip() +
    scale_x_discrete("variable") +
    scale_y_continuous("Estimate")

### The coefficients are similar enough; the differences are probably
### due to the regularization which comes from modelling the age and
### education terms as random intercepts

stop()
### This took way too long
the_samples <- sampling(the_model,
                        data = stan_data,
                        cores = 3,
                        chains = 3,
                        iter = 1000,
                        control = list(max_treedepth = 14),
                        pars = c("alpha", "beta",
                                 "beta_educ", "beta_age"))

