
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


