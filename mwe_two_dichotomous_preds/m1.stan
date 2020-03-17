
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


