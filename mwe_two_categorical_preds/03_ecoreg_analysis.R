
library(tidyverse)
## library(ecoreg)
source("../eco_reg_source.R", echo = FALSE)

psw <- readRDS("working/munged_psw.rds")
ind <- readRDS("working/ind_data.rds")
agg <- readRDS("working/agg_outcomes.rds")

### We need to change the default level on age: age16-19 causes problems
ind$age <- factor(ind$age)
psw$age <- factor(psw$age)
ind$age <- relevel(ind$age, "30-44")
psw$age <- relevel(psw$age, "30-44")

### Add on summary statistics to the aggregate outcomes
agg <- merge(agg,
             psw %>%
             group_by(group) %>%
             summarize(p.ab = unique(p.ab),
                       p.owns = unique(p.owns)),
             all = TRUE,
             sort = FALSE)

### Generate the crossing matrix that ecoreg needs
### Neither feature present, then feature 1 but not 2
### then feature 2 but not 1, then both
tmp <- expand.grid(age = levels(factor(ind$age)),
            education = levels(factor(ind$education)))

eco_cross <- sapply(unique(psw$group), function(g) { 
    tmp <- merge(tmp,
                 subset(psw, group == g),
                 by = c("age", "education"),
                 all.x = TRUE,
                 all.y = FALSE,
                 sort = FALSE)
    ## Return the weights
    tmp$w8
})
eco_cross <- t(eco_cross)
rownames(eco_cross) <- unique(psw$group)

stopifnot(all(rownames(eco_cross) == agg$group))

### Generate the categorical list-of-matrices that ecoreg needs
## ### categorical: An optional list of matrices or data frames.  Each element
##           corresponds to a categorical covariate.  Each element has the
##           same number of rows as the aggregate data, and number of
##           columns corresponding to the number of levels of the
##           categorical covariate.  The cells give the number or
##           proportion of individuals in the area in each category.
##           These will be modelled as individual-level predictors of the
##           response given in formula.
age_mat <- psw %>%
    group_by(group, age) %>%
    summarize(w8 = sum(w8)) %>%
    pivot_wider(names_from = age,
                values_from = w8)

age_mat <- as.matrix(age_mat[,-1])

edu_mat <- psw %>%
    group_by(group, education) %>%
    summarize(w8 = sum(w8)) %>%
    pivot_wider(names_from = education,
                values_from = w8)

edu_mat <- as.matrix(edu_mat[,-1])

eco_categorical <- list("age" = age_mat,
                        "education" = edu_mat)

glm_mod <- glm(y ~ p.ab + p.owns + age + education,
               family = binomial,
               data = ind)

system.time(eco_model <- eco(cbind(successes, count) ~
                     p.ab + p.owns,
                 categorical = eco_categorical,
                 iformula = y ~ p.ab + p.owns + age + education,
                 data = agg,
                 idata = ind,
                 cross = eco_cross,
                 pars = coef(glm_mod),
                 control = list(maxit = 1e6,
                 reltol = 1e-16)))

### Have we achieved convergence?
(eco_model$aux$res$convergence == 0)

system.time(eco_model2 <- eco(cbind(successes, count) ~
                     p.ab + p.owns,
                 categorical = eco_categorical,
                 iformula = y ~ p.ab + p.owns + age + education,
                 data = agg,
                 idata = ind,
                 cross = eco_cross,
                 pars = my_pars,
                 control = list(maxit = 1e6,
                 reltol = 1e-16)))

### How do the coefficients compare?
alpha <- -0.33
beta_age <- seq(-2, 2, length.out = 8)
### Make the fourth group zero
beta_age <- beta_age - beta_age[4]
beta_educ <- c(seq(0, 4, length.out = 4),
               -1, 0)
beta_ab <- 0.2
beta_owns <- 1.2

cbind(c(beta_age[-4], beta_educ[-1]),
      log(eco_model$ors.indiv))

cbind(c(alpha, beta_ab, beta_owns),
      log(eco_model$ors.ctx))

cbind(c(beta_age[-4], beta_educ[-1]),
      log(eco_model2$ors.indiv))

cbind(c(alpha, beta_ab, beta_owns),
      log(eco_model2$ors.ctx))

### Let's get the agg likelihood for specified values of alpha
predprob_wrapper <- function(pars, model) {
    mod <- model$mod
    U <- rep(0, 632)
    allgroups <- seq(length = 632)
    alpha.c <- pars[1:(1+mod$nctx)] # area-level covariate effects
    alpha <- if (mod$nbineffs>0) pars[(2+mod$nctx):(1+mod$nctx+mod$nbineffs)] else NULL # individual-level binary covariate effects
    sig <- 0
    d <- 0
    adata <- model$aux$adata
    y <- adata[,"y"]
    N <- adata[,"N"]
    q <- as.numeric(adata[,mod$ctx.labs,drop=FALSE] %*% alpha.c)
    q <- q + U[match(adata[,"group"], allgroups)]*sig # input U: one per group in allgroups.  replicate to length of adata
    q <- q + adata[,"off"]
    q <- outer(q, sapply(mod$whicha, function(x) sum(alpha[x])), "+")
    p <- rowSums(adata[,mod$phi.labs,drop=FALSE] * plogis(q))
    return(list(p = p, q = q,
                alpha.c = alpha.c,
                alpha = alpha))
}

### What about if we feed it the pars?
default_pars <- eco_model$aux$res$par
my_pars <- c(-0.33, beta_ab, beta_owns,
             beta_age[-4],
             beta_educ[-1])

p1 <- predprob_wrapper(default_pars, eco_model)$p
p2 <- predprob_wrapper(my_pars, eco_model)$p

a1 <- predprob_wrapper(default_pars, eco_model)$alpha
a2 <- predprob_wrapper(my_pars, eco_model)$alpha

adata <- eco_model$aux$adata

ll1 <- sum(dbinom(adata[,"y"], adata[,"N"], p1, log=TRUE))
ll2 <- sum(dbinom(adata[,"y"], adata[,"N"], p2, log=TRUE))

mean(abs((p1 - adata[,"y"]/adata[,"N"])))
mean(abs((p2 - adata[,"y"]/adata[,"N"])))

pars <- eco_model$aux$res$par
mod <- eco_model$mod
adata <- eco_model$aux$adata
idata <- eco_model$aux$idata
gh.points <- NULL

### With the recovered pars
ll1 <- ecoreg:::loglik.eco(pars, mod, adata, idata)

### With the known pars
ll2 <- ecoreg:::loglik.eco(my_pars, mod, adata, idata)

### With the initializing pars
### This is all -2LL, so higher values = worse
ll3 <- ecoreg:::loglik.eco(eco_model$aux$pars, mod, adata, idata)

### Which component is playing up?
U <- rep(0, 632)
allgroups <- seq(length = 632)
# alpha.c = area-level covariate effects
alpha.c <- pars[1:(1+mod$nctx)] # area-level covariate effects
# alpha = individual-level binary covariate effects
alpha <- if (mod$nbineffs>0) pars[(2+mod$nctx):(1+mod$nctx+mod$nbineffs)] else NULL # individual-level binary covariate effects
# beta = individual-level normal covariate effects
beta <- NULL
#
sig <- 0
d <- 0
agglik <- lik.agg(U, adata, mod,
                           allgroups,
                           alpha.c, alpha,
                           beta, sig, d)

ilik <- lik.indiv(U, idata, mod,
                           allgroups,
                           alpha.c, alpha, beta,
                           sig, d)

sum(ilik)
sum(agglik)

### Now feed in our alternative values
alpha.c <- c(-0.33, beta_ab, beta_owns)
alpha <- c(beta_age[-4], beta_educ[-1])
beta <- NULL

agglik2 <- lik.agg(U, adata, mod,
                           allgroups,
                           alpha.c, alpha,
                           beta, sig, d)

ilik2 <- lik.indiv(U, idata, mod,
                           allgroups,
                           alpha.c, alpha, beta,
                           sig, d)

### These are LL, so higher values closer to zero = better
sum(ilik2)
sum(agglik2)

logLik(glm_mod <- glm(y ~ p.ab + p.owns + age + education,
           data = ind,
           family = binomial))



alpha.c <- c(-0.33, beta_ab, beta_owns)
alpha <- c(beta_age[-4], beta_educ[-1])

q <- as.numeric(adata[,mod$ctx.labs,drop=FALSE] %*% alpha.c)
q <- q + U[match(adata[,"group"], allgroups)]*sig # input U: one per group in allgroups.  replicate to length of adata
q <- q + adata[,"off"]

q <- outer(q, sapply(mod$whicha, function(x) sum(alpha[x])), "+")
p1 <- rowSums(adata[,mod$phi.labs,drop=FALSE] * plogis(q), na.rm = TRUE)

plot(p1, adata[,"y"]/adata[,"N"])
abline(a = 0, b= 1)

cor1 <- cor(p1, adata[,"y"]/adata[,"N"])
sum(dbinom(adata[,"y"], adata[,"N"], p1, log=TRUE))

### Huh, check this shit against doing it manually
q[1,48]
tmp[48,]
##
## The figure for q should be
alt_q <- beta_age[8] + beta_educ[6] +
    sum(adata[,mod$ctx.labs][1,] * alpha.c)



### These are different!  Okay, so the issue is with whicha, which is
### picking out coefficients 7 and 10 in "alpha", which is the last
### age category, and Level 3 for education (which is wrong) and
### indeed,there is no element of whicha which has a value of more
### than ten, when the length of alpha is 12.

###
### The relevant code in eco reg source is
### cats <- c(cats,  nlevs)
## combs <- as.matrix(expand.grid(lapply(cats, function(x)(seq(length=x)-1))))  # matrix with one row for each cross-class category
## whicha <- t(apply(combs, 1, function(x) ifelse(x>0, aoff+x, NA)))

### cats is defined earlier: if !is.null(binary), then cats <- rep(2, mod$nbin), else NULL)
### So we only need to care about nlevs
nlevs <- sapply(categorical, ncol)

### So, we have
nlevs <- c(8, 6)
cats <- nlevs

cats_seq <- lapply(cats, function(x)(seq(length=x)-1))
combs <- as.matrix(expand.grid(cats_seq))
### So far, so good
### Now we get in to the construction of whicha
aoff <- NULL
aoff <- c(aoff, mod$nbin + cumsum(nlevs-1) - (nlevs[1]-1))
aoff[is.na(aoff)] <- 0

whicha <- t(apply(combs, 1, function(x) ifelse(x>0, aoff+x, NA)))

### I think the offending line is the above one ^^^
### I don't understand the purpose of adding on a vector aoff
### If I replace this with the following
          ## whicha <- t(apply(combs, 1, function(x) ifelse(x>0, cumsum(x), NA)))               # indices of parameter vector to use 
### Then the value of q is correct
### even if the log-likelihoods are incorrect.

alpha <- log(eco_model$ors.indiv)[,1]

q <- as.numeric(adata[,mod$ctx.labs,drop=FALSE] %*% alpha.c)
q <- q + U[match(adata[,"group"], allgroups)]*sig # input U: one per group in allgroups.  replicate to length of adata
q <- q + adata[,"off"]

q <- outer(q, sapply(mod$whicha, function(x) sum(alpha[x])), "+")
p2 <- rowSums(adata[,mod$phi.labs,drop=FALSE] * plogis(q), na.rm = TRUE)


plot(p, adata[,"y"]/adata[,"N"])
cor2 <- cor(p2, adata[,"y"]/adata[,"N"])
sum(dbinom(adata[,"y"], adata[,"N"], p2, log=TRUE))

### Huh, check this shit against doing it manually
q[1,48]
tmp[48,]
##
## The figure for q should be
alt_q <- beta_age[8] + beta_educ[6] +
    sum(adata[,mod$ctx.labs][1,] * alpha.c)
