library(tidyverse)
library(ecoreg)
## source("../eco_reg_source.R", echo = FALSE)

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

### Let's try again

### How do the coefficients compare?
my_alpha <- -0.33
beta_age <- seq(-2, 2, length.out = 8)
### Make the fourth group zero
beta_age <- beta_age - beta_age[4]
beta_educ <- c(seq(0, 4, length.out = 4),
               -1, 0)
beta_ab <- 0.2
beta_owns <- 1.2
my_pars <- c(my_alpha, beta_ab, beta_owns,
             beta_age, beta_educ)

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

cbind(c(beta_age[-4], beta_educ[-1]),
      log(eco_model$ors.indiv))

cbind(c(alpha, beta_ab, beta_owns),
      log(eco_model$ors.ctx))

cbind(c(beta_age[-4], beta_educ[-1]),
      log(eco_model2$ors.indiv))

cbind(c(my_alpha, beta_ab, beta_owns),
      log(eco_model2$ors.ctx))

### Huh, neither works. Why is this?
### What is the maximum parameter that can be picked out by whicha?
max(unlist(eco_model$mod$whicha))
### What is the length of our parameters?
length(eco_model$ors.indiv)
