library(tidyverse)
library(ecoreg)

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

### Generate the crossing matrix that ecoreg needs
### Neither feature present, then feature 1 but not 2
### then feature 2 but not 1, then both
tmp <- expand.grid(isFemale = levels(factor(ind$isFemale)),
            isOwner = levels(factor(ind$isOwner)))

eco_cross <- sapply(unique(psw$group), function(g) { 
    tmp <- merge(tmp,
                 subset(psw, group == g),
                 by = c("isFemale", "isOwner"),
                 all.x = TRUE,
                 all.y = FALSE,
                 sort = FALSE)
    ## Return the weights
    tmp$w8
})
eco_cross <- t(eco_cross)
rownames(eco_cross) <- unique(psw$group)

stopifnot(all(rownames(eco_cross) == agg$group))

eco_model <- eco(cbind(successes, count) ~
                     p.over65 + p.ab,
                 binary =~ isFemale + isOwner,
                 iformula = y ~ p.over65 + p.ab + isFemale + isOwner,
                 cross = eco_cross,
                 data = agg,
                 idata = ind)

### ecoreg reports odds ratios
### Compare these to our known coefs
### These are taken from 02_sim_outcome.R
alpha <- -0.33
beta_female <- -0.2
beta_owner <- 0.9
beta_ab <- 0.2
beta_over65 <- 1.2

known_indiv <- c(beta_female,
                 beta_owner)
cbind(known_indiv, log(eco_model$ors.indiv))

known_agg <- c(alpha, beta_over65, beta_ab)
cbind(known_agg, log(eco_model$ors.ctx))

### 
### What happens if we put in the parameters that we get through
### fixed?

### loglik.eco <- function(pars, mod, adata, idata, gh.points=NULL, ...)
pars <- eco_model$aux$res$par
mod <- eco_model$mod
adata <- eco_model$aux$adata
idata <- eco_model$aux$idata
gh.points <- NULL

ll1 <- ecoreg:::loglik.eco(pars, mod, adata, idata)

### Now what if we supply our pars
known_pars <- c(known_agg, known_indiv)
ll2 <- ecoreg:::loglik.eco(known_pars, mod, adata, idata)


### Change reltol from its default of 1e-8
eco_model2 <- eco(cbind(successes, count) ~
                     p.over65 + p.ab,
                 binary =~ isFemale + isOwner,
                 iformula = y ~ p.over65 + p.ab + isFemale + isOwner,
                 data = agg,
                 idata = ind,
                 cross = eco_cross,
                 control = list(maxit = 2000,
                 reltol = 1e-16))

### Did we achieve convergence?
(eco_model2$aux$res$convergence == 0)

cbind(known_indiv, log(eco_model2$ors.indiv))
cbind(known_agg, log(eco_model2$ors.ctx))
