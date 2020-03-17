# hrrtrials

This repository contains minimal working examples of the ecoreg package for R.

The ecoreg package implements hierarchical related regression, as described in:

> "Jackson, Christopher, Nicky Best, and Sylvia
> Richardson. "Hierarchical related regression for combining aggregate
> and individual data in studies of socio‐economic disease risk
> factors." *Journal of the Royal Statistical Society: Series A*
> (Statistics in Society) 171.1 (2008): 159-178."

Directory `mwe_two_dichotomous_preds` contains an example
hierarchial related regression using, as predictors of individual and
aggregate responses, two continuous area level covariates and two
dichotomous individual level covariates.

Directory `mwe_two_categorical_preds` contains an example
hierarchial related regression using, as predictors of individual and
aggregate responses, two continuous area level covariates and two
categorical individual level covariates.

These two directories are kept separate, because the use of
categorical covariates in ecoreg requires an additional argument which
requires a bit of data munging.

Both MWE draw on common post-stratification data originally published in:

> Hanretty, Chris, Benjamin E. Lauderdale, and Nick
> Vivyan. "Comparing strategies for estimating constituency opinion
> from national survey samples." *Political Science Research and
> Methods* 6.3 (2018): 571-591."
