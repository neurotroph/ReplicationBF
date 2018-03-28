# `ReplicationBF` -- R package to calculate Replication Bayes Factors

This package allows users to calculate Replication Bayes Factors for different
scenarios. See Verhagen & Wagenmakers (2014) and Harms (2016) for reference.

## Installation
The package is currently only available through this GitHub repository, so you
can install it through `devtools::install_github()`:

```R
devtools::install_github('neurotroph/ReplicationBF', dependencies=TRUE)
```

## Usage
The package provides two functions to calculate Replication Bayes factors for
t- (`RBF_ttest`) and F-tests (`RBF_Ftest`) respectively.

**t-Tests:**
```R
# Example 1 from Verhagen & Wagenmakers (2014)
# Using a Normal approximation to the original's posterior distribution
RBF_ttest(2.10, c(11, 11), 3.06, c(27, 28), method = "NormApprox")

# Using MCMC to draw samples from the original's posterior distribution
RBF_ttest(2.10, c(11, 11), 3.06, c(27, 28), method = "MCMC")
```

**F-Tests:**
```R
RBF_Ftest(27.0, c(3, 48), 52, 3.2, c(3, 33), 37)
```

## References
* Verhagen, J., & Wagenmakers, E.-J. (2014). Bayesian tests to quantify the result of a replication attempt. Journal of Experimental Psychology: General, 143(4), 1457–1475. http://doi.org/10.1037/a0036731
* Harms, C. (2018). A Bayes Factor for Replications of ANOVA Results. *arXiv pre-print*. Retrieved from http://arxiv.org/abs/1611.09341
