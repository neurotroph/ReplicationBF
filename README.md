# `ReplicationBF` -- R package to calculate Replication Bayes Factors

This package allows users to calculate Replication Bayes Factors for different
scenarios. See Verhagen & Wagenmakers (2014) for reference.

## Installation
The package is currently only available through this GitHub repository, so you
can install it through `devtools::install_github()`:

```R
library(devtools)
install_github('neurotroph/ReplicationBF', dependencies=TRUE)
```

## Usage
The package provides two functions to calculate Replication Bayes factors for
t- (`RBF_ttest`) and F-tests (`RBF_Ftest`) respectively.