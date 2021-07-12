# ExploratoryDataAnalysis

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DaymondLing.github.io/ExploratoryDataAnalysis.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DaymondLing.github.io/ExploratoryDataAnalysis.jl/dev)
[![Build Status](https://github.com/DaymondLing/ExploratoryDataAnalysis.jl/workflows/CI/badge.svg)](https://github.com/DaymondLing/ExploratoryDataAnalysis.jl/actions)
[![Coverage](https://codecov.io/gh/DaymondLing/ExploratoryDataAnalysis.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/DaymondLing/ExploratoryDataAnalysis.jl)

The phrase Exploratory Data Analysis can have many meanings,
in this context, it means assess the degree of association between
a target and a large number of predictors.

`ExploratoryDataAnalysis` calculates the following metrics from
the contigency table of predictor vs. target:

- **Mutual Information**: Kullback-Liebler Divergence of the observed probabilities
from the conditionally independent probabilities constructed from the
observed row and column marginals

- **Phi coefficient**: ϕ is sqrt(χ² / n)

If target is binary,

- **Information Value**: Symmetric Kullback-Liebler Divergence between the
Class 1 distribution and Class 0 distribution

Computationally, entropy based metrics such as Mutual Information and
Information Value need to take care of 0 probabilities as they result in
Infinite entropy.
Many literature and implementations add a small positive number
to the frequency table to avoid log of 0, this is because the
software isn't capable of dealing with infinities.
Julia, however, does handle infinities gracefully, thus
this package use Infinities when there are 0 probabilities.
To avoid infinities and also not artificially adjust probabillities,
re-bin the data so that there are no 0's in the frequency table.
