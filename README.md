
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mgCov

Multi-group Shared Subspace Covariance
Estimation

R package for the shared subspace model described in [http://jmlr.org/papers/v20/18-484.html](http://jmlr.org/papers/v20/18-484.html)

## Installation

You can install `mgcov` from github with:

``` r
# install.packages("devtools")
devtools::install_github("afranks86/mgCov")
```

## Example

The usage of `mgCov`

``` r
library(mgCov)
data(leukemia)

S <- getRank(data_list)

Vinit <- mgCov::subspaceInit(data_list, S)
EMFit <- subspaceEM(data_list, S=S)
#> [1] "Reached maximum iterations in line search."

Vfit <- EMFit$V ## inferred basis for the shared subspace
```

Now run (conditional) Bayesian covariance estimation using the inferred
subspace.

``` r
samples <- mgCov::fitBayesianSpike(V=Vfit, Ylist=data_list, 
                            niters=1000, nskip=10, verbose=FALSE)
```

``` r
groups_to_plot = c(1, 2, 4)
names(groups_to_plot) <- names(data_list)[groups_to_plot]
create_plots(V=Vfit, samples, group1=2, group2=4, to_plot = groups_to_plot, view=c(1, 2))
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

### Advanced features

TO DO
