
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spectralem

## Installation

If `RcppFaddeeva` is not available for the current R version use:

``` r
remotes::install_version("RcppFaddeeva", "0.2.2")
```

You can install the development version with:

``` r
remotes::install_github("jan-meissner/spectralem")
```

## Examples

``` r
library(spectralem)

data <- data.synthetic.stormyclouds(seed = 1, kp = 30, noise = 0.001)
x <- data$x
y <- data$y

res <- spectralem(x, y, max_peaks = 30, max_iter = 500)
#> ================================================================================

# get position of fitted peaks
res$fit_params$pos
#>  [1] 671.1593 502.8131 326.5182 274.6381 354.3956 252.1974 313.7891 801.7046
#>  [9] 706.4931 411.1171 771.9869 187.7639 547.1771 854.1433 726.4113 572.6754
#> [17] 691.0692 516.7862 236.2330 275.0527 811.3345 303.7524 772.5411 550.5499
#> [25] 203.6653 182.2159 802.0503 801.8416 801.8487 801.8521

# plot the fit
plot_helper(x, y, res$fit)
```

<img src="man/figures/README-unnamed-chunk-2-1.svg" width="100%" />
