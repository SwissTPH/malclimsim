# Plot Correlation of MCMC Samples

This function creates a correlation plot for MCMC samples, showing the
relationships between sampled parameters. The upper panels display
correlation values, and the lower panels show smoothed plots of
parameter relationships.

## Usage

``` r
plot_corr(results, title = NULL)
```

## Arguments

- results:

  list: Results from the `inf_run` function, containing MCMC samples.

- title:

  character: Optional title for the plot.

## Value

A `GGally` ggpairs object visualizing the correlation structure of MCMC
samples.
