# Sample Parameter Sets from MCMC Results

Randomly samples rows from a matrix or data frame containing MCMC
posterior samples. Useful for generating parameter sets to use in
posterior predictive simulations or scenario analyses.

## Usage

``` r
sample_mcmc_steps(mcmc_results, num_samples)
```

## Arguments

- mcmc_results:

  A matrix or data frame where each row represents a parameter set from
  the posterior distribution.

- num_samples:

  Integer specifying how many parameter sets to sample.

## Value

A matrix of dimension `num_samples` × `ncol(mcmc_results)`, containing
randomly selected parameter sets.

## Examples

``` r
# Simulate MCMC posterior with 100 samples and 3 parameters
posterior_samples <- matrix(rnorm(300), nrow = 100, ncol = 3)
colnames(posterior_samples) <- c("alpha", "beta", "gamma")

# Draw 10 random samples from the MCMC results
sampled <- sample_mcmc_steps(posterior_samples, 10)
print(sampled)
#>            alpha        beta       gamma
#>  [1,] -0.5351139  0.32330962  1.08950795
#>  [2,]  1.2646668 -0.87645548 -0.16686150
#>  [3,]  0.1496290 -0.06031555 -0.02083668
#>  [4,] -1.9245767 -0.22109813 -0.13654064
#>  [5,]  0.4500858  0.32809602  0.83228784
#>  [6,] -0.1057506 -2.42880752  0.75490016
#>  [7,] -0.8006491 -0.02241470 -0.50103796
#>  [8,]  0.1802437 -0.60255447 -1.60060346
#>  [9,] -0.8989590  0.40647951 -1.12716761
#> [10,]  0.2357880  0.16021413  0.70369819
```
