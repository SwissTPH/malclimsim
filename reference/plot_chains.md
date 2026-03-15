# Convert MCMC Output to List of Chains for Visualization

Converts the output of a multi-chain MCMC run into a
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) object,
separating chains for easier visualization and diagnostics (e.g., trace
plots or Gelman-Rubin diagnostics).

## Usage

``` r
plot_chains(mcmc_run)
```

## Arguments

- mcmc_run:

  An object returned by
  [`mcstate::pmcmc()`](https://rdrr.io/pkg/mcstate/man/pmcmc.html),
  containing posterior samples and a `chain` column indicating chain
  IDs.

## Value

An `mcmc.list` object (from the `coda` package), where each element is
an individual chain converted to an `mcmc` object.

## Details

This function assumes that `mcmc_run$pars` contains the posterior
samples and that there is a `chain` vector in `mcmc_run$chain`
specifying the chain each row belongs to.
