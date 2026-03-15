# Export Estimated and Fixed Parameters to LaTeX

Extracts posterior summaries for estimated parameters and writes both
estimated and fixed parameters as LaTeX tables.

## Usage

``` r
export_param_table_tex(
  results,
  params_to_estimate,
  out_dir,
  suffix = NULL,
  sigfig = 3
)
```

## Arguments

- results:

  A list containing inference results with `coda_pars` and
  `param_inputs`.

- params_to_estimate:

  Character vector of parameter names that were inferred.

- out_dir:

  Directory to save the LaTeX files.

- suffix:

  Optional suffix to append to filenames (e.g., dataset name).

- sigfig:

  Number of significant figures to use (default = 3).

## Value

No return value; writes LaTeX tables to disk.
