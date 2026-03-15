# Save LaTeX Table of Scenario Estimates

Writes a LaTeX-formatted table summarizing scenario estimates to file.

## Usage

``` r
save_scenario_summary_tex(
  summary_estimates,
  out_dir,
  file_name = "scenario_summary.tex"
)
```

## Arguments

- summary_estimates:

  A data frame with columns: Scenario, Lower_2.5, Median, Upper_97.5.

- out_dir:

  Directory to save the LaTeX file.

- file_name:

  Name of the .tex file (default = "scenario_summary.tex").

## Value

No return value. Side effect: saves LaTeX file to disk.
