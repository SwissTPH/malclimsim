# Loading an Odin Model

Loading an Odin Model

## Usage

``` r
load_model(name)
```

## Arguments

- name:

  name of the model written in Odin DSL to be loaded. This is a path to
  an R file,

## Value

the loaded model to be used for simulation, inference, etc

## Examples

``` r
load_model("model_new_R_with_FOI")
#> Loading required namespace: pkgbuild
#> Error: Unknown variable mu_SE
#>  update(EA) <- EA - mu_EI * EA +  mu_SE * SA - delta_d * EA + delta_a * EC # (line 43)
#>  update(EC) <- EC - mu_EI * EC +  mu_SE * SC - (delta_d + delta_a) * EC # (line 36)
```
