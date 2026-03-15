# Import an Odin Model into the malclimsim Package

Import an Odin Model into the malclimsim Package

## Usage

``` r
import_model(model_path, model_name)
```

## Arguments

- model_path:

  Full file path to the external Odin model (an R file) to be imported.

- model_name:

  Desired name (without ".R") for the model when stored inside the
  package.

## Value

None. The function copies the model file into the package's models
directory.

## Examples

``` r
import_model("path/to/external_model.R", "my_new_model")
#> [1] FALSE
```
