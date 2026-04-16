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
#> Unused equations: dt, size, size_1, size_2
#>  dt <- 1 / steps_per_day # (line 7)
#>  size <- user(10) # (line 178)
#>  size_1 <- user(10) # (line 179)
#>  size_2 <- user(10) # (line 180)
#> ℹ 21 functions decorated with [[cpp11::register]]
#> ✔ generated file cpp11.R
#> ✔ generated file cpp11.cpp
#> ℹ Re-compiling dust4d72123d
#> ── R CMD INSTALL ───────────────────────────────────────────────────────────────
#> * installing *source* package ‘dust4d72123d’ ...
#> ** this is package ‘dust4d72123d’ version ‘0.0.1’
#> ** using staged installation
#> ** libs
#> using C++ compiler: ‘g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0’
#> g++ -std=gnu++17 -I"/opt/R/4.5.3/lib/R/include" -DNDEBUG  -I'/home/runner/work/_temp/Library/cpp11/include' -I/usr/local/include   -I"/home/runner/work/_temp/Library/dust/include" -DHAVE_INLINE -fopenmp  -fpic  -g -O2  -Wall -pedantic -fdiagnostics-color=always  -c cpp11.cpp -o cpp11.o
#> g++ -std=gnu++17 -I"/opt/R/4.5.3/lib/R/include" -DNDEBUG  -I'/home/runner/work/_temp/Library/cpp11/include' -I/usr/local/include   -I"/home/runner/work/_temp/Library/dust/include" -DHAVE_INLINE -fopenmp  -fpic  -g -O2  -Wall -pedantic -fdiagnostics-color=always  -c dust.cpp -o dust.o
#> g++ -std=gnu++17 -shared -L/opt/R/4.5.3/lib/R/lib -L/usr/local/lib -o dust4d72123d.so cpp11.o dust.o -fopenmp -L/opt/R/4.5.3/lib/R/lib -lR
#> installing to /tmp/Rtmp07tfSS/devtools_install_20bd76b1021b/00LOCK-file20bd1f45d3c2/00new/dust4d72123d/libs
#> ** checking absolute paths in shared objects and dynamic libraries
#> * DONE (dust4d72123d)
#> ℹ Loading dust4d72123d
#> <dust> object generator
#>   Public:
#>     initialize: function (pars, time, n_particles, n_threads = 1L, seed = NULL, 
#>     name: function () 
#>     param: function () 
#>     run: function (time_end) 
#>     simulate: function (time_end) 
#>     run_adjoint: function () 
#>     set_index: function (index) 
#>     index: function () 
#>     ode_control: function () 
#>     ode_statistics: function () 
#>     n_threads: function () 
#>     n_state: function () 
#>     n_particles: function () 
#>     n_particles_each: function () 
#>     shape: function () 
#>     update_state: function (pars = NULL, state = NULL, time = NULL, set_initial_state = NULL, 
#>     state: function (index = NULL) 
#>     time: function () 
#>     set_stochastic_schedule: function (time) 
#>     reorder: function (index) 
#>     resample: function (weights) 
#>     info: function () 
#>     pars: function () 
#>     rng_state: function (first_only = FALSE, last_only = FALSE) 
#>     set_rng_state: function (rng_state) 
#>     has_openmp: function () 
#>     has_gpu_support: function (fake_gpu = FALSE) 
#>     has_compare: function () 
#>     real_size: function () 
#>     time_type: function () 
#>     rng_algorithm: function () 
#>     uses_gpu: function (fake_gpu = FALSE) 
#>     n_pars: function () 
#>     set_n_threads: function (n_threads) 
#>     set_data: function (data, shared = FALSE) 
#>     compare_data: function () 
#>     filter: function (time_end = NULL, save_trajectories = FALSE, time_snapshot = NULL, 
#>     gpu_info: function () 
#>     transform_variables: function (y) 
#>   Private:
#>     pars_: NULL
#>     pars_multi_: NULL
#>     index_: NULL
#>     info_: NULL
#>     n_threads_: NULL
#>     n_particles_: NULL
#>     n_particles_each_: NULL
#>     shape_: NULL
#>     ptr_: NULL
#>     gpu_config_: NULL
#>     ode_control_: NULL
#>     methods_: NULL
#>     param_: list
#>     reload_: list
#>   Parent env: <environment: namespace:dust4d72123d>
#>   Locked objects: TRUE
#>   Locked class: FALSE
#>   Portable: TRUE
```
