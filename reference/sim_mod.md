# Simulate the Malaria Model Over Time

Runs a stochastic malaria transmission model using an `odin`-generated
model object, storing state outputs across multiple time steps and
particles.

## Usage

``` r
sim_mod(odin_mod, pars, time_start, n_particles, sim_time)
```

## Arguments

- odin_mod:

  An `odin` model object (usually created with
  [`odin::odin()`](https://rdrr.io/pkg/odin/man/odin.html)),
  representing the compiled malaria model.

- pars:

  A named list of parameters to pass to the model.

- time_start:

  Numeric value specifying the starting time (e.g., 0).

- n_particles:

  Integer; number of particles used in the stochastic simulation (i.e.,
  number of replicate simulations).

- sim_time:

  Integer; total number of time steps (days) to simulate.

## Value

A list with two elements:

- x:

  A 3D array of simulation results with dimensions
  `[state, particle, time]`.

- model:

  The `odin` model object after simulation (can be used to inspect
  internals).

## Details

The model is initialized once using the `odin_mod$new()` constructor,
and then `model$run(t)` is called at each time step `t` to advance the
simulation.

## Examples

``` r
if (FALSE) { # \dontrun{
mod <- odin::odin("malaria_model.R")
results <- sim_mod(mod, pars = list(beta = 0.2), time_start = 0, n_particles = 10, sim_time = 100)
} # }
```
