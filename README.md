# malclimsim

## Installation
The package source code can be installed from GitHub using:

```
devtools::install_github("https://github.com/SwissTPH/malclimsim")
```

If there is trouble installing the package in this way, first installing two of the dependencies, *odin.dust* and *mcstate*, may help.

```
install.packages("odin.dust", repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
install.packages("mcstate", repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
devtools::install_github("https://github.com/SwissTPH/malclimsim")
```
## Purpose
The purpose of the package is to facilitate the simulation and calibration of a climate-driven model of malaria transmission as described in: 

**Nicholas Putney, Jessica Sayyad Hilario, Israel Ukawuba, Francesco
Grandesso, Saschveen Singh, Emilie Pothin, Elkoussing Djovouna,
Mahamat Saleh Issakha Diar, Clara Champagne, and Anton Camacho. *Modelling malaria routine surveillance data to inform seasonal malaria chemoprevention strategy in Moissala, Southern Chad*** (not yet published).

The general methodology described in the paper could be used to estimate SMC effectiveness and assess different strategies in other geographies. Furthermore, the package is flexibility enough to allow users to change the model and inference procedures to suit their application.

## Usage
Tutorials showing how to use the package are located at https://swisstph.github.io/malclimsim/.

## General Information
The package heavily relies on (is essentially a wrapper for) the suite of packages described in:

**FitzJohn RG, Knock ES, Whittles LK, Perez-Guzman PN, Bhatia S, Guntoro F, Watson OJ, Whittaker C, Ferguson NM, Cori A, Baguelin M, Lees JA. Reproducible parallel inference and simulation of stochastic state space models using odin, dust, and mcstate. Wellcome Open Res. 2021 Jun 10;5:288. doi: 10.12688/wellcomeopenres.16466.2. PMID: 34761122; PMCID: PMC8552050.**
