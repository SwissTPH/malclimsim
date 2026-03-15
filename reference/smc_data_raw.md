# Raw SMC Coverage Data from CPS

Raw seasonal malaria chemoprevention (SMC) campaign data from a cleaned
CPS Excel file covering 2018 to 2023. Includes total and card-confirmed
coverage with associated uncertainty bounds.

## Usage

``` r
smc_data_raw
```

## Format

A data frame with the following columns:

- date_start:

  Date the SMC round began (class: Date)

- epiweek_start:

  Epidemiological week the round began

- smc_couv_card:

  Card-confirmed SMC coverage

- smc_couv_card_lower:

  Lower bound of card-confirmed coverage

- smc_couv_card_upper:

  Upper bound of card-confirmed coverage

- smc_couv_tot:

  Total SMC coverage

- smc_couv_tot_lower:

  Lower bound of total coverage

- smc_couv_tot_upper:

  Upper bound of total coverage

- smc_round:

  SMC round number (1 through 4 or 5)

## Source

Extracted from a cleaned CPS Excel dataset for 2018 to 2023.
