# Convert model day index to calendar date assuming 360-day years

Converts a given day index (e.g. from a model using 360-day years) into
a calendar-like date assuming each year has 12 months of 30 days.

## Usage

``` r
model_date_360(day_index, start_year)
```

## Arguments

- day_index:

  Integer vector of day indices (starting from 1). Represents the number
  of days since the start of the model.

- start_year:

  Integer. The starting calendar year for day_index = 1.

## Value

A data frame with columns:

- `year`: The calendar year.

- `month`: The calendar month (1-12).

- `day`: The day of the month (1-30).
