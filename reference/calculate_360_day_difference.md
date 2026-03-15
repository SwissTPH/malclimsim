# Calculate the Day Difference Assuming a 360-Day Year

This function calculates the difference between two dates, assuming a
360-day year where each month is treated as having 30 days. This
approach is commonly used in financial calculations for simplicity.

## Usage

``` r
calculate_360_day_difference(date1, date2)
```

## Arguments

- date1:

  The start date as a character string or Date object.

- date2:

  The end date as a character string or Date object.

## Value

An integer representing the number of days between `date1` and `date2`
assuming each year has 360 days and each month has 30 days.

## Examples

``` r
# Example 1: Using character strings
calculate_360_day_difference("2023-01-01", "2023-12-31")
#> [1] 360

# Example 2: Using Date objects
start_date <- as.Date("2023-01-01")
end_date <- as.Date("2024-01-01")
calculate_360_day_difference(start_date, end_date)
#> [1] 360
```
