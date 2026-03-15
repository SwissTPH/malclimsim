# Filter a dataset by a range of years

This function filters a dataset to include only rows where the year
extracted from a specified date column falls within a given range of
years. The date column should be in the format "YYYY-MM-DD" as a string.

## Usage

``` r
filter_by_year(data, date_column, years_range)
```

## Arguments

- data:

  A data frame containing the dataset to filter.

- date_column:

  A string specifying the name of the column that contains dates as
  strings.

- years_range:

  A numeric vector specifying the range of years to include in the
  filtered dataset.

## Value

A filtered data frame containing only rows within the specified year
range.

## Examples

``` r
# Example dataset
met_360 <- data.frame(
  date = c("2012-01-01", "2015-06-15", "2020-03-22"),
  temp = c(26.3, 27.5, 28.1)
)

# Filter data for the years 2014 to 2023
years_analysis <- 2014:2023
met_360_filtered <- filter_by_year(met_360, "date", years_analysis)
```
