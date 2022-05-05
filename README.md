
<!-- README.md is generated from README.Rmd. Please edit that file -->

# conduits (CONDitional UI for Time Series normalisation)

<img src="man/figures/logo.png" align="right" height="200"/>
<!-- badges: start --> <!-- badges: end -->

Package `conduits` provides an user interface for conditionally
normalising a time series. This also facilitates functions to produce
conditional cross-correlations between two normalised time series at
different lags while providing some graphical tools for visualisation.

`conduits` can also be used to estimate the time delay between two
sensor locations in river systems.

## Installation

You can install `conduits` from github with:

``` r
# install.packages("devtools")
devtools::install_github("PuwasalaG/conduits")
```

## Example

This is a basic example which shows you how to use functions in
`conduits`:

`conduits` contains a set of water-quality variables measured by in-situ
sensors from the Pringle Creek located in Wise County, Texas. This is
one of the aquatic NEON field sites hosted by the US Forest Service.

This data contains water-quality variables such as, turbidity, specific
conductance, dissolved oxygen, pH and fDOM along with surface elevation
and surface temperature from two sites located about 200m apart. Data
are available from
![2019-07-01](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2019-07-01 "2019-07-01")
to
![2019-12-31](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2019-12-31 "2019-12-31")
at every 5 mintues.

In this example we choose turbidity from upstream and downstream sites
to calculate the cross-correlation while conditioning on level,
temperature and conductance from the upstream location.

Let us first prepare data as follows

``` r
library(conduits)
library(lubridate)
library(tidyverse)

data_z <- NEON_PRIN_5min_cleaned %>% 
  filter(Timestamp >= ymd("2019-10-01") 
         & Timestamp < ymd("2020-01-01"),
         site == "upstream") %>%
  select(Timestamp, level, conductance, temperature) %>% 
  rename("level_up" = level,
         "conductance_up" = conductance,
         "temperature_up" = temperature)

data_xy <- NEON_PRIN_5min_cleaned %>% 
  filter(Timestamp >= ymd("2019-10-01") 
         & Timestamp < ymd("2020-01-01")) %>% 
  select(Timestamp, site, turbidity) %>% 
  tidyr::spread(key = site, value = turbidity) %>% 
  rename("turbidity_up" = upstream,
         "turbidity_down" = downstream)

data_normalisation <- data_xy %>% 
  left_join(data_z, by = "Timestamp")

head(data_normalisation)
```

### Conditional normalisation

The following code shows how to use the `conditional_moments` function
to normalise turbidity from upstream sites.

``` r
cond_moments_x <- data_normalisation %>% 
  conditional_moments(x = turbidity_up,
                      z_numeric = c(level_up, temperature_up,
                                    conductance_up),
                      knots_mean = c(8,8,8),
                      knots_variance = c(7,7,7))
head(cond_moments_x$data_conditional_moments)
```

We can then visualise the fitted models for conditional means and
variances by passing the `cond_moments_x` object to the `autoplot`
method.

``` r
autoplot(cond_moments_x, type = "mean")
```
