
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tvCoef

The goal of `tvCoef` is to study time-varying coefficients models around
linear regression:

- piecewise regressions (using
  [`strucchange`](https://CRAN.R-project.org/package=strucchange) to
  detect the break dates) ;

- local regressions (using
  [`tvReg`](https://CRAN.R-project.org/package=tvReg)) ;

- state space models (using
  [`rjd3sts`](https://github.com/rjdverse/rjd3sts)).

## Installation

`tvCoef` relies on `rjd3sts` which need Java JRE 17 or later version.

To get the current stable version (from the latest release):

``` r
# install.packages("remotes")
remotes::install_github("InseeFrLab/tvCoef@*release")
```

To get the current development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("InseeFrLab/tvCoef")
```

If you have troubles installing Java, check the [installation
manual](https://github.com/jdemetra/rjdemetra/wiki/Installation-manual).

## Usage

``` r
library(tvCoef)
data_gdp <- window(gdp, start = 1980, end = c(2019, 4))
reg_lin <- lm(
  formula = growth_gdp ~ bc_fr_m1 + diff_bc_fr_m1,
  data = data_gdp
)
reg_lin
#> 
#> Call:
#> lm(formula = growth_gdp ~ bc_fr_m1 + diff_bc_fr_m1, data = data_gdp)
#> 
#> Coefficients:
#>   (Intercept)       bc_fr_m1  diff_bc_fr_m1  
#>      -1.60008        0.02047        0.04423
```

``` r
ssm <- ssm_lm(
  reg_lin, 
  # To estimate the variance of the coefficient of the intercept            
  fixed_var_intercept = FALSE, var_intercept = 0.01,
  # To estimate the variance of the explanatory variable 
  fixed_var_variables = FALSE, var_variables = 0.01)
ssm
#> Mean of time-varying estimated coefficients (smoothing): 
#>   (Intercept)      bc_fr_m1 diff_bc_fr_m1         noise 
#>       -1.6194        0.0207        0.0425        0.0000
```

``` r
summary(ssm)
#> Summary of time-varying estimated coefficients (smoothing): 
#>         (Intercept) bc_fr_m1 diff_bc_fr_m1      noise
#> Min.         -1.791  0.02069       0.01972 -1.425e+00
#> 1st Qu.      -1.745  0.02069       0.02480 -1.997e-01
#> Median       -1.572  0.02069       0.04388  1.565e-02
#> Mean         -1.619  0.02069       0.04247  1.643e-16
#> 3rd Qu.      -1.518  0.02069       0.05731  2.243e-01
#> Max.         -1.492  0.02069       0.07697  8.114e-01
```

For more details on the methods, see the associated article available at
<https://github.com/InseeFrLab/DT-tvcoef>.
