
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tvCoef

The goal of tvCoef is to study time-varying coefficients models around
linear regression:

- piecewise regressions (using
  [`strucchange`](https://CRAN.R-project.org/package=strucchange) to
  detect the break dates) ;

- local regressions (using
  [`tvReg`](https://CRAN.R-project.org/package=tvReg)) ;

- state space models (using
  [`rjd3sts`](https://github.com/palatej/rjd3sts)).

## Installation

tvCoef relies on rjd3sts which need Java JRE 17 or later version.

``` r
# install.packages("remotes")
remotes::install_github("palatej/rjd3toolkit")
remotes::install_github("palatej/rjd3sts")
remotes::install_github("AQLT/tvCoef")
```

If you have troubles installing Java, check the [installation
manual](https://github.com/jdemetra/rjdemetra/wiki/Installation-manual).
