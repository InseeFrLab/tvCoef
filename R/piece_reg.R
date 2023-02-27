#' @importFrom stats as.formula coef coefficients deltat end fitted frequency is.mts lm na.omit resid residuals start time ts ts.union tsp window
NULL

#' Break data
#'
#' @description Splits a database according to one (or more) date
#'
#' @param x a `mts` to split
#' @param break_dates the date(s) at which you want to divide the data
#' @param right `logical`. By default set to `TRUE`, i.e. the breakdate is the end date of each subcolumn
#'
#' @return a `mts` containing as many times more data columns than breakdates
#'
#' @export

break_data <- function(x, break_dates, right = TRUE, ...) {
  if (is.mts(x)) {
    name_x <- colnames(x)
  } else {
    name_x <- deparse(substitute(x))
  }
  range_time_x = tsp(x)
  break_dates = c(range_time_x[1] - deltat(x) * right,
                  break_dates,
                  range_time_x[2] + deltat(x) * !right)
  res <- do.call(ts.union, lapply(1:(length(break_dates) - 1), function(i) {
    window(x,
           start = (break_dates[i] + deltat(x) * right),
           end = (break_dates[i + 1] - deltat(x) * !right))
  }))
  res[is.na(res)] <- 0
  colnames(res) <- sprintf("%s_%s",
                           rep(name_x, length(break_dates) - 1),
                           rep(break_dates[-1], each = length(name_x))
  )
  res
}

#' Piecewise regression
#'
#' @description Computes one global linear regression, on splitted data
#'
#' @param x `lm` object, of which we will separate the data
#' @param break_dates optional, to indicate the breakdates if they are known. By default set to `NULL`
#'
#' @details Computes possible breakdates if not filled in. Uses function [break_data] and run a linear regression on the same splitted data.
#'
#' @return Returns an element of class `lm`
#'
#' @export

piece_reg <- function(x, break_dates = NULL, var_fixes = NULL, tvlm = FALSE, bw = NULL,...) {
  data <- get_data(x)
  intercept <- length(grep("Intercept", names(coef(x)))) > 0
  if (intercept) {
    data_ <- cbind(data[, 1], 1, data[, -1])
    colnames(data_) <- c("y", "(Intercept)", colnames(data)[-1])
    data <- data_
  } else {
    colnames(data) <- c("y", colnames(data)[-1])
  }
  formula <- sprintf("%s ~ 0 + .", colnames(data)[1])

  if (is.null(break_dates)) {
    break_dates <- breakdates(breakpoints(as.formula(formula), data = data))
  }
  if (all(is.na(break_dates))) {
    if (!tvlm) {
      return(x)
    } else {
      tv = tvLM(formula = formula, data = data, bw = bw, ...)
      return(tv)
    }
  }
  if(is.null(var_fixes)) {
    data_break <- break_data(data[,-1], break_dates = break_dates)
    data2 <- cbind(data[,1], data_break)
    colnames(data2) <- c(colnames(data)[1], colnames(data_break))
    if(!tvlm) {
      piecereg = dynlm(formula = as.formula(formula), data = data2)
    } else {
      piecereg = tvLM(formula = as.formula(formula), data = data2, bw = bw, ...)
    }
  } else {
    data_x = data[,-1]
    data_break <- break_data(data_x[,- var_fixes], break_dates = break_dates)
    data2 <- cbind(data[,1], data_x[,var_fixes], data_break)
    colnames(data2) <- c(colnames(data)[1], colnames(data_x)[var_fixes], colnames(data_break))
    if(!tvlm) {
      piecereg = dynlm(formula = as.formula(formula), data = data2)
    } else {
      piecereg = tvLM(formula = as.formula(formula), data = data2, bw = bw, ...)
    }
  }
  res <- list(
    model = piecereg,
    start = start(data),
    end = end(data),
    frequency = frequency(data),
    breakdates = break_dates,
    tvlm = tvlm
  )
  class(res) <- "piecereg"
  res
}


print.piecereg <- function(x, ...) {
  if(x$tvlm) {
    print(x$model)
  } else {
    print(summary(x$model))
  }
}
