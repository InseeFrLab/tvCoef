#' @importFrom stats as.formula coef coefficients deltat end fitted frequency is.mts lm na.omit resid residuals start time ts ts.union tsp window is.ts
#' @importFrom stats formula
#' @importFrom utils tail
NULL

#' Break data
#'
#' @description Splits a database according to one (or more) date
#'
#' @param x a `mts` to split
#' @param break_dates the date(s) at which you want to divide the data
#' @param left `logical`. By default set to `TRUE`, i.e. the breakdate is the end date of each subcolumn
#'
#' @param ... other unused arguments
#'
#' @return a `mts` containing as many times more data columns than breakdates
#'
#' @export

break_data <- function(x, break_dates, left = TRUE, ...) {
  if (is.mts(x)) {
    name_x <- colnames(x)
  } else {
    name_x <- deparse(substitute(x))
  }
  range_time_x = tsp(x)
  break_dates = c(range_time_x[1] - deltat(x) * left,
                  break_dates,
                  range_time_x[2] + deltat(x) * !left)
  res <- do.call(ts.union, lapply(1:(length(break_dates) - 1), function(i) {
    window(x,
           start = (break_dates[i] + deltat(x) * left),
           end = (break_dates[i + 1] - deltat(x) * !left))
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
#' @param x `lm` object. It is the global regression model
#' @param break_dates optional, to indicate the breakdates if they are known. By default set to `NULL`.
#' @param fixed_var fixed variables (not splitted using `break_dates`).
#' @param left `logical`. By default set to `TRUE`, i.e. the breakdate is the end date of each submodel
#' @param tvlm By default set to `FALSE`. Indicates which model will be run on each sub data. `FALSE` means a [lm] will be run.
#' @param bw bandwidth of the local regression (when `tvlm = TRUE`).
#' @param ... other arguments passed to [tvReg::tvLM()].
#'
#' @details Computes possible breakdates if not filled in. Uses function [break_data] and run a linear regression on the same splitted data.
#'
#' @return Returns an element of class `lm`
#'
#' @export

piece_reg <- function(x, break_dates = NULL, fixed_var = NULL, tvlm = FALSE, bw = NULL, left = TRUE, ...) {
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
    break_dates <- strucchange::breakdates(strucchange::breakpoints(as.formula(formula), data = data))
  }
  if (all(is.na(break_dates))) {
    if (!tvlm) {
      return(x)
    } else {
      tv = tvReg::tvLM(formula = formula, data = data, bw = bw, ...)
      return(tv)
    }
  }
  if(is.null(fixed_var)) {
    data_break <- break_data(data[,-1], break_dates = break_dates, left = left)
    data2 <- cbind(data[,1], data_break)
    colnames(data2) <- c(colnames(data)[1], colnames(data_break))
    if(!tvlm) {
      piecereg = dynlm::dynlm(formula = as.formula(formula), data = data2)
    } else {
      piecereg = tvReg::tvLM(formula = as.formula(formula), data = data2, bw = bw, ...)
    }
  } else {
    data_x = data[,-1]
    data_break <- break_data(data_x[,- fixed_var], break_dates = break_dates, left = left)
    data2 <- cbind(data[,1], data_x[,fixed_var], data_break)
    colnames(data2) <- c(colnames(data)[1], colnames(data_x)[fixed_var], colnames(data_break))
    if(!tvlm) {
      piecereg = dynlm::dynlm(formula = as.formula(formula), data = data2)
    } else {
      piecereg = tvReg::tvLM(formula = as.formula(formula), data = data2, bw = bw, ...)
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
