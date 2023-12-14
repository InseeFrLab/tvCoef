

#' Computes RMSE of different models
#'
#' @description
#' Computes 6 models: linear regression, piecewise regression (with linear regression and local regression), local regression (with [tvLM]), piecewise regression with some fixed coefficients and tvlm with fixed coefficients. Computes 6 rmse on their residuals and 6 others on the residuals of the predictions of these models.
#'
#' @param x an object of class `formula` or `lm`, that is the description of the model to be fitted, or the model itself.
#' @param data a `ts` object containing the variables in the model. Necessary only when `x` is a formula.
#' @param fixed_var which variables of the model should have their coefficients fixed for the models with fixed coefficients. Obtained thanks to the hansen test.
#' @param fixed_bw `logical`, by default set to \code{FALSE}. Indicates if the bandwidth has to be computed again in the prevision model, or if it takes the value of the bandwidth of the `tvlm` model.
#' @param ... additional arguments
#'
#' @details
#' In additional arguments, parameters as date and period can be informed. As in [oos_prev] they are by default respectively set to 28 and 1.
#'
#' To estimate prevision models, the function [oos_prev] is used.
#'
#' For the previsions of the two models with fixed coefficients, fixed coefficients are re-estimated at each date, before [bp_lm] or [tvLM] are run on moving variables.
#'
#'
#' @return Returns an object of class `prev` which is a list containing the following elements:
#' \item{model}{a list of the 6 explanatory models}
#' \item{prevision}{a list of the 6 predictions of the 6 previous models}
#' \item{rmse}{a list of the 2 computed rmse, in sample and out of sample}
#'
#' @export

rmse_prev <- function(x, data, fixed_var = NULL, fixed_bw = FALSE, ...){
  UseMethod("rmse_prev", x)
}

#' @rdname rmse_prev
#' @export
rmse_prev.formula <- function(x, data, fixed_var = NULL, fixed_bw = FALSE, ...) {
  formula <- paste(deparse(x), collapse = " ")
  x_lm <- dynlm::dynlm(formula = as.formula(formula), data = data)
  data = get_data(x_lm)

  intercept <- length(grep("Intercept", names(coef(x_lm)))) > 0
  if(intercept) {
    colnames(data) = c("y", names(x_lm$coefficients)[-1])
    formule <- sprintf("%s ~ .", colnames(data)[1])
  } else {
    colnames(data) = c("y", names(x_lm$coefficients))
    formule <- sprintf("%s ~ 0 + .", colnames(data)[1])
  }

  if (is.null(fixed_var)) {
    rmse_prev_stable(x_lm = x_lm, formule = formule, data = data, fixed_bw = fixed_bw, ...)
  } else {
    rmse_prev_instable(x_lm = x_lm, formule = formule, data = data, fixed_var = fixed_var, fixed_bw = fixed_bw, ...)
  }
}

#' @rdname rmse_prev
#' @export
rmse_prev.lm <- function(x, data, fixed_var = NULL, fixed_bw = FALSE, ...) {
  # data argument not used
  x_lm <- x
  data = get_data(x_lm)

  intercept <- length(grep("Intercept", names(coef(x_lm)))) > 0
  if(intercept) {
    colnames(data) = c("y", names(x_lm$coefficients)[-1])
    formule <- sprintf("%s ~ .", colnames(data)[1])
  } else {
    colnames(data) = c("y", names(x_lm$coefficients))
    formule <- sprintf("%s ~ 0 + .", colnames(data)[1])
  }

  if (is.null(fixed_var)) {
    rmse_prev_stable(x_lm = x_lm, formule = formule, data = data, fixed_bw = fixed_bw, ...)
  } else {
    rmse_prev_instable(x_lm = x_lm, formule = formule, data = data, fixed_var = fixed_var, fixed_bw = fixed_bw, ...)
  }
}

rmse_prev_instable <- function(x_lm, formule, data, fixed_var, fixed_bw = FALSE, date = 28, break_dates = NULL, ...) {
  x_tvlm <- tvLM(formula(formule), data = data, ...)
  if(fixed_bw){b = x_tvlm$bw} else{b = NULL}
  x_piecelm <- piece_reg(x_lm, tvlm = FALSE, fixed_var = NULL, break_dates = break_dates)
  x_piecetvlm <- piece_reg(x_lm, tvlm = TRUE, fixed_var = NULL, bw = b, break_dates = break_dates)
  fixed_coef <- lm_fixed_coeff(formula = formula(formule), data = data, fixed_var = fixed_var, ...)
  x_piecelm_fixe <- piece_reg(x_lm, tvlm = FALSE, fixed_var = fixed_var, break_dates = break_dates)
  x_piecetvlm_fixe <- piece_reg(x_lm, tvlm = TRUE, fixed_var = fixed_var, bw = b, break_dates = break_dates)
  x_tvlm_fixe <- fixed_coef$tv_reg
  resid_lm <- x_lm$residuals
  resid_tvlm <- x_tvlm$residuals
  resid_tvlm_fixe <- x_tvlm_fixe$residuals
  if(inherits(x_piecelm, "lm")) {
    resid_piecelm <- x_piecelm$residuals
    resid_piecetvlm <- x_piecetvlm$residuals
  } else {
    resid_piecelm <- x_piecelm$model$residuals
    resid_piecetvlm <- x_piecetvlm$model$residuals
  }
  if(inherits(x_piecelm_fixe, "lm")) {
    resid_piecelm_fixe <- x_piecelm_fixe$residuals
    resid_piecetvlm_fixe <- x_piecetvlm_fixe$residuals
  } else {
    resid_piecelm_fixe <- x_piecelm_fixe$model$residuals
    resid_piecetvlm_fixe <- x_piecetvlm_fixe$model$residuals
  }
  prev_x_lm <- oos_prev(x_lm, date = date, data = data, ...)
  prev_x_tvlm <- oos_prev(x_tvlm, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, date = date, ...)
  prev_x_piecelm <- oos_prev(x_piecelm, date = date, data = get_data(x_piecelm), ...)
  prev_x_piecetvlm <- oos_prev(x_piecetvlm, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, date = date, ...)
  prev_x_piecelm_fixe <- oos_prev(x_piecelm_fixe, date = date, data = get_data(x_piecelm_fixe), ...)
  prev_x_piecetvlm_fixe <- oos_prev(x_piecetvlm_fixe, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, date = date, ...)
  data_est <- lapply(prev_x_lm$model, resid_lm_fixed, fixed_var = fixed_var)
  prev_x_tvlm_fixe <- oos_prev(x_tvlm_fixe, data_est = data_est, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, date = date, ...)
  prev_fixed_coef <- prev_lm_fixed(x = prev_x_lm, fixed_var = fixed_var, ...)
  y <- prev_x_lm$model[[length(prev_x_lm$model)]]$model[, 1]
  y <- y[-seq_len(prev_x_lm$debut)]
  prev_x_tvlm_fixe$prevision <- prev_x_tvlm_fixe$prevision + prev_fixed_coef$prevision
  prev_x_tvlm_fixe$residuals <- y - prev_x_tvlm_fixe$prevision
  resid_prev_lm <- prev_x_lm$residuals
  resid_prev_piecelm <- prev_x_piecelm$residuals
  resid_prev_piecetvlm <- prev_x_piecetvlm$residuals
  resid_prev_tvlm <- prev_x_tvlm$residuals
  resid_prev_tvlm_fixe <- prev_x_tvlm_fixe$residuals
  resid_prev_piecelm_fixe <- prev_x_piecelm_fixe$residuals
  resid_prev_piecetvlm_fixe <- prev_x_piecetvlm_fixe$residuals
  residus_past <- list(resid_lm, resid_piecelm, resid_piecelm_fixe, resid_tvlm, resid_piecetvlm, resid_piecetvlm_fixe, resid_tvlm_fixe)
  residus_fut <- data.frame(ts.union(
    resid_prev_lm,
    resid_prev_piecelm,
    resid_prev_piecelm_fixe,
    resid_prev_tvlm,
    resid_prev_piecetvlm,
    resid_prev_piecetvlm_fixe,
    resid_prev_tvlm_fixe
  ))
  rmse_past <- sapply(residus_past, rmse)
  rmse_fut <- sapply(residus_fut, rmse)
  names(rmse_past) <- names(rmse_fut) <-
    c("lm", "piece_lm", "piece_lm fixed coeff", "TvLM", "piece_tvlm", "piece_tvlm fixed coeff", "TvLM fixed coeff")
  res <- list(
    model = list(
      "lm" = x_lm,
      "piece_lm" = x_piecelm,
      "piece_lm_fixe" = x_piecelm_fixe,
      "tvlm" = x_tvlm,
      "piece_tvlm" = x_piecetvlm,
      "piece_tvlm_fixe" = x_piecetvlm_fixe,
      "tvlm_fixe" = x_tvlm_fixe
    ),
    prevision = list(
      "prev_lm" = prev_x_lm,
      "prev_piece_lm" = prev_x_piecelm,
      "prev_piece_lm_fixe" = prev_x_piecelm_fixe,
      "prev_tvlm" = prev_x_tvlm,
      "prev_piece_tvlm" = prev_x_piecetvlm,
      "prev_piece_tvlm_fixe" = prev_x_piecetvlm_fixe,
      "prev_tvlm_fixe" = prev_x_tvlm_fixe
    ),
    rmse = list(rmse_past, rmse_fut)
  )
  class(res) <- "prev"
  res
}



rmse_prev_stable <- function(x_lm, formule, data, fixed_bw = FALSE, date = 28, break_dates = NULL, ...) {
  x_tvlm <- tvLM(formula(formule), data = data, ...)
  if(fixed_bw){b = x_tvlm$bw} else{b = NULL}
  x_piecelm <- piece_reg(x_lm, tvlm = FALSE, fixed_var = NULL, break_dates = break_dates)
  x_piecetvlm <- piece_reg(x_lm, tvlm = TRUE, fixed_var = NULL, bw = b, break_dates = break_dates)
  resid_lm <- x_lm$residuals
  resid_tvlm <- x_tvlm$residuals
  if(inherits(x_piecelm, "lm")) {
    resid_piecelm <- x_piecelm$residuals
    resid_piecetvlm <- x_piecetvlm$residuals
  } else {
    resid_piecelm <- x_piecelm$model$residuals
    resid_piecetvlm <- x_piecetvlm$model$residuals
  }
  prev_x_lm <- oos_prev(x_lm, date = date, data = data, ...)
  prev_x_tvlm <- oos_prev(x_tvlm, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, date = date, ...)
  prev_x_piecelm <- oos_prev(x_piecelm, date = date, data = get_data(x_piecelm), ...)
  prev_x_piecetvlm <- oos_prev(x_piecetvlm, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, date = date, ...)
  resid_prev_lm <- prev_x_lm$residuals
  resid_prev_piecelm <- prev_x_piecelm$residuals
  resid_prev_piecetvlm <- prev_x_piecetvlm$residuals
  resid_prev_tvlm <- prev_x_tvlm$residuals
  residus_past <- list(resid_lm, resid_piecelm, resid_tvlm, resid_piecetvlm)
  residus_fut <- data.frame(ts.union(
    resid_prev_lm,
    resid_prev_piecelm,
    resid_prev_tvlm,
    resid_prev_piecetvlm
  ))
  rmse_past <- sapply(residus_past, rmse)
  rmse_fut <- sapply(residus_fut, rmse)
  names(rmse_past) <- names(rmse_fut) <- c("lm", "piece_lm", "TvLM", "piece_tvlm")
  res <- list(
    model = list(
      "lm" = x_lm,
      "piece_lm" = x_piecelm,
      "tvlm" = x_tvlm,
      "piece_tvlm" = x_piecetvlm
    ),
    prevision = list(
      "prev_lm" = prev_x_lm,
      "prev_piece_lm" = prev_x_piecelm,
      "prev_tvlm" = prev_x_tvlm,
      "prev_piece_tvlm" = prev_x_piecetvlm
    ),
    rmse = list(rmse_past, rmse_fut)
  )
  class(res) <- "prev"
  res
}

#' @export
print.prev <- function(x, ...) {
  df <- data.frame(
    RMSE_in_sample = x$rmse[[1]],
    RMSE_out_of_sample = x$rmse[[2]]
  )
  print(df)
}






