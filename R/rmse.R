

#' Computes RMSE of different models
#'
#' @description
#' Computes 6 models : linear regression, piecewise regression (with linear regression and local regression), local regression (with [tvLM]), piecewise regression with some fixed coefficients and tvlm with fixed coefficients. Computes 6 rmse on their residuals and 6 others on the residuals of the predictions of these models.
#'
#' @param formula an object of class `formula`, that is the description of the model to be fitted.
#' @param data a `ts` object containing the variables in the model
#' @param var_fixes which variables of the model should have their coefficients fixed for the models with fixed coefficients. Obtained thanks to the hansen test.
#' @param fixed_bw `logical`, by default set to \code{FALSE}. Indicates if the bandwidth has to be computed again in the prevision model, or if it takes the value of the bandwidth of the `tvlm` model.
#' @param ... additional arguments
#'
#' @details
#' In additional arguments, parameters as date and period can be informed. As in [soos_prev] they are by default respectively set to 10 and 1.
#'
#' To estimate prevision models, the function [soos_prev] is used.
#'
#' For the previsions of the two models with fixed coefficients, fixed coefficients are re-estimated at each date, before [bp.lms] or [tvLM] are run on oving variables.
#'
#'
#' @return Returns an object of class `prev` which is a list containing the following elements :
#' \item{model}{a list of the 6 explanatory models}
#' \item{prevision}{a list of the 6 predictions of the 6 previous models}
#' \item{rmse}{a list of the 2 computed rmse, in sample and out of sample}
#'
#' @export
#'

rmse_prev <- function(formula, data, var_fixes = NULL, fixed_bw = FALSE, ...) {
  if (is.null(var_fixes)) {
    rmse_prev_stable(formula = formula, data = data, fixed_bw = fixed_bw, ...)
  } else {
    rmse_prev_instable(formula = formula, data = data, var_fixes = var_fixes, fixed_bw = fixed_bw, ...)
  }
}

rmse_prev_instable <- function(formula, data, var_fixes, fixed_bw = FALSE, ...) {
  formula <- paste(deparse(formula), collapse = " ")
  x_lm <- dynlm(formula = formula(formula), data = data)
  formule <- sprintf("%s ~ .", colnames(x_lm$model)[1])
  x_piecelm <- piece_reg(x_lm, tvlm = FALSE, var_fixes = NULL, ...)
  x_piecetvlm <- piece_reg(x_lm, tvlm = TRUE, var_fixes = NULL, ...)
  x_tvlm <- tvLM(formula(formule), data = x_lm$model, ...)
  fixed_coef <- lm_fixed_coeff(formula = formula(formula), data = data, var_fixes = var_fixes, ...)
  x_piecelm_fixe <- piece_reg(x_lm, tvlm = FALSE, var_fixes = var_fixes, ...)
  x_piecetvlm_fixe <- piece_reg(x_lm, tvlm = TRUE, var_fixes = var_fixes, ...)
  x_tvlm_fixe <- fixed_coef$tv_reg
  print("ok1")
  resid_lm <- x_lm$residuals
  resid_tvlm <- x_tvlm$residuals
  resid_tvlm_fixe <- x_tvlm_fixe$residuals
  if(inherits(x_piecelm, "lm")) {
    resid_piecelm <- x_piecelm$residuals
    resid_piecetvlm <- x_piecetvlm$residuals
    prev_x_piecelm <- soos_prev(x_piecelm, ...)
    prev_x_piecetvlm <- soos_prev(x_piecetvlm, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, ...)
  } else {
    resid_piecelm <- x_piecelm$model$residuals
    resid_piecetvlm <- x_piecetvlm$model$residuals
    prev_x_piecelm <- soos_prev(x_piecelm$model, ...)
    prev_x_piecetvlm <- soos_prev(x_piecetvlm$model, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, ...)
  }
  if(inherits(x_piecelm_fixe, "lm")) {
    resid_piecelm_fixe <- x_piecelm_fixe$residuals
    resid_piecetvlm_fixe <- x_piecetvlm_fixe$residuals
    prev_x_piecelm_fixe <- soos_prev(x_piecelm_fixe, ...)
    prev_x_piecetvlm_fixe <- soos_prev(x_piecetvlm_fixe, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, ...)
  } else {
    resid_piecelm_fixe <- x_piecelm_fixe$model$residuals
    resid_piecetvlm_fixe <- x_piecetvlm_fixe$model$residuals
    prev_x_piecelm_fixe <- soos_prev(x_piecelm_fixe$model, ...)
    prev_x_piecetvlm_fixe <- soos_prev(x_piecetvlm_fixe$model, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, ...)
  }
  print("ok2")
  prev_x_lm <- soos_prev(x_lm, ...)
  prev_x_tvlm <- soos_prev(x_tvlm, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, ...)
  print("ok3")
  data_est <- lapply(prev_x_lm$model, resid_lm_fixed, var_fixes = var_fixes)
  prev_x_tvlm_fixe <- soos_prev(x_tvlm_fixe, data_est = data_est, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, ...)
  prev_fixed_coef <- prev_lm_fixed(x = prev_x_lm, var_fixes = var_fixes, ...)
  print("ok4")
  y <- prev_x_lm$model[[length(prev_x_lm$model)]]$model[, 1]
  y <- y[-seq_len(prev_x_lm$debut)]
  prev_x_tvlm_fixe$prevision <- prev_x_tvlm_fixe$prevision + prev_fixed_coef$prevision
  prev_x_tvlm_fixe$residuals <- y - prev_x_tvlm_fixe$prevision
  print("ok5")
  resid_prev_lm <- prev_x_lm$residuals
  resid_prev_piecelm <- prev_x_piecelm$residuals
  resid_prev_piecetvlm <- prev_x_piecetvlm$residuals
  resid_prev_tvlm <- prev_x_tvlm$residuals
  resid_prev_tvlm_fixe <- prev_x_tvlm_fixe$residuals
  resid_prev_piecelm_fixe <- prev_x_piecelm_fixe$residuals
  resid_prev_piecetvlm_fixe <- prev_x_piecetvlm_fixe$residuals
  print("ok6")
  residus_past <- list(resid_lm, resid_piecelm, resid_piecetvlm, resid_tvlm, resid_piecelm_fixe, resid_piecetvlm_fixe, resid_tvlm_fixe)
  residus_fut <- data.frame(ts.union(
    resid_prev_lm,
    resid_prev_piecelm,
    resid_prev_piecetvlm,
    resid_prev_tvlm,
    resid_prev_piecelm_fixe,
    resid_prev_piecetvlm_fixe,
    resid_prev_tvlm_fixe
  ))
  print("ok7")
  rmse_past <- sapply(residus_past, rmse_res)
  rmse_fut <- sapply(residus_fut, rmse_res)
  names(rmse_past) <- c("lm", "piece_lm", "piece_tvlm", "TvLM", "piece_lm fixed coeff", "piece_tvlm fixed coeff", "TvLM fixed coeff")
  names(rmse_fut) <- c("lm", "piece_lm", "piece_tvlm", "TvLM", "piece_lm fixed coeff", "piece_tvlm fixed coeff", "TvLM fixed coeff")
  print("ok8")
  res <- list(
    model = list(
      "lm" = x_lm,
      "piece_lm" = x_piecelm,
      "piece_tvlm" = x_piecetvlm,
      "tvlm" = x_tvlm,
      "piece_lm_fixe" = x_piecelm_fixe,
      "piece_tvlm_fixe" = x_piecetvlm_fixe,
      "tvlm_fixe" = x_tvlm_fixe
    ),
    prevision = list(
      "prev_lm" = prev_x_lm,
      "prev_piece_lm" = prev_x_piecelm,
      "prev_piece_tvlm" = prev_x_piecetvlm,
      "prev_tvlm" = prev_x_tvlm,
      "prev_piece_lm_fixe" = prev_x_piecelm_fixe,
      "prev_piece_tvlm_fixe" = prev_x_piecetvlm_fixe,
      "prev_tvlm_fixe" = prev_x_tvlm_fixe
    ),
    rmse = list(rmse_past, rmse_fut)
  )
  class(res) <- "prev"
  print("ok9")
  res
}


rmse_prev_stable <- function(formula, data, fixed_bw = FALSE, ...) {
    formula <- paste(deparse(formula), collapse = " ")
    x_lm <- dynlm(formula = formula(formula), data = data)
    formule <- sprintf("%s ~ .", colnames(x_lm$model)[1])
    x_piecelm <- piece_reg(x_lm, tvlm = FALSE, var_fixes = NULL)
    x_piecetvlm <- piece_reg(x_lm, tvlm = TRUE, var_fixes = NULL, ...)
    x_tvlm <- tvLM(formula(formule), data = x_lm$model, ...)
    print("ok1")
    resid_lm <- x_lm$residuals
    resid_tvlm <- x_tvlm$residuals
    print("ok2")
    prev_x_lm <- soos_prev(x_lm, ...)
    prev_x_tvlm <- soos_prev(x_tvlm, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, ...)
    if(inherits(x_piecelm, "lm")) {
      resid_piecelm <- x_piecelm$residuals
      resid_piecetvlm <- x_piecetvlm$residuals
      prev_x_piecelm <- soos_prev(x_piecelm, ...)
      prev_x_piecetvlm <- soos_prev(x_piecetvlm, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, ...)
    } else {
      resid_piecelm <- x_piecelm$model$residuals
      resid_piecetvlm <- x_piecetvlm$model$residuals
      prev_x_piecelm <- soos_prev(x_piecelm$model, ...)
      prev_x_piecetvlm <- soos_prev(x_piecetvlm$model, end = end(data), frequency = frequency(data), fixed_bw = fixed_bw, ...)
    }
    print("ok3")
    print("ok4")
    print("ok5")
    resid_prev_lm <- prev_x_lm$residuals
    resid_prev_piecelm <- prev_x_piecelm$residuals
    resid_prev_piecetvlm <- prev_x_piecetvlm$residuals
    resid_prev_tvlm <- prev_x_tvlm$residuals
    print("ok6")
    residus_past <- list(resid_lm, resid_piecelm, resid_piecetvlm, resid_tvlm)
    residus_fut <- data.frame(ts.union(
      resid_prev_lm,
      resid_prev_piecelm,
      resid_prev_piecetvlm,
      resid_prev_tvlm
    ))
    print("ok7")
    rmse_past <- sapply(residus_past, rmse_res)
    rmse_fut <- sapply(residus_fut, rmse_res)
    names(rmse_past) <- c("lm", "piece_lm", "piece_tvlm", "TvLM")
    names(rmse_fut) <- c("lm", "piece_lm", "piece_tvlm", "TvLM")
    print("ok8")
    res <- list(
      model = list(
        "lm" = x_lm,
        "piece_lm" = x_piecelm,
        "piece_tvlm" = x_piecetvlm,
        "tvlm" = x_tvlm
      ),
      prevision = list(
        "prev_lm" = prev_x_lm,
        "prev_piece_lm" = prev_x_piecelm,
        "prev_piece_tvlm" = prev_x_piecetvlm,
        "prev_tvlm" = prev_x_tvlm
      ),
      rmse = list(rmse_past, rmse_fut)
    )
    class(res) <- "prev"
    print("ok9")
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






