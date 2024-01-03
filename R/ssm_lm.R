#' State space model
#' @description
#' Computes state space model with one equation. Starting with a simple `lm` model, build the all state space model and run it.
#'
#' Use `rjd3sts` packages.
#'
#' @param x a `lm` or `dynlm` model
#' @param trend trend
#' @param var_intercept,var_slope variance of the intercept and the slope (used if `trend = TRUE`).
#' @param var_variables variance of the other variables: can be either a single value (same variance for all the variables) or a vector specifying each variance.
#' @param fixed_var_intercept,fixed_var_trend,fixed_var_variables logical indicating if the variance are fixed or estimated.
#'
#'
#' @param ... other arguments used in [rjd3sts::estimate()].
#' @param remove_last_dummies boolean indicating if current dummies (i.e.: only 0 and 1 at the last date) should be removed.
#'
#' @return Returns a `list` containing:
#' \item{smoothed_states}{\eqn{E[a_t|y_0,\dots,y_n]}}
#' \item{smoothed_stdev}{\eqn{\sqrt{V[a_t|y_0,\dots,y_n]}}}
#' \item{filtering_states}{\eqn{E[a_t|y_0,\dots,y_{t-1}]}}
#' \item{filtering_stdev}{\eqn{\sqrt{V[a_t|y_0,\dots,y_{t-1}]}}}
#' \item{parameters}{some estimation parameters}
#' \item{data}{data used in the original model}
#'
#' @export

ssm_lm <- function(x, trend = FALSE,
                   var_intercept = 0,
                   var_slope = 0,
                   var_variables = 0,
                   fixed_var_intercept = TRUE,
                   fixed_var_trend = TRUE,
                   fixed_var_variables = TRUE, ...,
                   remove_last_dummies = FALSE,
                   intercept = TRUE) {
  UseMethod("ssm_lm", x)
}
#' @export
ssm_lm.lm <- function(x, trend = FALSE,
                      var_intercept = 0,
                      var_slope = 0,
                      var_variables = 0,
                      fixed_var_intercept = TRUE,
                      fixed_var_trend = TRUE,
                      fixed_var_variables = TRUE, ...,
                      remove_last_dummies = FALSE,
                      intercept = has_intercept(x)) {
  ssm_lm(get_data(x),
         intercept = intercept,
         trend = trend,
         var_intercept = var_intercept,
         var_slope = var_slope,
         var_variables = var_variables,
         fixed_var_intercept = fixed_var_intercept,
         fixed_var_trend = fixed_var_trend,
         fixed_var_variables = fixed_var_variables,
         remove_last_dummies = remove_last_dummies,
         ...)
}

#' @export
ssm_lm.default <- function(x,
                           trend = FALSE,
                           var_intercept = 0,
                           var_slope = 0,
                           var_variables = 0,
                           fixed_var_intercept = TRUE,
                           fixed_var_trend = TRUE,
                           fixed_var_variables = TRUE, ...,
                           remove_last_dummies = FALSE,
                           intercept = TRUE) {
  data <- x

  if(length(var_variables) == 1)
    var_variables <- rep(var_variables, ncol(data[,-1, drop = FALSE]))[1:ncol(data[,-1, drop = FALSE])]
  if(length(fixed_var_variables) == 1)
    fixed_var_variables <- rep(fixed_var_variables, ncol(data[,-1, drop = FALSE]))[1:ncol(data[,-1, drop = FALSE])]
  names(var_variables) <- names(fixed_var_variables) <- colnames(data)[-1]

  if (remove_last_dummies) {
    # on enleve la derniere ligne car on utilise le filtering et il faut au moins 1 obs
    # apres la premiere date de l indicatrice
    data_0 = apply(data[-nrow(data),],2, function(x) all(x==0))
    data = data[, !data_0]
  }
  jmodel <- rjd3sts::model()
  # jeq <- rjd3sts::equation("eq1",variance = 0, fixed = TRUE)  #ne pas modifier

  if (trend) {
    rjd3sts::add(jmodel, rjd3sts::locallineartrend("Trend",
                                                   levelVariance = var_intercept,
                                                   fixedLevelVariance = fixed_var_intercept,
                                                   slopevariance = var_slope, fixedSlopeVariance = fixed_var_trend))
    # rjd3sts::add_equation(jeq, "Trend", coeff = 1, fixed = TRUE) #ne pas modifier
  } else if (intercept) {
    rjd3sts::add(jmodel, rjd3sts::locallevel("(Intercept)", variance = var_intercept, fixed = fixed_var_intercept))
    # rjd3sts::add_equation(jeq, "(Intercept)", coeff = 1, fixed = TRUE) #ne pas modifier
  }

  for (nom_var in colnames(data)[-1]) {
    if (fixed_var_variables[nom_var] & var_variables[nom_var] == 0){
      rjd3sts::add(jmodel, rjd3sts::reg(nom_var, x = data[, nom_var], var = NULL))
    } else {
      rjd3sts::add(jmodel, rjd3sts::reg(nom_var, x = data[, nom_var], var = var_variables[nom_var], fixed = fixed_var_variables[nom_var]))
    }
    # rjd3sts::add_equation(jeq, nom_var, coeff = 1, fixed = TRUE) #ne pas modifier
  }

  rjd3sts::add(jmodel, rjd3sts::noise("noise", variance = 0.1, fixed = FALSE)) #ne pas modifier
  # rjd3sts::add_equation(jeq, "noise", coeff = 1, fixed = TRUE) #ne pas modifier

  # rjd3sts::add(jmodel, jeq)

  jmodestimated <- rjd3sts::estimate(jmodel, data = data[,1], ...)

  cmp_names <- rjd3toolkit::result(jmodestimated, "ssf.cmpnames")

  if (trend) {
    cmp_names <- c("(Intercept)", cmp_names)
  }
  default_matrix <- matrix(NA, nrow = nrow(data), ncol = length(cmp_names))
  smoothed_states <- tryCatch(rjd3sts::smoothed_states(jmodestimated), error = function(e) default_matrix)
  smoothed_stdev <- tryCatch(rjd3sts::smoothed_states_stdev(jmodestimated), error = function(e) default_matrix)
  filtering_states <- tryCatch(rjd3sts::filtering_states(jmodestimated), error = function(e) default_matrix)
  filtering_stdev <- tryCatch(rjd3sts::filtering_states_stdev(jmodestimated), error = function(e) default_matrix)


  colnames(smoothed_states) <- colnames(smoothed_stdev) <-
    colnames(filtering_states) <-
    colnames(filtering_stdev) <-
    cmp_names

  if (is.ts(data)) {
    smoothed_states <- ts(smoothed_states, end = end(data), frequency = frequency(data))
    smoothed_stdev <- ts(smoothed_stdev, end = end(data), frequency = frequency(data))
    filtering_states <- ts(filtering_states, end = end(data), frequency = frequency(data))
    filtering_stdev <- ts(filtering_stdev, end = end(data), frequency = frequency(data))
  }
  col_to_keep <- seq_len(ncol(filtering_states) - 1)
  if (trend)
    col_to_keep <- col_to_keep[-2]
  X = data[, -1, drop = FALSE]
  if (trend | intercept)
    X = cbind(1, X)
  fitted_filtering <- rowSums(X * filtering_states[, col_to_keep])
  fitted_smoothed <- rowSums(X * smoothed_states[, col_to_keep])
  fitted <- cbind(fitted_smoothed, fitted_filtering)
  colnames(fitted) <- c("smoothed", "filtering")
  if (is.ts(data)) {
    fitted <- ts(fitted, end = end(data), frequency = frequency(data))
  }

  if (remove_last_dummies && any(data_0)) {
    smoothed_states <- add_removed_var(smoothed_states, data_0, intercept, trend)
    smoothed_stdev <- add_removed_var(smoothed_stdev, data_0, intercept, trend)
    filtering_states <- add_removed_var(filtering_states, data_0, intercept, trend)
    filtering_stdev <- add_removed_var(filtering_stdev, data_0, intercept, trend)
  }

  parameters <- rjd3toolkit::result(jmodestimated, "parameters")
  names(parameters) <- rjd3toolkit::result(jmodestimated, "parametersnames")
  scaling_factor <- rjd3toolkit::result(jmodestimated, "scalingfactor")

  res <- list(smoothed_states = smoothed_states,
              smoothed_stdev = smoothed_stdev,
              filtering_states = filtering_states,
              filtering_stdev = filtering_stdev,
              fitted = fitted,
              parameters = list(parameters = parameters,
                                scaling = scaling_factor,
                                ll = rjd3toolkit::result(jmodestimated, "likelihood.ll"),
                                lser= rjd3toolkit::result(jmodestimated, "likelihood.ser"),
                                ncmp = rjd3toolkit::result(jmodestimated, "ssf.ncmps")),
              data = data)
  class(res) <- "ssm_lm"
  res
}

#' @export
fitted.ssm_lm <- function(object, ...) {
  object$fitted
}

#' @export
residuals.ssm_lm <- function(object, ...) {
  object$data[,1] - object$fitted
}

#' @export
summary.ssm_lm <- function(object, digits = max(3, getOption("digits") - 3),
                           ...) {
  cat("Summary of time-varying estimated coefficients (smoothing):", "\n")
  coef <- object$smoothed_states
  noise <- grep("^noise$", colnames(coef))
  if (length(noise) > 0) {
    coef <- coef[,-noise, drop = FALSE]
  }
  print(apply(object$smoothed_states, 2, summary), digits = digits)
  invisible(object)
}

#' @export
print.ssm_lm <- function(object, digits = max(3, getOption("digits") - 3),
                           ...) {
  cat("Mean of time-varying estimated coefficients (smoothing):", "\n")
  coef <- object$smoothed_states
  noise <- grep("^noise$", colnames(coef))
  if (length(noise) > 0) {
    coef <- coef[,-noise, drop = FALSE]
  }
  print(round(apply(object$smoothed_states, 2, mean, na.rm = TRUE), digits), digits = digits)
  invisible(object)
}


#' @export
coef.ssm_lm <- function(object, digits = max(3, getOption("digits") - 3),
                           ...) {
  coef <- object$smoothed_states
  noise <- grep("^noise$", colnames(coef))
  if (length(noise) > 0) {
    coef <- coef[,-noise, drop = FALSE]
  }
  coef
}

# internal function

add_removed_var <- function(x, data_0, intercept, trend = FALSE) {
  new_x <- cbind(x, matrix(0, ncol = sum(data_0), nrow = nrow(x)))
  colnames(new_x) <- c(colnames(x), names(which(data_0)))
  if (trend) {
    intercept_name <- colnames(x)[1:2]
  } else if (intercept) {
    intercept_name <- colnames(x)[1]
  } else {
    intercept_name <- NULL
  }
  new_x[, c(intercept_name, names(data_0)[-1])]
}



#' Out of sample forecast of state space model
#'
#' @description
#' Computes out of sample forecasts of a given state space model. Unlike [ssm_lm] it can manage dummies.
#'
#' @inheritParams ssm_lm
#' @inheritParams oos_prev
#'
#' @return
#' Returns all coefficients of all variables and the residual
#'
#' @export
ssm_lm_oos <- function(x,
                       trend = FALSE,
                       var_intercept = 0,
                       var_slope = 0,
                       var_variables = 0,
                       fixed_var_intercept = TRUE,
                       fixed_var_trend = TRUE,
                       fixed_var_variables = TRUE,
                       date = 28, ...) {
  data <- get_data(x)
  intercept <- has_intercept(x)
  est_data <- lapply(time(data)[-(1:date)], function(end_date) {
    window(data, end = end_date)
  })
  est_models <- lapply(est_data, ssm_lm,
                       trend = trend,
                       intercept = intercept,
                       var_intercept = var_intercept,
                       var_slope = var_slope,
                       var_variables = var_variables,
                       fixed_var_intercept = fixed_var_intercept,
                       fixed_var_trend = fixed_var_trend,
                       fixed_var_variables = fixed_var_variables,
                       remove_last_dummies = TRUE)
  oos_f <- ts(t(sapply(est_models, function(x) tail(x$filtering_states,
                                                    1))),
              end = end(data), frequency = frequency(data))
  oos_f <- oos_f[,-ncol(oos_f)]
  oos_noise <- data[,1] -
    ts(sapply(est_models, function(x) tail(x$fitted[,"filtering"],
                                           1)),
       end = end(data), frequency = frequency(data))
  colnames(oos_f) <- colnames(est_models[[1]]$filtering_states)[1:(length(colnames(est_models[[1]]$filtering_states))-1)]
  forecast = data[,1] - oos_noise
  res = list(oos_filtering = oos_f,
             oos_noise = oos_noise,
             forecast = forecast,
             all_models = est_models)
  res
}


#' Check if model has intercept
#'
#' @param x a model
#'
#' @export
has_intercept <- function(x) {
  UseMethod("has_intercept", x)
}
#' @export
has_intercept.lm <- function(x) {
  length(grep("Intercept", names(coef(x)))) > 0
}

predict.ssm_lm <- function(object, newdata){

}
