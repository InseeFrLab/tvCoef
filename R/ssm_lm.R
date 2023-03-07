#' State space model
#' @description
#' Computes state space model with one equation. Starting with a simple `lm` model, build the all state space model and run it.
#'
#' Use `rjd3sts` packages.
#'
#' @param model a `lm` or `dynlm` model
#' @param trend trend
#' @param var_intercept,var_slope variance of the intercept and the slope (used if `trend = TRUE`).
#' @param var_variables ??
#' @param fixed_intercept ??
#' @param fixed_variables ??.
#'
#' @return Returns a `list` containig :
#' \item{smoothed_states}{\eqn{E[a_t|y_0,\dots,y_n]}}
#' \item{smoothed_stdev}{\eqn{V[a_t|y_0,\dots,y_n]}}
#' \item{filtering_states}{\eqn{E[a_t|y_0,\dots,y_{t-1}]}}
#' \item{filtering_stdev}{\eqn{V[a_t|y_0,\dots,y_{t-1}]}}
#' \item{data}{data used in the original model}
#'
#' @export

ssm_lm <- function(x, trend = FALSE,
                   var_intercept = 0,
                   var_trend = 0,
                   var_slope = 0,
                   var_variables = 0,
                   fixed_intercept = TRUE,
                   fixed_trend = TRUE,
                   fixed_variables = TRUE, ...,
                   remove_last_dummies = FALSE) {
  UseMethod("ssm_lm", x)
}
#' @export
ssm_lm.lm <- function(x, trend = FALSE,
                      var_intercept = 0,
                      var_trend = 0,
                      var_slope = 0,
                      var_variables = 0,
                      fixed_intercept = TRUE,
                      fixed_trend = TRUE,
                      fixed_variables = TRUE, ...,
                      remove_last_dummies = FALSE) {
  ssm_lm(get_data(x),
         intercept = length(grep("Intercept", names(coef(x)))) > 0,
         trend = trend,
         var_intercept = var_intercept,
         var_trend = var_trend,
         var_slope = var_slope,
         var_variables = var_variables,
         fixed_intercept = fixed_intercept,
         fixed_trend = fixed_trend,
         fixed_variables = fixed_variables,
         remove_last_dummies = remove_last_dummies,
         ...)
}

#' @export
ssm_lm.default <- function(x,
                           intercept = TRUE,
                           trend = FALSE,
                           var_intercept = 0,
                           var_trend = 0,
                           var_slope = 0,
                           var_variables = 0,
                           fixed_intercept = TRUE,
                           fixed_trend = TRUE,
                           fixed_variables = TRUE, ...,
                           remove_last_dummies = FALSE) {
  data <- x

  if(length(var_variables) == 1)
    var_variables <- rep(var_variables, ncol(data[,-1, drop = FALSE]))[1:ncol(data[,-1, drop = FALSE])]
  if(length(fixed_variables) == 1)
    fixed_variables <- rep(fixed_variables, ncol(data[,-1, drop = FALSE]))[1:ncol(data[,-1, drop = FALSE])]
  names(var_variables) <- names(fixed_variables) <- colnames(data)[-1]

  if (remove_last_dummies) {
    # on enleve la derniere ligne car on utilise le filtering et il faut au moins 1 obs
    # apres la premiere date de l indicatrice
    data_0 = apply(data[-nrow(data),],2, function(x) all(x==0))
    data = data[, !data_0]
  }
  jmodel <- rjd3sts::model()
  jeq <- rjd3sts::equation("eq1",variance = 0, fixed = TRUE)  #ne pas modifier

  if (trend) {
    rjd3sts::add(jmodel, rjd3sts::locallineartrend("Trend",
                                                   levelVariance = var_intercept,
                                                   fixedLevelVariance = fixed_intercept,
                                                   slopevariance = var_trend, fixedSlopeVariance = fixed_trend))
    rjd3sts::add.equation(jeq, "Trend", coeff = 1, fixed = TRUE) #ne pas modifier
  } else if (intercept) {
    rjd3sts::add(jmodel, rjd3sts::locallevel("(Intercept)", variance = var_intercept, fixed = fixed_intercept))
    rjd3sts::add.equation(jeq, "(Intercept)", coeff = 1, fixed = TRUE) #ne pas modifier
  }

  for (nom_var in colnames(data)[-1]) {
    if (fixed_variables[nom_var] & var_variables[nom_var] == 0){
      rjd3sts::add(jmodel, rjd3sts::reg(nom_var, x = data[, nom_var], var = NULL, fixed = fixed_variables[nom_var]))
    } else {
      rjd3sts::add(jmodel, rjd3sts::reg(nom_var, x = data[, nom_var], var = var_variables[nom_var], fixed = fixed_variables[nom_var]))
    }
    rjd3sts::add.equation(jeq, nom_var, coeff = 1, fixed = TRUE) #ne pas modifier
  }

  rjd3sts::add(jmodel, rjd3sts::noise("noise", variance = 0.1, fixed = FALSE)) #ne pas modifier
  rjd3sts::add.equation(jeq, "noise", coeff = 1, fixed = TRUE) #ne pas modifier

  rjd3sts::add(jmodel, jeq)

  jmodestimated <- rjd3sts::estimate(jmodel, data = data[,1])

  cmp_names <- rjd3toolkit::result(jmodestimated, "ssf.cmpnames")

  if (trend) {
    cmp_names <- c("(Intercept)", cmp_names)
  }
  default_matrix <- matrix(NA, nrow = nrow(data), ncol = length(cmp_names))
  smoothed_states <- tryCatch(rjd3sts::smoothedstates(jmodestimated), error = function(e) default_matrix)
  smoothed_stdev <- tryCatch(rjd3sts::smoothedstatesstdev(jmodestimated), error = function(e) default_matrix)
  filtering_states <- tryCatch(rjd3sts::filteringstates(jmodestimated), error = function(e) default_matrix)
  filtering_stdev <- tryCatch(rjd3sts::filteringstatesstdev(jmodestimated), error = function(e) default_matrix)


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
  col_to_remove = ncol(filtering_states)
  if(trend)
    col_to_remove = c(2, col_to_remove)
  X = data[, -1, drop = FALSE]
  if (trend | intercept)
    X = cbind(1, X)
  fitted_filtering <- rowSums(X * filtering_states[,-col_to_remove])
  fitted_smoothed <- rowSums(X * smoothed_states[,-col_to_remove])
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

  res <- list(smoothed_states = smoothed_states,
              smoothed_stdev = smoothed_stdev,
              filtering_states = filtering_states,
              filtering_stdev = filtering_stdev,
              fitted = fitted,
              # parameters = list(scalingfactor = rjd3toolkit::result(jmodestimated, "scalingfactor"),
              #                   ll = rjd3toolkit::result(jmodestimated, "likelihood.ll"),
              #                   lser= rjd3toolkit::result(jmodestimated, "likelihood.ser"),
              #                   ncmp = rjd3toolkit::result(jmodestimated, "ssf.ncmps")),
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
  new_x[, c(intercept_name, names(data_0)[-1], "noise")]
}



#' Out of sample prevision of state space model
#'
#' @description
#' Computes out of sample previsions of a given state space model. Unlike [ssm_lm] it can manage dummies.
#'
#' @param model a `lm` or `dynlm` model
#' @param trend trend
#' @param var_intercept ??
#' @param var_variables ??
#' @param fixed_intercept `logical`
#' @param fixed_variables `logical`
#'
#' @return
#' Returns all coefficients of all variables and the residual
#'
#' @export
ssm_lm_oos <- function(model,
                       trend = FALSE,
                       var_intercept = 0,
                       var_trend = 0,
                       var_variables = 0,
                       fixed_intercept = TRUE,
                       fixed_variables = TRUE,
                       date = 28, ...) {
  data <- get_data(model)
  intercept <- length(grep("Intercept", names(coef(model)))) > 0
  est_data <- lapply(time(data)[-(1:date)], function(end_date) {
    window(data, end = end_date)
  })
  est_models <- lapply(est_data, ssm_lm,
                       trend = trend,
                       intercept = intercept,
                       var_intercept = var_intercept,
                       var_trend = var_trend,
                       var_variables = var_variables,
                       fixed_intercept = fixed_intercept,
                       fixed_variables = fixed_variables,
                       remove_last_dummies = TRUE)
  oos_f <- ts(t(sapply(est_models, function(x) tail(x$filtering_states,
                                                    1))),
              end = end(data), frequency = frequency(data))
  oos_f <- oos_f[, c(1:(ncol(oos_f)-1))]
  oos_noise <- data[,1] -
    ts(sapply(est_models, function(x) tail(x$fitted[,"filtering"],
                                           1)),
       end = end(data), frequency = frequency(data))
  colnames(oos_f) <- colnames(est_models[[1]]$filtering_states)[1:(length(colnames(est_models[[1]]$filtering_states))-1)]
  prevision = data[,1] - oos_noise
  res = list(oos_filtering = oos_f,
             oos_noise = oos_noise,
             prevision = prevision,
             all_models = est_models)
  res
}


#' Best state space model
#'
#' @description
#' Computes 25 state space models, using [ssm_lm] with different parameters : `fixed_intercept` and `fixed_variables` either `TRUE` or `FALSE`, and `var_intercept` and `var_variables` taking 0.01, 100 or 0 (only when associated parameters is set to `TRUE`)
#'
#' @param model a `lm` or `dynlm` model
#'
#' @return
#' Return an object of class `best_ssm`, a list containing :
#' \item{rmse_best_model_smoothed}{gives the rmse of the best model for smoothed states}
#' \item{rmse_best_model_filtering}{gives the rmse of the best model for filtering states}
#' \item{best_model_smoothed}{the entire best model for smoothed states}
#' \item{best_model_filtering}{the entire best model for filtering states}
#' \item{model}{all 25 models}
#' \item{rmse}{a `list` of 3 : all computed rmse}
#'
#' @details
#' It can give different models for smoothed and filtering states.
#'
#' @export
ssm_lm_best <- function(model) {
  var_move_1_int_fixe_0 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                          fixed_variables = FALSE,
                                          var_intercept = 0,
                                          fixed_intercept = TRUE),
                                   error = function(e) NA)
  var_move_100_int_fixe_0 = tryCatch(ssm_lm(model, var_variables = 100,
                                            fixed_variables = FALSE,
                                            var_intercept = 0,
                                            fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_1_int_fixe_0 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                          fixed_variables = TRUE,
                                          var_intercept = 0,
                                          fixed_intercept = TRUE),
                                   error = function(e) NA)
  var_fixe_100_int_fixe_0 = tryCatch(ssm_lm(model, var_variables = 100,
                                            fixed_variables = TRUE,
                                            var_intercept = 0,
                                            fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_0_int_fixe_0 = tryCatch(ssm_lm(model, var_variables = 0,
                                          fixed_variables = TRUE,
                                          var_intercept = 0,
                                          fixed_intercept = TRUE),
                                   error = function(e) NA)

  var_move_1_int_fixe_1 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                          fixed_variables = FALSE,
                                          var_intercept = 0.01,
                                          fixed_intercept = TRUE),
                                   error = function(e) NA)
  var_move_100_int_fixe_1 = tryCatch(ssm_lm(model, var_variables = 100,
                                            fixed_variables = FALSE,
                                            var_intercept = 0.01,
                                            fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_1_int_fixe_1 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                          fixed_variables = TRUE,
                                          var_intercept = 0.01,
                                          fixed_intercept = TRUE),
                                   error = function(e) NA)
  var_fixe_100_int_fixe_1 = tryCatch(ssm_lm(model, var_variables = 100,
                                            fixed_variables = TRUE,
                                            var_intercept = 0.01,
                                            fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_0_int_fixe_1 = tryCatch(ssm_lm(model, var_variables = 0,
                                          fixed_variables = TRUE,
                                          var_intercept = 0.01,
                                          fixed_intercept = TRUE),
                                   error = function(e) NA)

  var_move_1_int_fixe_100 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                            fixed_variables = FALSE,
                                            var_intercept = 100,
                                            fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_move_100_int_fixe_100 = tryCatch(ssm_lm(model, var_variables = 100,
                                              fixed_variables = FALSE,
                                              var_intercept = 100,
                                              fixed_intercept = TRUE),
                                       error = function(e) NA)
  var_fixe_1_int_fixe_100 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                            fixed_variables = TRUE,
                                            var_intercept = 100,
                                            fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_100_int_fixe_100 = tryCatch(ssm_lm(model, var_variables = 100,
                                              fixed_variables = TRUE,
                                              var_intercept = 100,
                                              fixed_intercept = TRUE),
                                       error = function(e) NA)
  var_fixe_0_int_fixe_100 = tryCatch(ssm_lm(model, var_variables = 0,
                                            fixed_variables = TRUE,
                                            var_intercept = 100,
                                            fixed_intercept = TRUE),
                                     error = function(e) NA)

  var_move_1_int_move_1 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                          fixed_variables = FALSE,
                                          var_intercept = 0.01,
                                          fixed_intercept = FALSE),
                                   error = function(e) NA)
  var_move_100_int_move_1 = tryCatch(ssm_lm(model, var_variables = 100,
                                            fixed_variables = FALSE,
                                            var_intercept = 0.01,
                                            fixed_intercept = FALSE),
                                     error = function(e) NA)
  var_fixe_1_int_move_1 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                          fixed_variables = TRUE,
                                          var_intercept = 0.01,
                                          fixed_intercept = FALSE),
                                   error = function(e) NA)
  var_fixe_100_int_move_1 = tryCatch(ssm_lm(model, var_variables = 100,
                                            fixed_variables = TRUE,
                                            var_intercept = 0.01,
                                            fixed_intercept = FALSE),
                                     error = function(e) NA)
  var_fixe_0_int_move_1 = tryCatch(ssm_lm(model, var_variables = 0,
                                          fixed_variables = TRUE,
                                          var_intercept = 0.01,
                                          fixed_intercept = FALSE),
                                   error = function(e) NA)

  var_move_1_int_move_100 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                            fixed_variables = FALSE,
                                            var_intercept = 100,
                                            fixed_intercept = FALSE),
                                     error = function(e) NA)
  var_move_100_int_move_100 = tryCatch(ssm_lm(model, var_variables = 100,
                                              fixed_variables = FALSE,
                                              var_intercept = 100,
                                              fixed_intercept = FALSE),
                                       error = function(e) NA)
  var_fixe_1_int_move_100 = tryCatch(ssm_lm(model, var_variables = 0.01,
                                            fixed_variables = TRUE,
                                            var_intercept = 100,
                                            fixed_intercept = FALSE),
                                     error = function(e) NA)
  var_fixe_100_int_move_100 = tryCatch(ssm_lm(model, var_variables = 100,
                                              fixed_variables = TRUE,
                                              var_intercept = 100,
                                              fixed_intercept = FALSE),
                                       error = function(e) NA)
  var_fixe_0_int_move_100 = tryCatch(ssm_lm(model, var_variables = 0,
                                            fixed_variables = TRUE,
                                            var_intercept = 100,
                                            fixed_intercept = FALSE),
                                     error = function(e) NA)
  mod = list(var_move_1_int_fixe_0 = var_move_1_int_fixe_0,
             var_move_100_int_fixe_0 = var_move_100_int_fixe_0,
             var_fixe_1_int_fixe_0 = var_fixe_1_int_fixe_0,
             var_fixe_100_int_fixe_0 = var_fixe_100_int_fixe_0,
             var_fixe_0_int_fixe_0 = var_fixe_0_int_fixe_0,

             var_move_1_int_fixe_1 = var_move_1_int_fixe_1,
             var_move_100_int_fixe_1 = var_move_100_int_fixe_1,
             var_fixe_1_int_fixe_1 = var_fixe_1_int_fixe_1,
             var_fixe_100_int_fixe_1 = var_fixe_100_int_fixe_1,
             var_fixe_0_int_fixe_1 = var_fixe_0_int_fixe_1,

             var_move_1_int_fixe_100 = var_move_1_int_fixe_100,
             var_move_100_int_fixe_100 = var_move_100_int_fixe_100,
             var_fixe_1_int_fixe_100 = var_fixe_1_int_fixe_100,
             var_fixe_100_int_fixe_100 = var_fixe_100_int_fixe_100,
             var_fixe_0_int_fixe_100 = var_fixe_0_int_fixe_100,

             var_move_1_int_move_1 = var_move_1_int_move_1,
             var_move_100_int_move_1 = var_move_100_int_move_1,
             var_fixe_1_int_move_1 = var_fixe_1_int_move_1,
             var_fixe_100_int_move_1 = var_fixe_100_int_move_1,
             var_fixe_0_int_move_1 = var_fixe_0_int_move_1,

             var_move_1_int_move_100 = var_move_1_int_move_100,
             var_move_100_int_move_100 = var_move_100_int_move_100,
             var_fixe_1_int_move_100 = var_fixe_1_int_move_100,
             var_fixe_100_int_move_100 = var_fixe_100_int_move_100,
             var_fixe_0_int_move_100 = var_fixe_0_int_move_100)

  mod_not_na <- sapply(mod, function(x) !all(is.na(x)))
  mod <- mod[mod_not_na]

  resid <- lapply(mod, fitted)

  smoothed_res = lapply(resid, function(x) {
    x[,"smoothed"]
  })
  filtering_res = lapply(resid, function(x) {
    x[,"filtering"]
  })

  names(smoothed_res) = names(filtering_res) = names(mod)

  rmse_smoothed = lapply(smoothed_res, rmse_res)
  rmse_filtering = lapply(filtering_res, rmse_res)
  min_rmse_smoothed = min(unlist(rmse_smoothed), na.rm = TRUE)
  min_rmse_filtering = min(unlist(rmse_filtering), na.rm = TRUE)
  rmse_best_model_smoothed = rmse_smoothed[c(match(min_rmse_smoothed, rmse_smoothed))]
  rmse_best_model_filtering = rmse_filtering[c(match(min_rmse_filtering, rmse_filtering))]
  best_model_smoothed = mod[names(rmse_best_model_smoothed)]
  best_model_filtering = mod[names(rmse_best_model_filtering)]

  res = list(rmse_best_model_smoothed = rmse_best_model_smoothed,
             rmse_best_model_filtering = rmse_best_model_filtering,
             best_model_smoothed = best_model_smoothed,
             best_model_filtering = best_model_filtering,
             model = mod,
             rmse = list(rmse_smoothed = rmse_smoothed,
                         rmse_filtering = rmse_filtering)
  )
  class(res) <- "best_ssm"
  res
}

#' Best out of sample state space model
#'
#' @description
#' Computes 25 state space models, using [ssm_lm_oos] with different parameters : `fixed_intercept` and `fixed_variables` either `TRUE` or `FALSE`, and `var_intercept` and `var_variables` taking 0.01, 100 or 0 (only when associated parameters is set to `TRUE`).
#'
#' Is is the same 25 conbinations of parameters than in [ssm_lm_best]
#'
#' @param model a `lm` or `dynlm` model
#'
#' @return
#' \item{rmse_best_model_filtering}{gives the rmse of the best model for out of sample prevision}
#' \item{best_model_filtering}{the best model}
#' \item{model}{all 25 models}
#' \item{rmse_filtering}{all 25 rmse}
#'
#' @details
#' It can give different models for smoothed and filtering states.
#'
#' @export

ssm_lm_best_oos <- function(model) {
  var_move_1_int_fixe_0 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                              fixed_variables = FALSE,
                                              var_intercept = 0,
                                              fixed_intercept = TRUE),
                                   error = function(e) NA)
  var_move_100_int_fixe_0 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                fixed_variables = FALSE,
                                                var_intercept = 0,
                                                fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_1_int_fixe_0 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                              fixed_variables = TRUE,
                                              var_intercept = 0,
                                              fixed_intercept = TRUE),
                                   error = function(e) NA)
  var_fixe_100_int_fixe_0 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                fixed_variables = TRUE,
                                                var_intercept = 0,
                                                fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_0_int_fixe_0 = tryCatch(ssm_lm_oos(model, var_variables = 0,
                                              fixed_variables = TRUE,
                                              var_intercept = 0,
                                              fixed_intercept = TRUE),
                                   error = function(e) NA)

  var_move_1_int_fixe_1 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                              fixed_variables = FALSE,
                                              var_intercept = 0.01,
                                              fixed_intercept = TRUE),
                                   error = function(e) NA)
  var_move_100_int_fixe_1 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                fixed_variables = FALSE,
                                                var_intercept = 0.01,
                                                fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_1_int_fixe_1 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                              fixed_variables = TRUE,
                                              var_intercept = 0.01,
                                              fixed_intercept = TRUE),
                                   error = function(e) NA)
  var_fixe_100_int_fixe_1 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                fixed_variables = TRUE,
                                                var_intercept = 0.01,
                                                fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_0_int_fixe_1 = tryCatch(ssm_lm_oos(model, var_variables = 0,
                                              fixed_variables = TRUE,
                                              var_intercept = 0.01,
                                              fixed_intercept = TRUE),
                                   error = function(e) NA)

  var_move_1_int_fixe_100 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                                fixed_variables = FALSE,
                                                var_intercept = 100,
                                                fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_move_100_int_fixe_100 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                  fixed_variables = FALSE,
                                                  var_intercept = 100,
                                                  fixed_intercept = TRUE),
                                       error = function(e) NA)
  var_fixe_1_int_fixe_100 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                                fixed_variables = TRUE,
                                                var_intercept = 100,
                                                fixed_intercept = TRUE),
                                     error = function(e) NA)
  var_fixe_100_int_fixe_100 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                  fixed_variables = TRUE,
                                                  var_intercept = 100,
                                                  fixed_intercept = TRUE),
                                       error = function(e) NA)
  var_fixe_0_int_fixe_100 = tryCatch(ssm_lm_oos(model, var_variables = 0,
                                                fixed_variables = TRUE,
                                                var_intercept = 100,
                                                fixed_intercept = TRUE),
                                     error = function(e) NA)

  var_move_1_int_move_1 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                              fixed_variables = FALSE,
                                              var_intercept = 0.01,
                                              fixed_intercept = FALSE),
                                   error = function(e) NA)
  var_move_100_int_move_1 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                fixed_variables = FALSE,
                                                var_intercept = 0.01,
                                                fixed_intercept = FALSE),
                                     error = function(e) NA)
  var_fixe_1_int_move_1 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                              fixed_variables = TRUE,
                                              var_intercept = 0.01,
                                              fixed_intercept = FALSE),
                                   error = function(e) NA)
  var_fixe_100_int_move_1 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                fixed_variables = TRUE,
                                                var_intercept = 0.01,
                                                fixed_intercept = FALSE),
                                     error = function(e) NA)
  var_fixe_0_int_move_1 = tryCatch(ssm_lm_oos(model, var_variables = 0,
                                              fixed_variables = TRUE,
                                              var_intercept = 0.01,
                                              fixed_intercept = FALSE),
                                   error = function(e) NA)

  var_move_1_int_move_100 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                                fixed_variables = FALSE,
                                                var_intercept = 100,
                                                fixed_intercept = FALSE),
                                     error = function(e) NA)
  var_move_100_int_move_100 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                  fixed_variables = FALSE,
                                                  var_intercept = 100,
                                                  fixed_intercept = FALSE),
                                       error = function(e) NA)
  var_fixe_1_int_move_100 = tryCatch(ssm_lm_oos(model, var_variables = 0.01,
                                                fixed_variables = TRUE,
                                                var_intercept = 100,
                                                fixed_intercept = FALSE),
                                     error = function(e) NA)
  var_fixe_100_int_move_100 = tryCatch(ssm_lm_oos(model, var_variables = 100,
                                                  fixed_variables = TRUE,
                                                  var_intercept = 100,
                                                  fixed_intercept = FALSE),
                                       error = function(e) NA)
  var_fixe_0_int_move_100 = tryCatch(ssm_lm_oos(model, var_variables = 0,
                                                fixed_variables = TRUE,
                                                var_intercept = 100,
                                                fixed_intercept = FALSE),
                                     error = function(e) NA)
  filtering = list(var_move_1_int_fixe_0 = var_move_1_int_fixe_0,
                   var_move_100_int_fixe_0 = var_move_100_int_fixe_0,
                   var_fixe_1_int_fixe_0 = var_fixe_1_int_fixe_0,
                   var_fixe_100_int_fixe_0 = var_fixe_100_int_fixe_0,
                   var_fixe_0_int_fixe_0 = var_fixe_0_int_fixe_0,

                   var_move_1_int_fixe_1 = var_move_1_int_fixe_1,
                   var_move_100_int_fixe_1 = var_move_100_int_fixe_1,
                   var_fixe_1_int_fixe_1 = var_fixe_1_int_fixe_1,
                   var_fixe_100_int_fixe_1 = var_fixe_100_int_fixe_1,
                   var_fixe_0_int_fixe_1 = var_fixe_0_int_fixe_1,

                   var_move_1_int_fixe_100 = var_move_1_int_fixe_100,
                   var_move_100_int_fixe_100 = var_move_100_int_fixe_100,
                   var_fixe_1_int_fixe_100 = var_fixe_1_int_fixe_100,
                   var_fixe_100_int_fixe_100 = var_fixe_100_int_fixe_100,
                   var_fixe_0_int_fixe_100 = var_fixe_0_int_fixe_100,

                   var_move_1_int_move_1 = var_move_1_int_move_1,
                   var_move_100_int_move_1 = var_move_100_int_move_1,
                   var_fixe_1_int_move_1 = var_fixe_1_int_move_1,
                   var_fixe_100_int_move_1 = var_fixe_100_int_move_1,
                   var_fixe_0_int_move_1 = var_fixe_0_int_move_1,

                   var_move_1_int_move_100 = var_move_1_int_move_100,
                   var_move_100_int_move_100 = var_move_100_int_move_100,
                   var_fixe_1_int_move_100 = var_fixe_1_int_move_100,
                   var_fixe_100_int_move_100 = var_fixe_100_int_move_100,
                   var_fixe_0_int_move_100 = var_fixe_0_int_move_100)

  filtering_not_na <- sapply(filtering, function(x) !all(is.na(x)))
  filtering <- filtering[filtering_not_na]
  filtering_res = lapply(filtering, residuals)
  names(filtering_res) = names(filtering)
  rmse_filtering = lapply(filtering_res, rmse_res)
  min_rmse_filtering = min(unlist(rmse_filtering), na.rm = TRUE)
  rmse_best_model_filtering = rmse_filtering[c(match(min_rmse_filtering, rmse_filtering))]
  best_model_filtering = filtering[names(rmse_best_model_filtering)]
  res = list(rmse_best_model_filtering = rmse_best_model_filtering,
             best_model_filtering = best_model_filtering,
             model = filtering,
             rmse_filtering = rmse_filtering)
  res
}

#' @export
print.best_ssm <- function(x, ...) {
  print(rbind("rmse_best_model_smoothed" = names(x$rmse_best_model_smoothed),
              "rmse_best_model_filtering" = names(x$rmse_best_model_filtering)))
}

#' Plot best_ssm
#'
#' @description
#' Plot explained variable against fitted values computed with smoothed or filtering states.
#'
#' @param x a `best_ssm` object
#' @param choice choose between `smoothed` or `filtering`
#'
#' @export
plot.best_ssm <- function(x, choice = c("smoothed", "filtering"), ...) {
  choice = match.arg(tolower(choice)[1],
                     choices = c("smoothed", "filtering"))
  best_mod = x[[sprintf("best_model_%s", choice)]]
  best_mod = best_mod[[names(best_mod)]]
  fitted =  fitted(best_mod)
  plot(cbind(best_mod$data[,1], fitted), plot.type = "s", col = c("black", "red"))
}




