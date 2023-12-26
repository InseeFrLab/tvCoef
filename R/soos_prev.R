
#' Out of sample forecast (or simulated out of sample)
#'
#' @param model an object used to select a method
#' @param date choose when we want to start the revision process after the start date. By default set to 28 periods.
#' @param period choose by how many values we want to move forward. By default set to 1.
#' @param ... other arguments
#'
#' @return
#' oos_prev returns an object of class `revision`, only for models of class [lm] and [tvlm]. For an object of class `bplm` it returns the same forecasts and residuals as below.
#' An object of class `revision` is a list containing the following elements:
#' \item{model}{all models used to forecast}
#' \item{debut}{same as date chosen earlier}
#' \item{intervalle}{same as period chosen earlier}
#' \item{end_dates}{a vector of all end date of each models}
#' \item{frequency}{the frequency of the data}
#' \item{forecast}{the forecast}
#' \item{residuals}{the errors of the forecast}
#'
#' @export

oos_prev <- function(model, date = 28, period = 1, ...) {
  UseMethod("oos_prev", model)
}

#' @rdname oos_prev
#' @export

oos_prev.lm <- function(model, date = 28, period = 1, data, ...) {
  # formula <- get_formula(model)
  if (missing(data)) {
    data <- get_data(model)
  }
  intercept_coef <- length(grep("Intercept", names(coef(model)))) > 0
  intercept_data <- length(grep("Intercept", colnames(data))) > 0
  if (intercept_coef & intercept_data) {
    formule <- sprintf("%s ~ 0 + .", colnames(data)[1])
  } else if (intercept_coef & !intercept_data) {
    formule <- sprintf("%s ~ .", colnames(data)[1])
  } else if (!intercept_coef & !intercept_data) {
    formule <- sprintf("%s ~ 0+ .", colnames(data)[1])
  } else {
    formule <- sprintf("%s ~ .", colnames(data)[1])
  }
  est_dates <- time(data)
  est_dates <- est_dates[seq(date, length(est_dates))]
  est_dates <- est_dates[seq(1, length(est_dates), by = period)]
  if (est_dates[length(est_dates)] != time(data)[length(time(data))]) {
    est_dates <- c(est_dates, time(data)[length(time(data))])
  }
  dataf <- lapply(est_dates, function(end) {
    suppressWarnings(window(data,
                            end = end
    ))
  })
  model <- lapply(dataf, function(data, formula) {
    lm(data = data, formula = as.formula(formula))
  },
  formula = formule
  )
  res <- list(
    model = model,
    debut = date,
    intervalle = period,
    end_dates = est_dates,
    frequency = frequency(data)
  )
  prev <- prev_lm(res)
  result <- list(
    model = res$model,
    debut = res$debut,
    intervalle = res$intervalle,
    end_dates = res$end_dates,
    frequency = res$frequency,
    forecast = prev$forecast,
    residuals = prev$residuals
  )
  class(result) <- "revision"
  result
}

#' @rdname oos_prev
#' @export

oos_prev.tvlm <- function(model, date = 28, period = 1, data_est = NULL, fixed_bw = FALSE, bw = NULL, end, frequency, ...) {
  # formula = get_formula(model)
  est <- model$est
  if (fixed_bw & is.null(bw)) {
    bw <- model$bw
  } else if (fixed_bw & !is.null(bw)) {
    bw <- bw
  } else {
    bw <- NULL
  }
  data <- get_data(model, end = end, frequency = frequency)
  intercept_coef <- length(grep("Intercept", colnames(coef(model)))) > 0
  intercept_data <- length(grep("Intercept", colnames(data))) > 0

  if (intercept_coef & intercept_data) {
    formule <- sprintf("%s ~ 0 + .", colnames(data)[1])
  } else if (intercept_coef & !intercept_data) {
    formule <- sprintf("%s ~ .", colnames(data)[1])
  } else if (!intercept_coef & !intercept_data) {
    formule <- sprintf("%s ~ 0+ .", colnames(data)[1])
  } else {
    formule <- sprintf("%s ~ .", colnames(data)[1])
  }
  est_dates <- time(data)
  est_dates <- est_dates[seq(date, length(est_dates))]
  est_dates <- est_dates[seq(1, length(est_dates), by = period)]
  if (est_dates[length(est_dates)] != time(data)[length(time(data))]) {
    est_dates <- c(est_dates, time(data)[length(time(data))])
  }
  if (is.null(data_est)) {
    data_est <- lapply(est_dates, function(end) {
      suppressWarnings(window(data,
                              end = end
      ))
    })
  } else {
    data_est <- lapply(data_est, ts, end = end, frequency = frequency)
    formule <- get_formula(model)
  }

  model <- lapply(data_est, function(data, formula) {
    tryCatch(tvReg::tvLM(data = data, formula = as.formula(formula), bw = bw, est = est),
             error = function(e) {
               tvReg::tvLM(data = data, formula = as.formula(formula), bw = NULL, est = est)
             })
  },
  formula = formule
  )
  res <- list(
    model = model,
    debut = date,
    intervalle = period,
    end_dates = est_dates,
    frequency = frequency(data)
  )
  prev <- prev_tvlm(res)
  result <- list(
    model = res$model,
    debut = res$debut,
    intervalle = res$intervalle,
    end_dates = res$end_dates,
    frequency = res$frequency,
    forecast = prev$forecast,
    residuals = prev$residuals
  )
  class(result) <- "revision"
  result
}


#' @rdname oos_prev
#' @export

oos_prev.bp_lm <- function(model, date = 28, period = 1, data_est = NULL, data, fixed_bw = FALSE, bw = NULL, ...) {
  data <- get_data(model)

  est_dates <- sapply(data, time)
  est_dates <- lapply(est_dates, function(x) {
    x <- x[seq(date, length(x))]
    x <- x[seq(1, length(x), by = period)]
    x
  })
  if(model$tvlm){
    formule <- "y ~."
    if (fixed_bw & is.null(bw)) {
      bw <- lapply(model$model, `[[`, "bw")
    } else if (fixed_bw & !is.null(bw)) {
      bw <-  rep(as.list(bw), length(est_dates))
    } else {
      bw <- rep(list(NULL), length(est_dates))
    }
  } else {
    formule <- sprintf("%s ~ .", colnames(model$model[[1]]$model)[1])
  }
  if (!is.null(data_est)) {
    line_to_remove <- c(date + 1, cumsum(sapply(data, nrow)))
    line_to_remove <- c(0, cumsum(sapply(data, nrow)) - date + 1)
    data_est <- lapply(seq_along(data), function(i) {
      data <- data_est[seq(line_to_remove[i] + 1, line_to_remove[i + 1])]
      if (i > 1) {
        data <- lapply(data, function(x) {
          x[-seq_len(line_to_remove[i] + date - 1), ]
        })
        data <- data[-seq_len(date - 1)]
      }
      data
    })
  }
  res <- lapply(seq_along(est_dates), function(i) {
    if (is.null(data_est)) {
      dataf <- lapply(est_dates[[i]], function(end) {
        data.frame(suppressWarnings(window(data[[i]],
                                           end = end
        )))
      })
    } else {
      dataf <- data_est[[i]]
    }
    if(model$tvlm) {
      model <- lapply(dataf, function(data) {
        tryCatch(tvReg::tvLM(data = data, formula = as.formula(formule), bw = bw[[i]], ...),
                 error = function(e) {
                   tvReg::tvLM(data = data, formula = as.formula(formule), bw = NULL, ...)
                 })
      })
    } else {
      model <- lapply(dataf, function(data) {
        lm(data = data, formula = as.formula(formule))
      })
    }
  })
  if(model$tvlm) {
    results <- lapply(seq_along(res), function(i) {
      prev_tvlm(list(
        model = res[[i]],
        debut = date,
        intervalle = period,
        end_dates = est_dates[[i]],
        frequency = frequency(data[[1]])
      ))
    })
  } else {
    results <- lapply(seq_along(res), function(i) {
      prev_lm(list(
        model = res[[i]],
        debut = date,
        intervalle = period,
        end_dates = est_dates[[i]],
        frequency = frequency(data[[1]])
      ))
    })
  }
  for (i in seq_len(length(res) - 1)) {
    first_date <- time(results[[i + 1]][[1]])[1]
    first_date <- first_date - 1 / frequency(results[[i + 1]][[1]])
    results[[i]] <- lapply(results[[i]], window, end = first_date, extend = TRUE)
  }
  results_ <- list(
    forecast = unlist(lapply(results, `[[`, "forecast")),
    residuals = unlist(lapply(results, `[[`, "residuals"))
  )
  results_ <- lapply(results_, ts, start = start(results[[1]][[1]]), frequency = frequency(results[[1]][[1]]))
  results_
  resultat <- list(
    model = res,
    debut = date,
    intervalle = period,
    end_dates = est_dates,
    frequency = frequency(data[[1]]),
    forecast = results_$forecast,
    residuals = results_$residuals
  )
  class(resultat) <- "revision"
  resultat
}

#' @rdname oos_prev
#' @export
oos_prev.piece_reg <- function(model, date = 28, period = 1, ...) {
  oos_prev(model$model, date = date, period = period, ...)
}

#' @export

print.revision <- function(x, ...) {
  print(list(forecast = x$forecast, residuals = x$residuals))
}
