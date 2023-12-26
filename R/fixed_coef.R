#' Residual Effect Regression
#'
#' Checks there is no residual effect in a model
#'
#' @param x `lm` object
#' @param var the variables on which the residuals are to be regressed. By default use them all and cancel the explained variable
#'
#' @export


lm_residual_effect <- function(x, var = c(-1)) {
  resid_effect <- lm(resid(x) ~ ., data = as.data.frame(x$model[, var]))
  class(resid_effect) <- "lm"
  resid_effect
}

#' Fixed Window Regression
#'
#' @param formula a `formula` object.
#' @param data time series data.
#'
#' @returns
#'
#' Return an object of class "lmffixe".
#' Return all models, from which we can extract the usual coefficients, residuals, and fitted.values. And the divisor chosen by the function (arbitrary the middle one), the period, i.e. the length of each sub models, and the frequency of the data.
#'
#' @export


lm_fenetre_fixe <- function(formula, data, nbw = 1) {
  borne <- trunc(seq(1, nrow(data), length.out = nbw + 1))
  dataf <- lapply(1:(length(borne) - 1), function(i) {
    window(data,
      start = time(data)[borne[i]],
      end = time(data)[borne[i + 1] - deltat(data)],
      frequency = frequency(data)
    )
  })
  model <- lapply(dataf, function(data, formula) {
    dynlm::dynlm(data = data, formula = as.formula(formula))
  },
  formula = deparse(formula)
  )
  res <- list(
    model = model,
    frequency = frequency(data)
  )
  class(res) <- "lmffixe"
  res
}

#' @export

print.lmffixe <- function(x, ...) {
  colname <- sapply(1:length(x$model), function(i) {
    start <- time(x$model[[i]])[1]
    end <- time(x$model[[i]])[length(time(x$model[[i]]))]
    c(start, end)
  })
  colname2 <- sapply(1:(length(colname) / 2), function(i) {
    paste(colname[1, i], colname[2, i], sep = "-")
  })
  data <- sapply(1:length(x$model), function(i) {
    coef(x$model[[i]])
  })
  colnames(data) <- colname2
  data
}


prev_lm <- function(x) {
  predict <- sapply(x$model, predict, x$model[[length(x$model)]]$model)
  prev <- res <- vector("list", length(x$model) - 1)
  y <- x$model[[length(x$model)]]$model[, 1]
  i <- 1
  while ((i <= length(x$model)) & (x$debut + (i - 1) * x$intervalle + 1 <= nrow(predict))) {
    prev[i] <- predict[x$debut + (i - 1) * x$intervalle + 1, i]
    res[i] <- y[x$debut + (i - 1) * x$intervalle + 1] - predict[x$debut + (i - 1) * x$intervalle + 1, i]
    i <- i + 1
  }
  prev <- unlist(prev)
  prev <- ts(prev, start = x$end_dates[1] + 1 / x$frequency, frequency = x$frequency / x$intervalle)
  res <- unlist(res)
  res <- ts(res, start = x$end_dates[1] + 1 / x$frequency, frequency = x$frequency / x$intervalle)
  list(
    forecast = prev,
    residuals = res
  )
}

prev_lm_fixed <- function(x, fixed_var) {
  data_predict <- x$model[[length(x$model)]]$model
  data_predict[, -fixed_var] <- data_predict[, -fixed_var] * 0
  predict <- sapply(x$model, predict, x$model[[length(x$model)]]$model)
  predict_const <- sapply(x$model, function(mod) {
    const <- coef(mod)["(Intercept)"]
    if (all(is.na(const))) {
      const <- 0
    }
    matrix(const, nrow = nrow(data_predict))
  })
  predict <- predict - predict_const
  prev <- res <- vector("list", length(x$model) - 1)
  y <- x$model[[length(x$model)]]$model[, 1]
  i <- 1
  while ((i <= length(x$model)) & (x$debut + (i - 1) * x$intervalle + 1 <= nrow(predict))) {
    prev[i] <- predict[x$debut + (i - 1) * x$intervalle + 1, i]
    res[i] <- y[x$debut + (i - 1) * x$intervalle + 1] - predict[x$debut + (i - 1) * x$intervalle + 1, i]
    i <- i + 1
  }
  prev <- unlist(prev)
  prev <- ts(prev, start = x$end_dates[1] + 1 / x$frequency, frequency = x$frequency / x$intervalle)
  res <- unlist(res)
  res <- ts(res, start = x$end_dates[1] + 1 / x$frequency, frequency = x$frequency / x$intervalle)
  list(
    forecast = prev,
    residuals = res
  )
}

prev_tvlm <- function(x) {
  coefs <- sapply(x$model, last_coef)
  coefs[is.na(coefs)] <- 0
  mod <- x$model[[length(x$model)]]$x
  dataf <- yf <- vector("list", length(x$model) - 1)
  y <- x$model[[length(x$model)]]$y
  i <- 1
  while ((i <= length(x$model)) & (x$debut + (i - 1) * x$intervalle + 1 <= nrow(mod))) {
    dataf[[i]] <- mod[x$debut + (i - 1) * x$intervalle + 1, ]
    yf[[i]] <- y[x$debut + (i - 1) * x$intervalle + 1]
    i <- i + 1
  }
  if(is.null(ncol(coefs))){
    predict <- sapply(seq_along(dataf), function(i) {
      dataf[[i]] %*% coefs[i]
    })
  } else {
    predict <- sapply(seq_along(dataf), function(i) {
      dataf[[i]] %*% coefs[, i]
    })
  }
  res <- unlist(yf) - predict
  prev <- ts(predict, start = x$end_dates[1] + 1 / x$frequency, frequency = x$frequency / x$intervalle)
  res <- ts(res, start = x$end_dates[1] + 1 / x$frequency, frequency = x$frequency / x$intervalle)
  list(
    forecast = prev,
    residuals = res
  )
}


#' Extract data frame for lm_fixed_coeff
#'
#' @description
#' According to parameter `fixed_var`, computes a new explained variable, which is the explained variable minus the product between estimated coefficients and values of the fixed variables.
#'
#' @param x `lm` model
#' @param fixed_var list of variables that don't vary through time according to [hansen_test]
#'
#' @return
#' A new environment where the explained variable is named "fixed".
#' @export

resid_lm_fixed <- function(x, fixed_var) {
  intercept <- length(grep("Intercept", names(coef(x)))) > 0
  if (intercept) {
    data <- cbind(x$model[, 1], 1, x$model[, -1])
    colnames(data) <- c(colnames(x$model)[1], "(Intercept)", colnames(x$model)[-1])
  } else {data = x$model}
  data_y <- data[,1]
  data <- data[,-1]
  data_fixes <- data[, fixed_var, drop = FALSE]
  data_variables <- data[, -c(fixed_var), drop = FALSE]
  coef_invariant <- coef(x)[fixed_var]
  variables_fixes <- sapply(seq_len(ncol(data_fixes)), function(i) {
    data_fixes[, i] * coef_invariant[i]
  })
  fixes <- data_y - rowSums(variables_fixes, na.rm = TRUE)
  res <- cbind(fixes, data_variables)
  colnames(res) <- c("fixes", colnames(data_variables))
  res
}


#' Fixed and variable coefficients regressions
#'
#' @description
#' Computes different types of regressions with some coefficients fixed and others allowed to vary
#'
#' @param formula a `formula` object.
#' @param data time series data.
#' @param fixed_var chosen variables whose coefficients aren't allowed to vary through time
#'
#' @param ... further arguments passed to [tvReg::tvLM()]
#'
#' @return
#' \item{global_model}{the simple `lm` model}
#' \item{linear_reg}{simple `lm` model with fixed coefficients}
#' \item{piecewise_reg}{`bp_lm` model with fixed coefficients}
#' \item{tv_reg}{`tvlm` model with fixed coefficients}
#'
#' @export

lm_fixed_coeff <- function(formula, data, fixed_var, ...) {
  formula <- deparse(formula)
  x <- dynlm::dynlm(formula = formula(paste(formula, collapse = " ")), data = data)
  data_variables <- resid_lm_fixed(x, fixed_var = fixed_var)
  y_lm <- dynlm::dynlm(formula = fixes ~ -1 + ., data = data_variables)
  y_tvlm <- tvReg::tvLM(fixes ~ -1 + ., data = data_variables, ...)
  y_bplm <- bp_lm(y_lm, data_variables)
  res <- list(
    global_model = x,
    linear_reg = y_lm,
    piecewise_reg = y_bplm,
    tv_reg = y_tvlm
  )
  res
}


#' Extract Last Coefficients
#'
#' @description
#' Get last coefficients of lm or tvLM models.
#'
#' @param x a 'tvlm' or 'lm' object
#'
#' @export

last_coef <- function(x) {
  if (inherits(x, "tvlm")) {
    coef(x)[nrow(coef(x)), ]
  } else {
    coef(x)
  }
}

#' Root mean squarred error
#'
#' @param resid the residuals vector on which rmse will be calculated
#'
#' @export

rmse <- function(resid) {
  sqrt(mean(resid^2, na.rm = TRUE))
}
