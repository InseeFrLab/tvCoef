#' Extract Formula From model
#'
#' @param x the model.
#'
#' @export
get_formula <- function(x) {
  UseMethod("get_formula", x)
}

#' @rdname get_formula
#' @export

get_formula.default <- function(x) {
  as.formula(x[["call"]][["formula"]])
}

#' Get data function
#'
#' @description Retrieves the data used in the model
#'
#' @param model the model whose data we want

#' @export

get_data <- function(model, ...) {
  UseMethod("get_data", model)
}

#' @rdname get_data
#' @export
get_data.lm <- function(model, start, ...) {
  if (missing(start))
    start <- as.numeric(rownames(model$model)[1])
  ts(model$model, start = start,...)
}

#' @rdname get_data
#' @export
get_data.dynlm <- function(model, ...) {
  dates <- sapply(strsplit(rownames(model$model), " "), function(x) {
    gsub("\\D", "", x)
  })
  frequency <- max(table(dates[1, ]))
  last_year <- dates[1, ncol(dates)]
  last_period <- dates[2, ncol(dates)]
  ts(model$model, end = as.numeric(c(last_year, last_period)), frequency = frequency)
}

#' @rdname get_data
#' @export
get_data.tvlm <- function(model, end, frequency, ...) {
  data <- cbind(model$y, model$x)
  # gsub(" ~.*", "", deparse(model$call$formula)) # plutot que y
  colnames(data) <- c("y", colnames(model$x))
  if (colnames(data)[2] == "(Intercept)") {
    data <- data[, -2]
  }
  ts(data, end = end, frequency = frequency)
}

#' @rdname get_data
#' @export
get_data.bp.lms <- function(model, ...) {
  if (model$tvlm) {
    data <- sapply(seq_along(model$model), function(i) {
      ts(data.frame(model$model[[i]]$y, model$model[[i]]$x),
         end = model$breakdates[i + 1],
         frequency = model$frequency)
    })
    for(i in seq_along(data)){
      colnames(data[[i]]) = c("y", colnames(model$model[[i]]$x))
    }
    if (colnames(data[[1]])[2] == "(Intercept)") {
      for(i in seq_along(data_x)) {
        data[[i]] <- data[[i]][, -2]
      }
    }
  } else {
    data <- sapply(seq_along(model$model), function(i) {
      ts(model$model[[i]]$model, end = model$breakdates[i + 1], frequency = model$frequency)
    })
  }
  data
}

#' @rdname get_data
#' @export
get_data.piecereg <- function(model, ...) {
  get_data(model$model, ...)
}

#' Extract Full Transformed Exogoneous matrix
#'
#' @param model the model
#' @param ... other unused parameters
#' @importFrom stats model.matrix na.pass
#' @export
full_exogeneous_matrix <- function(model, ...){
  new_diff <- new.env()
  new_diff$diff <- diff_na
  new_model <- eval(model$call, envir = new_diff)
  res <- model.matrix(new_model,
                      na.action = na.pass)
  data <- eval(model$call$data)
  if (is.ts(data))
    res <- ts(res, start = start(data), frequency = frequency(data))

  res
}

diff_na <- function(x, lag = 1, differences = 1, ...){
  res <- diff(x, lag = differences, differences = differences, ...)
  if (is.matrix(x)){
    res <- rbind(matrix(NA, ncol = ncol(res),
                        nrow = nrow(x) - nrow(res)),
                 res)
  } else {
    res <- c(rep(NA, length(x) - length(res)), res)
  }
  if (is.ts(x))
    res <- ts(res, start = start(x), frequency = frequency(x))
  res
}
