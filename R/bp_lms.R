#' Breakdates regression
#'
#' Computes as many regressions as breakup dates within a global model
#'
#' @param x `lm` object. It is the global regression model
#' @param data 	a data frame, list or environment containing the variables in the model
#' @param right `logical`. By default set to `TRUE`, i.e. the breakdate is the end date of each submodel
#' @param break_dates optional, to indicate the breakup dates if they are known
#' @param tvlm By default set to `FALSE`. Indicates which model will be run on each sub data. FALSE means a [lm] will be run.
#'
#' @return Returns an element of class `bp.lms`. It is a list containing the following elements :
#' \item{model}{all computed models, each of class `lm` or `tvlm` according to the parameter specified above}
#' \item{start}{start date of the time serie}
#' \item{end}{end date of the time serie}
#' \item{frequency}{frequency of the time serie}
#' \item{breakdates}{a list of the breakup dates}
#' \item{right}{same as the parameter specified above}
#' \item{tvlm}{same as the parameter specified above}
#'
#' @export

bp.lms <- function(x, data, right = TRUE, break_dates, tvlm = FALSE, ...) {
  formula <- sprintf("%s ~ .", colnames(x$model)[1])
  .data <- ts(x$model, end = end(data), frequency = frequency(data))
  if (missing(break_dates)) {
    break_dates <- breakdates(breakpoints(as.formula(formula), data = x$model))
  }
  if (all(is.na(break_dates))) {
    if (!tvlm) {
      return(x)
    } else {
      tv = tvLM(formula = formula, data = x$model, ...)
      return(tv)
    }
  } else {
    break_dates <- c(time(.data)[1] - deltat(.data) * right, break_dates, time(.data)[length(time(.data))] + deltat(.data) * !right)
    if (!tvlm) {
      model <- lapply(1:(length(break_dates) - 1), function(i) {
        lm(formula = formula, window(.data,
                                     start = (break_dates[i] + deltat(.data) * right),
                                     end = (break_dates[i + 1] - deltat(.data) * !right)
        ))
      })
    } else {
      model <- lapply(1:(length(break_dates) - 1), function(i) {
        tvLM(formula = formula,
             data = window(.data,
                           start = (break_dates[i] + deltat(.data) * right),
                           end = (break_dates[i + 1] - deltat(.data) * !right)
             ),
             ...)
      })
    }
  }
  res <- list(
    model = model,
    start = start(.data),
    end = end(.data),
    frequency = frequency(.data),
    breakdates = break_dates,
    right = right,
    tvlm = tvlm
  )
  class(res) <- "bp.lms"
  res
}


#' @export
print.bp.lms <- function(x, ...) {
  if (is.null(x$breakdates)) {
    print("No breakdate")
  } else {
    start_date <- x$start[1] + (x$start[2] - 1) / x$frequency
    end_date <- x$end[1] + (x$end[2] - 1) / x$frequency
    v1 <- c("Start date", rep("Breakdates", time = length(x$breakdates) - 2), "End date", "Frequency", "Right")
    v2 <- c(start_date, x$breakdates[-c(1, length(x$breakdates))], end_date, x$frequency, x$right)
    datav <- kable(rbind(v2), col.names = v1, row.names = FALSE)
    print(datav)
  }
  print(kable(sapply(x$model, coef)))
}



#' @export
plot.bp.lms <- function(x, y, plot.type = c("single", "multiple"), ...) {
  plot(cbind(fitted(x) + resid(x), fitted(x)),
       plot.type = plot.type[1],
       ...
  )
}

#' @export
fitted.bp.lms <- function(object, ...) {
  res <- unlist(lapply(object$model, fitted, ...))
  names(res) <- 1:length(res)
  ts(res,
     start = object$start,
     frequency = object$frequency
  )
}

#' @export
residuals.bp.lms <- function(object, ...) {
  res <- unlist(lapply(object$model, residuals, ...))
  names(res) <- 1:length(res)
  ts(res,
     start = object$start,
     frequency = object$frequency
  )
}

#' @export
coef.bp.lms <- function(object, ...) {
  if (inherits(object$model[[1]], "lm")) {
    res <- do.call(rbind, lapply(1:(length(object$breakdates) - 1), function(i) {
      ts(matrix(coef(object$model[[i]]), nrow = 1),
         start = object$breakdates[i] + 1 / object$frequency * (object$right),
         end = object$breakdates[i + 1] - 1 / object$frequency * !(object$right),
         frequency = object$frequency,
         names = rownames(sapply(object$model, coef))
      )
    }))
    ts(res,
       start = object$start,
       frequency = object$frequency
    )
  } else {
    coef_tvlm = sapply(object$model, coef)
    res <- sapply(seq_along(coef_tvlm), function(i){
      ts(coef_tvlm[[i]],
         start = object$breakdates[i] + 1 / object$frequency * (object$right),
         end = object$breakdates[i + 1] - 1 / object$frequency * !(object$right),
         frequency = object$frequency,
         names = colnames(coef_tvlm[[i]])
      )
    })
    res
  }
}
