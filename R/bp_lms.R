#' Pie regression
#'
#' Computes as many regressions as breakup dates within a global model
#'
#' @param x `lm` object. It is the global regression model.
#' @param left `logical`. By default set to `TRUE`, i.e. the breakdate is the end date of each submodel
#' @param break_dates optional, to indicate the breakup dates if they are known.
#' @param tvlm By default set to `FALSE`. Indicates which model will be run on each sub data. `FALSE` means a [lm] will be run.
#' @param ... other arguments passed to [tvReg::tvLM()].
#'
#'
#' @return Returns an element of class `bp_lm`. It is a list containing the following elements:
#' \item{model}{all computed models, each of class `lm` or `tvlm` according to the parameter specified above}
#' \item{start}{start date of the time serie}
#' \item{end}{end date of the time serie}
#' \item{frequency}{frequency of the time serie}
#' \item{breakdates}{a list of the breakup dates}
#' \item{left}{same as the parameter specified above}
#' \item{tvlm}{same as the parameter specified above}
#'
#' @export
#' @importFrom tvReg tvLM
#' @importFrom strucchange breakdates breakpoints
bp_lm <- function(x, left = TRUE, break_dates, tvlm = FALSE, ...) {
  if (has_intercept(x)) {
    formula <- sprintf("%s ~ -1 + .", colnames(x$model)[1])
  } else {
    formula <- sprintf("%s ~ .", colnames(x$model)[1])
  }
  .data <- get_data(x, ...)
  if (missing(break_dates)) {
    break_dates <- strucchange::breakdates(strucchange::breakpoints(as.formula(formula), data = .data))
  }
  if (all(is.na(break_dates))) {
    if (!tvlm) {
      return(x)
    } else {
      tv = tvReg::tvLM(formula = formula, data = .data, ...)
      return(tv)
    }
  } else {
    break_dates <- c(time(.data)[1] - deltat(.data) * left, break_dates, time(.data)[length(time(.data))] + deltat(.data) * !left)
    if (!tvlm) {
      model <- lapply(1:(length(break_dates) - 1), function(i) {
        lm(formula = formula, window(.data,
                                     start = (break_dates[i] + deltat(.data) * left),
                                     end = (break_dates[i + 1] - deltat(.data) * !left)
        ))
      })
    } else {
      model <- lapply(1:(length(break_dates) - 1), function(i) {
        tvLM(formula = formula,
             data = window(.data,
                           start = (break_dates[i] + deltat(.data) * left),
                           end = (break_dates[i + 1] - deltat(.data) * !left)
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
    left = left,
    tvlm = tvlm
  )
  class(res) <- "bp_lm"
  res
}

#' @importFrom tvReg tvLM
#' @importFrom strucchange breakdates breakpoints


#' @export
print.bp_lm <- function(x, ...) {
  if (is.null(x$breakdates)) {
    print("No breakdate")
  } else {
    start_date <- x$start[1] + (x$start[2] - 1) / x$frequency
    end_date <- x$end[1] + (x$end[2] - 1) / x$frequency
    v1 <- c("Start date", rep("Breakdates", time = length(x$breakdates) - 2), "End date", "Frequency", "left")
    v2 <- c(start_date, x$breakdates[-c(1, length(x$breakdates))], end_date, x$frequency, x$left)
    datav <- rbind(v2)
    colnames(datav) <- v1
    print(datav)
  }
  print(sapply(x$model, coef))
}



#' @export
plot.bp_lm <- function(x, y, plot.type = c("single", "multiple"), ...) {
  plot(cbind(fitted(x) + resid(x), fitted(x)),
       plot.type = plot.type[1],
       ...
  )
}

#' @export
fitted.bp_lm <- function(object, ...) {
  res <- unlist(lapply(object$model, fitted, ...))
  names(res) <- 1:length(res)
  ts(res,
     start = object$start,
     frequency = object$frequency
  )
}

#' @export
fitted.piece_reg <- function(object, ...) {
  ts(fitted(object$model, ...),
     start = object$start,
     frequency = object$frequency
  )
}


#' @export
residuals.bp_lm <- function(object, ...) {
  res <- unlist(lapply(object$model, residuals, ...))
  names(res) <- 1:length(res)
  ts(res,
     start = object$start,
     frequency = object$frequency
  )
}
#' @export
residuals.piece_reg <- function(object, ...) {
  ts(residuals(object$model, ...),
     start = object$start,
     frequency = object$frequency
  )
}

#' @export
coef.bp_lm <- function(object, ...) {
  if (inherits(object$model[[1]], "lm")) {
    res <- do.call(rbind, lapply(1:(length(object$breakdates) - 1), function(i) {
      ts(matrix(coef(object$model[[i]]), nrow = 1),
         start = object$breakdates[i] + 1 / object$frequency * (object$left),
         end = object$breakdates[i + 1] - 1 / object$frequency * !(object$left),
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
         start = object$breakdates[i] + 1 / object$frequency * (object$left),
         end = object$breakdates[i + 1] - 1 / object$frequency * !(object$left),
         frequency = object$frequency,
         names = colnames(coef_tvlm[[i]])
      )
    })
    res
  }
}
#'@export
coef.piece_reg <- function(object, ...) {

  coef <- coef(object$model)
  if (object$tvlm) {
    ts(coef,
       start = object$start,
       end = object$end,
       frequency = object$frequency
    )
  } else {
    end_date <- object$end[1] + (object$end[2] - 1) / object$frequency
    all_n <- gsub(paste0("(_",
                         c(object$breakdates, end_date),
                         ")", collapse = "|"),
                  "", names(coef))
    unique_names <- unique(all_n, names(coef))

    int <- ts(1,
              start = object$start,
              end = object$end,
              frequency = object$frequency)
    break_int <- break_data(int, break_dates = object$breakdates, left = object$left_breaks)
    all_coefs <- do.call(cbind, lapply(unique_names, function(name) {
      possible_var <- which(name == all_n)
      if (length(possible_var) == 1) {
        res <- coef[possible_var]
      } else {
        res <- break_int %*% matrix(coef[possible_var], ncol = 1)
      }

      ts(res,
         start = object$start,
         end = object$end,
         frequency = object$frequency
      )
    }))
    colnames(all_coefs) <- gsub("`", "", unique_names)
    all_coefs
  }
}
#'@export
breakpoints.lm <- function(obj, ...) {
  if (has_intercept(obj)) {
    formula <- sprintf("%s ~ .", colnames(obj$model)[1])
  } else {
    formula <- sprintf("%s ~ -1 + .", colnames(obj$model)[1])
  }
  strucchange::breakpoints(formula = as.formula(formula),
                           data = obj$model, ...)
}
