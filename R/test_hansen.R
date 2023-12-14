#' Hansen Test
#'
#' Performs Hansen test
#'
#' @param x `lm` object.
#' @param var variables used for the joint test. By default all the variable are used.
#' @param sigma `logical` indicating if the joint test should include the variance
#'
#' @details
#' Perform Hansen test, which indicates if the variance of a model, a global model and the coefficients of the variable within this model are likely to be unstable over time.
#'
#' HO: the coefficient/model is stable over time.
#'
#' @references Bruce E Hansen "Testing for parameter instability in linear models". Journal of policy Modeling (1992)
#'
#'
#' @export

hansen_test <- function(x, var, sigma = FALSE) {
  if (!inherits(x, "lm")) {
    stop('x must be a "lm" object')
  }
  e_t <- residuals(x)
  intercept <- length(grep("Intercept", names(coef(x)))) > 0
  if (intercept) {
    x_reg <- cbind(1, x$model[, -1])
  } else {
    x_reg <- x$model[, -1]
  }
  if (missing(var)) {
    var <- 1:length(coef(x))
  }
  # On ajoute la constante
  if (sigma) {
    var <- c(var, length(coef(x)) + 1)
  }
  var <- unique(var) # to remove duplicated

  sigma2 <- mean(e_t^2, na.rm = TRUE)
  f <- cbind(apply(x_reg, 2, `*`, e_t), e_t^2 - sigma2) # i en colonne, t en ligne
  colnames(f) <- c(names(coef(x)), "sigma2")

  S <- apply(f, 2, cumsum)
  V_i <- colSums(f^2)
  L <- colMeans(S^2) / V_i
  L_c <- tryCatch(
    {
      # Todo si une seule var selectionnee ne marche pas
      f <- f[, var]
      S <- S[, var, drop = FALSE]
      V <- lapply(1:nrow(f), function(i) matrix(f[i, ], ncol = 1) %*% matrix(f[i, ], nrow = 1))
      V <- Reduce(`+`, V)
      L_c <- lapply(1:nrow(S), function(i) matrix(S[i, ], nrow = 1) %*% solve(V, matrix(S[i, ], ncol = 1)))
      L_c <- Reduce(`+`, L_c)
      c(L_c / nrow(S))
    },
    error = function(e) {
      warning("Test joint impossible: supprimer les indicatrices")
      NA
    }
  )
  tests <- c(L, L_c)
  res <- list(
    L = L,
    L_c = L_c,
    selected_var = var
  )
  class(res) <- "hansen_test"
  return(res)
}

#' @export
hansen.test <- function (...) {
  .Deprecated("hansen_test")
  hansen_test(...)
}

#' @export
print.hansen_test <- function(x, a = c(5, 1, 2.5, 7.5, 10, 20), digits = 4, ...) {
  cat("\n")
  cat("Variable                 ", "L       ", "Stat    ", "Conclusion", "\n")
  cat("______________________________________________________________", "\n")
  k <- length(x$L)
  b <- paste0(a, "%")
  b <- match.arg(b[1], choices = c("1%", "2.5%", "5%", "7.5%", "10%", "20%"))
  sigma <- length(grep("sigma2", names(x$L))) > 0
  if (sigma) {
    u <- length(x$selected_var) + 1
  } else {
    u <- length(x$selected_var)
  }
  rejet <- c(
    x$L >= hansen_table[1, b],
    x$L_c >= hansen_table[u, b]
  )
  if (sigma) {
    names <- rbind(as.matrix(names(x$L)[-length(names(x$L))]), "Variance", "Joint Lc")
  } else {
    names <- rbind(as.matrix(names(x$L)), "Joint Lc")
  }
  datnames <- format(names, digits = 4)
  stat <- c(rep(hansen_table[1, b], times = k), hansen_table[u + 1, b])
  l <- format(c(x$L, x$L_c), digits = 4)
  rejet <- format(rejet, digits = 4)
  for (j in 1:nrow(names)) {
    cat(datnames[j], " ", l[j], " ", stat[j], " ", rejet[j], " ", "\n")
  }
  cat("\n")
  cat("\n")
  cat(sprintf("Lecture: True means reject H0 at level %s", b), "\n")
}

#' Detect Fixed or Moving Coefficients
#'
#' Functions to test if any coefficient is fixed or moving according to the Hansen test ([hansen_test()])
#'
#' @inheritParams hansen_test
#' @param a level
#' @param intercept boolean indicating if the intercept should be consider as a moving coefficient when at least one other variable is moving.
#'
#' @return `NULL` if no variable selected, otherwise the order of the variables.
#' @export

moving_coefficients <- function(x, a = c(5, 1, 2.5, 7.5, 10, 20), sigma = FALSE, intercept = TRUE) {
  b <- paste0(a, "%")
  b <- match.arg(b[1], choices = c("1%", "2.5%", "5%", "7.5%", "10%", "20%"))
  test <- hansen_test(x, sigma = sigma)
  uni_tests <- test$L
  if (!sigma)
    uni_tests <- uni_tests[-grep("sigma2", names(uni_tests))]

  uni_var = which(uni_tests >= hansen_table[1, b])

  if(length(uni_var) == 0) {
    if (test$L_c >= hansen_table[length(test$selected_var), b])
      warning("Result not conform with joint test")
    return(NULL)
  }

  has_unique_intercept <- length(grep("Intercept", names(coef(x)))) == 1

  if(length(uni_var) == 1) {
    if (intercept && has_unique_intercept)
      uni_var <- unique(1, uni_var)
    return(uni_var)
  }
  joint_test <- hansen_test(x, sigma = sigma, var = uni_var)
  if (is.na(joint_test$L_c)) {
    warning("Joint test impossible, check dummies")
  } else if (joint_test$L_c < hansen_table[length(uni_var), b]) {
    warning("Result not conform with joint test")
  }
  if (intercept && has_unique_intercept)
    uni_var <- unique(1, uni_var)
  return(uni_var)
}

#' @name moving_coefficients
#' @export

fixed_coefficients <- function(x, a = c(5, 1, 2.5, 7.5, 10, 20), sigma = FALSE, intercept = TRUE) {
  b <- paste0(a, "%")
  b <- match.arg(b[1], choices = c("1%", "2.5%", "5%", "7.5%", "10%", "20%"))
  test <- hansen_test(x, sigma = sigma)
  uni_tests <- test$L
  if (!sigma)
    uni_tests <- uni_tests[-grep("sigma2", names(uni_tests))]

  uni_var = which(uni_tests < hansen_table[1, b])

  if(length(uni_var) == 0) {
    if (test$L_c < hansen_table[length(test$selected_var), b])
      warning("Result not conform with joint test")
    return(NULL)
  }

  has_unique_intercept <- length(grep("Intercept", names(coef(x)))) == 1


  if(length(uni_var) == 1) {
    if (intercept && has_unique_intercept && uni_var == 1)
      uni_var <- NULL
    return(uni_var)
  }
  joint_test <- hansen_test(x, sigma = sigma, var = uni_var)
  if (is.na(joint_test$L_c)) {
    warning("Joint test impossible, check dummies")
  } else if (joint_test$L_c >= hansen_table[length(uni_var), b]) {
    warning("Result not conform with joint test")
  }
  if (intercept && has_unique_intercept)
    uni_var <- uni_var[uni_var != 1]
  uni_var
}


# usethis::use_data(hansen_table)
# hansen_table <- structure(
#   list(
#     `Degrees of freedom (m + 1)` = 1:20,
#     `1%` = c(
#       0.748, 1.07, 1.35, 1.6, 1.88, 2.12, 2.35, 2.59, 2.82, 3.05, 3.27, 3.51,
#       3.69, 3.9, 4.07, 4.3, 4.51, 4.73, 4.92, 5.13
#     ),
#     `2.5%` = c(
#       0.593, 0.898, 1.16, 1.39, 1.63, 1.89, 2.1, 2.33, 2.55, 2.76, 2.99, 3.18,
#       3.39, 3.6, 3.81, 4.01, 4.21, 4.4, 4.6, 4.79
#     ),
#     `5%` = c(
#       0.47, 0.749, 1.01, 1.24, 1.47, 1.68, 1.9, 2.11, 2.32, 2.54, 2.75, 2.96,
#       3.15, 3.34, 3.54, 3.75, 3.95, 4.14, 4.33, 4.52
#     ),
#     `7.5%` = c(
#       0.398, 0.67, 0.913, 1.14, 1.36, 1.58, 1.78, 1.99, 2.19, 2.4, 2.6, 2.81,
#       3, 3.19, 3.38, 3.58, 3.77, 3.96, 4.16, 4.36
#     ),
#     `10%` = c(
#       0.353, 0.61, 0.846, 1.07, 1.28, 1.49, 1.69, 1.89, 2.1, 2.29, 2.49, 2.69,
#       2.89, 3.08, 3.26, 3.46, 3.64, 3.83, 4.03, 4.22
#     ),
#     `20%` = c(
#       0.243, 0.469, 0.679, 0.883, 1.08, 1.28, 1.46, 1.66, 1.85, 2.03, 2.22,
#       2.41, 2.59, 2.77, 2.95, 3.14, 3.32, 3.5, 3.69, 3.86
#     )
#   ),
#   class = "data.frame", row.names = c(NA, -20L)
# )
