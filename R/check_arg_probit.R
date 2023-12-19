#' Error handling for multi.mtd_probit()
#'
#' @param s matrix of categorical data sequences
#' @param initial numerical vector of initial values
#' @param nummethod Numerical maximisation method, currently either "NR" (for Newton-Raphson), "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), "BFGSR" (for the BFGS algorithm implemented in R), "BHHH" (for Berndt-Hall-Hall-Hausman), "SANN" (for Simulated ANNealing), "CG" (for Conjugate Gradients), or "NM" (for Nelder-Mead). Lower-case letters (such as "nr" for Newton-Raphson) are allowed. The default method is "BFGS". For more details see maxLik.
#'
#' @return the function returns an error if the formats of the arguments are not correct
#' @keywords internal
#' @noRd
#' @examples
#' data(stockreturns)
#' s <- cbind(stockreturns$sp500, stockreturns$djia)
#' check_arg_probit(s, initial = c(1, 1, 1), nummethod = "bfgs")
check_arg_probit <- function(s, initial, nummethod) {
  if (!is.numeric(s)) {
    stop("Argument 's' is not numerical.")
  } else if (any(s %% 1 != 0)) {
    stop("Argument 's' contains non-discrete elements.")
  } else if (length(initial) != ncol(s) + 1) {
    stop("Initial values inserted do not
         have the same size as the models parameters.")
  } else if (!(tolower(nummethod) %in% c(
    "nr", "bfgsr", "bfgs", "bhhh",
    "sann", "cg", "nm"
  ))) {
    stop("Maximisation method not allowed.")
  }
}
