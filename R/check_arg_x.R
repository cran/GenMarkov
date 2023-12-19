#' Error handling for mmcx() function
#'
#' @param s matrix of categorical data sequences
#' @param x matrix of covariates
#' @param initial numerical vector of initial values.
#'
#' @return the function returns an error if the formats of the arguments are not correct
#' @keywords internal
#' @noRd
#'
#' @examples
#' data(stockreturns)
#' s <- cbind(stockreturns$sp500, stockreturns$djia)
#' x <- stockreturns$spread_1
#' check_arg_x(s, x, initial = c(1, 1))
check_arg_x <- function(s, x, initial) {
  if (!is.numeric(s)) {
    stop("Argument 's' is not numerical.")
  } else if (any(s %% 1 != 0)) {
    stop("Argument 's' contains non-discrete elements.")
  } else if (nrow(s) != length(x)) {
    stop("Argument 's' and 'x' do not have same size.")
  } else if (length(initial) != ncol(s)) {
    stop("Initial values inserted do not
         have the same size as the models parameters.")
  }
}
