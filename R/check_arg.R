#' Error handling for multi.mtd()
#'
#' @param s matrix of categorical data sequences
#' @param is_constrained flag indicating whether the function will consider the usual set of constraints (usual set: TRUE, new set of constraints: FALSE).
#'
#' @return the function returns an error if the formats of the arguments are not correct
#' @keywords internal
#' @noRd
#'
#' @examples
#' data(stockreturns)
#' s <- cbind(stockreturns$sp500, stockreturns$djia)
#' check_arg(s, is_constrained = TRUE)
check_arg <- function(s, is_constrained) {
  if (!is.numeric(s)) {
    stop("Argument 's' is not numerical.")
  } else if (any(s %% 1 != 0)) {
    stop("Argument 's' contains non-discrete elements.")
  } else if (!is.logical(is_constrained)) {
    stop("Argument 'is_constrained' is not logical.")
  }
}
