#' First Partial derivatives of MTD
#'
#' Calculate first partial derivatives for numerical maximization - (This function is adapted from march::march.mtd.construct())
#'
#' @param ni numerical matrix of number of transitions between states and data sequences
#' @param qi numerical matrix of transitions probabilities between states and data sequences
#' @param lambda numerical vector
#'
#' @return numerical vector with partial derivates values of each lambda
#' @keywords internal
#' @noRd
#'
#' @examples
#' data(stockreturns)
#' s <- cbind(stockreturns$sp500, stockreturns$djia)
#' f <- ProbMatrixF(s)
#' q <- ProbMatrixQ(s)
#' ni_exmp <- apply(f[, , 1:2], 3, function(x) matrixcalc::vec(x))
#' qi_exmp <- apply(q[, , 1:2], 3, function(x) matrixcalc::vec(x))
#' lambda_exmp <- c(0.5, 0.5)
#'
#' PartialDerivatives(ni_exmp, qi_exmp, lambda_exmp)
PartialDerivatives <- function(ni, qi, lambda) {
  den <- apply(qi, 1, FUN = function(q) {
    q %*% lambda
  })

  pd_lambda <- colSums(ni * (qi / den))
  pd_lambda[is.infinite(pd_lambda)] <- 0

  return(pd_lambda)
}
