#' Inference function to compute standard errors, t-values and p-values of MTD model (of the multi.mtd.R function)
#'
#' @param ni numerical matrix of transitions frequencies
#' @param qi numerical matrix of transitions probabilities
#' @param lambda numerical vector
#'
#' @return list with standard error, z-statistics and pvalues associated with each estimate
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
#' Inference(ni_exmp, qi_exmp, lambda_exmp)
Inference <- function(ni, qi, lambda) {
  den <- apply(qi, 1, FUN = function(q) {
    (q %*% lambda)^2
  })

  num <- apply(qi, 2, FUN = function(qi_c) {
    (qi_c * qi)
  }, simplify = FALSE)

  hess_list <- apply(array(unlist(num), dim = c(nrow(ni), length(lambda), length(lambda))), 3,
    FUN = function(num) {
      colSums(ni * (-(num / den)))
    }, simplify = FALSE
  )

  hess <- matrix(unlist(hess_list), nrow = length(lambda), ncol = length(lambda), byrow = TRUE)

  if (matrixcalc::is.singular.matrix(hess)) {
    flag <- 1
    hessinv <- matrix(NA, nrow = length(lambda), ncol = length(lambda), byrow = TRUE)
    var <- diag(hessinv)
    se <- sqrt(var)
    zstat <- se
    pvalue <- se
  } else {
    flag <- 0
    hessinv <- solve(-hess)
    var <- diag(hessinv)
    se <- sqrt(var)
    zstat <- lambda / se
    pvalue <- 2 * (1 - stats::pnorm(abs(zstat)))
  }


  return(l = list(
    flag = flag,
    se = se,
    zstat = zstat,
    pvalue = pvalue
  ))
}
