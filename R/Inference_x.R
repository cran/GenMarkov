#' Inference for Generalized Multivariate Markov Chains mmcx()
#'
#' @param hess numerical matrix
#' @param lambda numerical vector
#'
#' @return list with standard error, z-statistics and pvalues associated with each estimate
#' @keywords internal
#' @noRd
#'
#' @examples
#' lambda_exmp <- c(0.5, 0.5)
#' hess_exmp <- matrix(c(0.5, 0.25, 0.25, 0.5), ncol = 2, byrow = TRUE)
#' Inference_x(hess_exmp, lambda_exmp)
Inference_x <- function(hess, lambda) {
  if (matrixcalc::is.singular.matrix(hess, tol = 1e-05) == FALSE) {
    hessinv <- solve(-hess)

    var <- diag(hessinv)

    ifelse(any(var < 0),
      return(l = list(
        warning = 1,
        se = ".",
        zstat = ".",
        pvalue = "."
      )),
      return(l = list(
        warning = 0,
        se = sqrt(var),
        zstat = lambda / sqrt(var),
        pvalue = 2 * (1 - stats::pnorm(abs(lambda / sqrt(var))))
      ))
    )
  } else {
    return(l = list(
      warning = 1,
      se = ".",
      zstat = ".",
      pvalue = "."
    ))
  }
}
