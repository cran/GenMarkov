#' Output table
#'
#' @param estimates numerical vector of estimates
#' @param se numerical vector of standard error of estimates
#' @param zstat numerical vector of z-statistics of estimates
#' @param pvalue numerical vector of pvalue
#' @param ll numerical value of loglikelihood value
#' @param eq numerical value of number of equation
#' @param flag numerical value of warning in optimization
#'
#' @return prints the results table
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
#' inf <- Inference(ni_exmp, qi_exmp, lambda_exmp)
#' output.table(lambda_exmp, inf$se, inf$zstat, inf$pvalue, ll = 2, eq = 1)
output.table <- function(estimates, se, zstat, pvalue, ll, eq, flag = 0) {
  stars <- rep("", length(pvalue))

  if (!is.character(se)) {
    stars[pvalue <= 0.01] <- "***"
    stars[pvalue > 0.01 & pvalue <= 0.05] <- "**"
    stars[pvalue > 0.05 & pvalue <= 0.1] <- "*"

    se <- formatC(se, digits = 6, format = "f")
    zstat <- formatC(zstat, digits = 3, format = "f")
    pvalue <- formatC(pvalue, digits = 3, format = "f")
    stars <- format(stars, justify = "left")
  } else {
    se <- rep(".", length(pvalue))
    zstat <- rep(".", length(pvalue))
    pvalue <- rep(".", length(pvalue))
  }

  estimates <- formatC(estimates, digits = 6, format = "f")
  results <- data.frame(cbind(estimates, se, zstat, pvalue, stars), row.names = NULL)

  # Print output table and warnings messages (if necessary)
  if (flag == 2) {
    cat("--------------------------------------------\n")
    cat("Equation", eq, "\n")
    cat("Algorithm did not reach a solution with the constraints imposed.")
    cat("\n")
    cat("--------------------------------------------\n")
  } else if (flag == 1) {
    cat("--------------------------------------------\n")
    cat("Equation", eq, "\n")
    print(results)
    cat("\n")
    cat("Log-Likelihood:", ll, "\n")
    cat("--------------------------------------------\n")
    cat("Hessian is singular, cannot compute standard errors.")
  } else {
    namcol <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
    colnames(results) <- namcol
    cat("--------------------------------------------\n")
    cat("Equation", eq, "\n")
    print(results)
    cat("\n")
    cat("Log-Likelihood:", ll, "\n")
    cat("--------------------------------------------\n")
  }
}
