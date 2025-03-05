#' @title Non-homogeneous Multivariate Markov Chains
#'
#' @description Estimates Multivariate Markov Chains that depend on a exogeneous variables.
#' The model is based on the Mixture Transition Distribution model, and considers non-homogeneous Markov Chains, instead of homogeneous Markov Chains as in Raftery (1985).
#'
#' @param y matrix of categorical data sequences
#' @param x matrix of covariates
#' @param initial numerical vector of initial values.
#' @param ... additional arguments to be passed down to \code{\link[alabama:auglag]{auglag()}}.
#'
#' @return The function returns a list with the parameter estimates, standard-errors, z-statistics, p-values and the value of the log-likelihood function, for each equation.
#' @export
#'
#' @author Carolina Vasconcelos and Bruno Damásio
#'
#' @seealso Optimization is done through \code{\link[alabama:auglag]{auglag()}}.
#'
#' @examples
#' data(stockreturns)
#' s <- cbind(stockreturns$sp500, stockreturns$djia)
#' x <- stockreturns$spread_1
#' mmcx(s, x, initial = c(1, 1))
#'
#' @references
#' Raftery, A. E. (1985). A Model for High-Order Markov Chains. Journal of the Royal Statistical Society. Series B (Methodological), 47(3), 528-539. \url{http://www.jstor.org/stable/2345788}
#'
#' Ching, W. K., E. S. Fung, and M. K. Ng (2002). A multivariate Markov chain model for categorical data sequences and its applications in demand predictions. IMA Journal of Management Mathematics, 13(3), 187-199. \doi{10.1093/imaman/13.3.187}
mmcx <- function(y, x, initial, ...) {
  check_arg_x(s = y, x = x, initial = initial)


  m1 <- max(y)
  result <- list(NA)

  q <- ProbValuesXDependent(s = y, x = x) # Calculate P(St | St-1, x) for all equations

  # Define indices to select appropriate values, for each equation
  index_vector <- 1:(ncol(y) * ncol(y))

  split_result <- lapply(
    split(1:length(index_vector), rep(1:ncol(y), each = ncol(y))),
    function(indices) index_vector[indices]
  )

  j <- 1
  for (i in split_result) {
    LogLikelihood_x <- function(lambda, qi = q[, i]) {
      qi_lambda <- qi %*% lambda

      qi_lambda[qi_lambda < 0] <- 1

      ll <- sum(-log(qi_lambda))

      ll
    }

    # Impose restrictions
    he <- function(lambda) {
      h <- rep(NA, 1)
      h[1] <- sum(lambda) - 1
      h
    }

    hi <- function(lambda) {
      w <- rep(NA, 1)
      for (i in 1:length(lambda)) {
        w[i] <- lambda[i]
      }
      w
    }

    # Jacobians of restrictions to improve optimization in auglag()
    he.jac <- function(lambda) {
      j <- matrix(NA, 1, length(lambda))
      j[1, ] <- rep(1, length(lambda))
      j
    }

    hi.jac <- function(lambda) {
      j <- diag(1, length(lambda), length(lambda))
      j
    }

    # Optimization through auglag() function
    opt <-
      alabama::auglag(
        par = initial,
        fn = LogLikelihood_x,
        hin = hi,
        heq = he,
        heq.jac = he.jac,
        hin.jac = hi.jac,
        control.outer = list(trace = FALSE),
        ...
      )

    # Only performs inference if the model reaches convergence
    if (any(opt$ineq < 0) | any(round(opt$par, 1) == 1.0)) {
      ll <- "."
      hessian <- "."
      lambdahat <- "."
      inf <- list(warning = 2)
    } else {
      hessian <- -opt$hessian
      ll <- -opt$value
      lambdahat <- opt$par

      inf <- Inference_x(hess = hessian, lambda = lambdahat)
    }

    output.table(lambdahat, inf$se, inf$zstat, inf$pvalue, ll, j, flag = inf$warning)

    table <- data.frame(cbind(lambdahat, inf$se, inf$zstat, inf$pvalue), row.names = NULL)
    colnames(table) <- c("Estimates", "Std. Error", "z value", "p-value")

    result[[j]] <- list(table, ll)
    names(result[[j]]) <- c("Estimation", "Log-likelihood")

    j <- j + 1
  }
  names(result) <- paste("Equation", 1:ncol(y), sep = " ")
  invisible(unlist(result, recursive = FALSE))
}
