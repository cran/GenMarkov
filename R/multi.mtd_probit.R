#' @title Estimation of Multivariate Markov Chains: MTD - Probit Model
#'
#' @description Estimation of Multivariate Markov Chains through the proposed model by Nicolau (2014).
#' This model presents two attractive features: it is completely free of constraints, thereby facilitating the estimation procedure,
#' and it is more precise at estimating the transition probabilities of a multivariate or higher-order Markov chain than the Raftery's MTD model.
#'
#' @param y matrix of categorical data sequences
#' @param initial numerical vector of initial values
#' @param nummethod Numerical maximisation method, currently either "NR" (for Newton-Raphson), "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), "BFGSR" (for the BFGS algorithm implemented in R), "BHHH" (for Berndt-Hall-Hall-Hausman), "SANN" (for Simulated ANNealing), "CG" (for Conjugate Gradients), or "NM" (for Nelder-Mead). Lower-case letters (such as "nr" for Newton-Raphson) are allowed. The default method is "BFGS". For more details see \code{\link[maxLik:maxLik]{maxLik()}}.
#'
#' @return The function returns a list with the parameter estimates, standard-errors, z-statistics, p-values and the value of the log-likelihood function, for each equation.
#' @export
#'
#' @author Carolina Vasconcelos and Bruno Damásio
#'
#' @examples
#' data(stockreturns)
#' s <- cbind(stockreturns$sp500, stockreturns$djia)
#' multi.mtd_probit(s, initial = c(1, 1, 1), nummethod = "bfgs")
#'
#' @references
#' Nicolau, J. (2014). A new model for multivariate markov chains. Scandinavian Journal of Statistics, 41(4), 1124-1135.\doi{10.1111/sjos.12087}
multi.mtd_probit <- function(y, initial, nummethod = "bfgs") {
  # Check arguments
  check_arg_probit(y, initial, nummethod)

  # Define number of equations to estimate
  nr.eq <- ncol(y)

  # Calculate frequency transition matrices and vectorize
  ni <- apply(ProbMatrixF(y), 3, function(x) matrixcalc::vec(x))

  # Calculate probability transition matrices and vectorize
  qi <- apply(ProbMatrixQ(y), 3, function(x) matrixcalc::vec(x))

  # Define indices to select appropriate columns of ni and qi, for each equation
  index_vector <- 1:ncol(ni)

  split_result <- lapply(
    split(1:length(index_vector), rep(1:nr.eq, each = nr.eq)),
    function(indices) index_vector[indices]
  )


  # Obtain etas for each equation
  results <- lapply(split_result, FUN = function(i) {
    neq <- ni[, i]
    qeq <- qi[, i]

    # Define log-likelihood
    LogLikelihood_p <- function(eta, n = neq, qn = qeq) {
      ll <- 0

      denom <- unlist(lapply(split_result, function(i) {
        stats::pnorm(sum(cbind(rep(1, nrow(qi[, i])), qi[, i]) %*% eta))
      }))

      eta_mat <- log(stats::pnorm(cbind(rep(1, nrow(qn)), qn) %*% eta) / sum(denom))

      ll <- sum(t(n) %*% eta_mat)

      return(ll)
    }

    # Maximize log-likelihood
    otim <- maxLik::maxLik(LogLikelihood_p, start = initial, method = nummethod)

    res <- summary(otim)

    return(res)
  })

  # Return results
  names(results) <- paste("Equation", 1:nr.eq)

  return(results)
}
