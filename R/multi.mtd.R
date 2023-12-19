#' @title Estimation of Multivariate Markov Chains - MTD model
#'
#' @description
#' This function estimates the Mixture Distribution Model (Raftery (1985)) for Multivariate Markov Chains.
#' It considers Berchtold (2001) optimization algorithm for the parameters and estimates the probabilities transition matrices as proposed in Ching (2002).
#'
#'
#' @param y matrix of categorical data sequences
#' @param deltaStop value below which the optimization phases of the parameters stop
#' @param is_constrained flag indicating whether the function will consider the usual set of constraints (usual set: \emph{TRUE}, new set of constraints: \emph{FALSE}).
#' @param delta the amount of change to increase/decrease in the parameters for each iteration of the optimization algorithm.
#'
#' @return The function returns a list with the parameter estimates, standard-errors, z-statistics, p-values and the value of the log-likelihood function, for each equation.
#' @export
#'
#' @note
#' See details of the optimization procedure in \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/1467-9892.00231}{Berchtold (2001)}.
#'
#' @references
#'
#' Raftery, A. E. (1985). A Model for High-Order Markov Chains. Journal of the Royal Statistical Society. Series B (Methodological), 47(3), 528-539. \url{http://www.jstor.org/stable/2345788}
#'
#' Berchtold, A. (2001). Estimation in the Mixture Transition Distribution Model. Journal of Time Series Analysis, 22(4), 379-397.\doi{10.1111/1467-9892.00231}
#'
#' Ching, W. K., E. S. Fung, and M. K. Ng (2002). A multivariate Markov chain model for categorical data sequences and its applications in demand predictions. IMA Journal of Management Mathematics, 13(3), 187-199. \doi{10.1093/imaman/13.3.187}
#'
#' @examples
#' data(stockreturns)
#' s <- cbind(stockreturns$sp500, stockreturns$djia)
#' multi.mtd(s)
multi.mtd <-
  function(y,
           deltaStop = 0.0001,
           is_constrained = TRUE,
           delta = 0.1) {
    # Check arguments
    check_arg(y, is_constrained)

    results <- list(NA)

    # Log-likelihood (Berchtold (2001))
    LogLikelihood <- function(lambda, ni, qi) {
      ll <- sum(t(ni) %*% log(qi %*% lambda))

      return(ll)
    }

    # Numerical maximization: Algorithm from Berchtold (2001) - (This function is adapted from march::march.mtd.construct())
    OptimizeLambda <-
      function(lambda,
               delta,
               ni,
               qi,
               is_constrained,
               delta_stop) {
        delta_it <- delta

        ll <- 0
        ll <- LogLikelihood(lambda = lambda, ni = ni, qi = qi)
        pd_lambda <- 0
        pd_lambda <- PartialDerivatives(ni = ni, qi = qi, lambda = lambda)

        i_inc <- which.max(pd_lambda)
        i_dec <- which.min(pd_lambda)
        lambda_r <- lambda
        par_inc <- lambda_r[i_inc]
        par_dec <- lambda_r[i_dec]

        if (is_constrained) {
          if (par_inc == 1) {
            return(list(
              lambda = lambda,
              ll = ll,
              delta = delta
            ))
          }
          if (par_inc + delta_it > 1) {
            delta_it <- 1 - par_inc
          }
          if (par_dec == 0) {
            pd_lambda_sorted <- sort(pd_lambda, index.return = TRUE)
            i_dec <-
              pd_lambda_sorted$ix[min(which(lambda[pd_lambda_sorted$ix] > 0))]
            par_dec <- lambda_r[i_dec]
          }
          if (par_dec - delta_it < 0) {
            delta_it <- 1 - par_dec
          }
        }

        # Otimization loop
        while (TRUE) {
          if (is_constrained) {
            delta_it <- min(c(delta_it, 1 - par_inc, par_dec))
          }

          new_lambda <- lambda_r
          new_lambda[i_inc] <- par_inc + delta_it
          new_lambda[i_dec] <- par_dec - delta_it
          new_ll <- LogLikelihood(new_lambda, ni, qi)
          if (new_ll > ll) {
            return(list(
              lambda = new_lambda,
              ll = new_ll,
              delta = delta
            ))
          } else {
            if (delta_it <= delta_stop) {
              delta <- 2 * delta
              return(list(
                lambda = new_lambda,
                ll = new_ll,
                delta = delta
              ))
            }
            delta_it <- delta_it / 2
          }
        }
      }


    nr.eq <- ncol(y)

    # Define indices to select appropriate values, for each equation
    index_vector <- 1:(ncol(y) * ncol(y))

    split_result <- lapply(
      split(1:length(index_vector), rep(1:ncol(y), each = ncol(y))),
      function(indices) index_vector[indices]
    )

    f <- ProbMatrixF(y)
    q <- ProbMatrixQ(y)

    ni <- apply(f, 3, function(x) matrixcalc::vec(x))
    qi <- apply(q, 3, function(x) matrixcalc::vec(x))

    lambda <- CalculateInitialValues(f, split_result = split_result)

    j <- 1
    for (i in split_result) {
      n <- ni[, i]
      q <- qi[, i]
      l <- lambda[i]

      opt <- OptimizeLambda(l, delta, n, q, is_constrained, deltaStop)

      lambdaOptim <- opt$lambda
      ll <- opt$ll

      inf <- Inference(n, q, lambdaOptim)

      output.table(lambdaOptim, inf$se, inf$zstat, inf$pvalue, ll, j, flag = inf$flag)

      table <- data.frame(cbind(lambdaOptim, inf$se, inf$zstat, inf$pvalue), row.names = NULL)
      colnames(table) <- c("Estimates", "Std. Error", "z value", "p-value")

      results[[j]] <- list(table, ll)
      names(results[[j]]) <- c("Estimation", "Log-likelihood")

      j <- j + 1
    }
    names(results) <- paste("Equation", 1:ncol(y), sep = " ")
    results <- unlist(results, recursive = FALSE)
    invisible(results)
  }
