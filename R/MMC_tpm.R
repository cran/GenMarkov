#' @title Transition Probability Matrices
#'
#' @description
#' This functions allows to obtain the transition probability matrices for a specific value of x, considering the estimates obtained from \code{\link[GenMarkov:mmcx]{mmcx()}}.
#'
#' @param s numerical matrix with categorical data sequences
#' @param x exogeneous variable
#' @param value fixed value of x
#' @param result result from the function \code{\link[GenMarkov:mmcx]{mmcx()}}
#'
#' @return The function returns a numerical array with the probability transition matrices for each equation
#' @export
#'
#' @author Carolina Vasconcelos and Bruno Damásio
#'
#' @examples
#' data(stockreturns)
#' s <- cbind(stockreturns$sp500, stockreturns$djia)
#' x <- stockreturns$spread_1
#' res <- mmcx(s, x, initial = c(1, 1))
#' tpm <- MMC_tpm(s, x, value = max(x), result = res)
MMC_tpm <- function(s, x, value = max(x), result) {
  m1 <- max(s)
  n <- nrow(s)

  ## Extract estimates from the result
  lambda <- unlist(lapply(
    Filter(function(x) is.data.frame(x), result),
    function(df) as.numeric(df[[1]])
  ))

  # Create matrix with dummies for each state
  dummies_list <-
    apply(s, 2, function(x) {
      fastDummies::dummy_cols(x, remove_selected_columns = TRUE)
    })
  dummies <- matrix(unlist(dummies_list),
    ncol = m1 * ncol(s),
    nrow = nrow(s)
  )

  # Create all possible combinations of column indices
  combinations <- expand.grid(1:ncol(s), 1:ncol(dummies))
  # Order by the first variable
  combinations <- combinations[order(combinations$Var1), ]

  # Extract columns from S and S_L based on the combinations
  combined_list <- lapply(1:nrow(combinations), function(i, x) {
    cbind(s[, combinations[i, 1]], x, dummies[, combinations[i, 2]])
  }, x = x)

  if (m1 == 2) {
    tpm_x <- sapply(combined_list, function(data) {
      # Define dependent variable
      y <- ifelse(data[, 1] == 1, 1, 0)

      # Define lagged St
      s_l <- Hmisc::Lag(data[, 3])

      # Estimate logistic regression
      res <- stats::glm(y[s_l == 1] ~ data[, "x"][s_l == 1], family = stats::binomial(link = "logit"))

      # Extract coefficients
      estim <- stats::coefficients(res)

      num <- exp(sum(estim * c(1, value)))
      den <- exp(sum(estim * c(1, value)))
      c(1 / (1 + num), den / (1 + num))
    })

    index <- split(1:ncol(tpm_x), combinations$Var1)

    tpm_final <- lapply(index, function(index) {
      l <- rep(lambda, each = ncol(s))[index]
      apply(tpm_x[, index], 1, function(x) rowSums(matrix(x * l, ncol = ncol(s))))
    })
  } else {
    # If m1>2, we have more than two states, we will resort to a multinomial logistic regression
    tpm_x <- sapply(combined_list, function(data) {
      # Define dependent variable
      y <- factor(data[, 1], levels = 1:max(data[, 1]))

      # Define lagged St
      s_l <- Hmisc::Lag(data[, 3])

      # Estimate multinomial logistic regression
      res <- nnet::multinom(y[s_l == 1] ~ data[, "x"][s_l == 1], trace = FALSE)

      # Extract coefficients
      estim <- stats::coefficients(res)

      # Calculate row of each Transition Probability Matrix
      num <- sum(apply(estim, 1, function(x) exp(sum(x * c(1, value)))))
      den <- apply(estim, 1, function(x) exp(sum(x * c(1, value))))
      c(1 / (1 + num), den / (1 + num))
    })

    # Change to friendlier format
    index <- split(1:ncol(tpm_x), rep(1:(ncol(s) * ncol(s)), each = m1))

    tpm_final <- lapply(index, function(index) {
      apply(tpm_x[, index], 1, function(x) t(x))
    })

    # Calculate TPM for fixed value of x
    tpm_final <- sapply(1:(ncol(s) * ncol(s)), function(x) {
      tpm_final[[x]] * lambda[x]
    }, simplify = "array")

    index <- split(1:(ncol(s) * ncol(s)), rep(1:ncol(s), each = ncol(s)))

    tpm_final <- array(data = as.numeric(sapply(index, function(x) {
      apply(tpm_final[, , x], c(1, 2), sum, simplify = FALSE)
    })), dim = c(m1, m1, ncol(s)))
  }

  return(tpm_final)
}
