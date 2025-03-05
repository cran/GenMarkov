#' Calculate Initial Values
#'
#' Initial Values p. 387 de Berchtold (2001) - (This function is adapted from march::march.mtd.construct())
#'
#' @param f a numerical array
#' @param split_result a numerical vector, to split indexes of array
#'
#' @return a numerical vector
#' @keywords internal
#' @noRd
#'
#' @examples
#' data(stockreturns)
#' s <- cbind(stockreturns$sp500, stockreturns$djia)
#' f <- ProbMatrixF(s)
#' CalculateInitialValues(f, c(1, 2)) ## for first equation
#' CalculateInitialValues(f, c(3, 4)) ## for second equation
CalculateInitialValues <- function(f, split_result = split_result) {
  u <- apply(f, 3, FUN = function(cg) {
    tc <- sum(cg)
    sr <- rowSums(cg)
    sc <- colSums(cg)

    num_in <- cg * log2(sc) + t(t(cg) * log2(sr)) - cg * log2(cg) - cg * log2(tc)

    num_in[is.nan(num_in)] <- 0

    num <- sum(num_in)

    den_in <- sc * log2(sc) - sc * log2(tc)

    den_in[is.nan(den_in)] <- 0

    den <- sum(den_in)

    u <- num / den
    u[is.infinite(u)] <- 0

    return(u)
  })

  lambda <- lapply(split_result, FUN = function(x) u[x] / sum(u[x]))

  return(lambda = unlist(lambda))
}
