#' @title Stock returns data
#'
#' @description Data from 5-week-day daily stock returns (rt = 100 x log(Pt/Pt-1), where Pt is the adjusted close price) of two indexes, S&P500 and DJIA, from November 11th 2011 to September 1st 2021.
#' The dataset also includes the interest rate spread, the 10-Year Treasury Constant Maturity Minus 3-Month Treasury Constant Maturity.
#' The data was retrieved from FRED.
#'
#' @format
#'
#' A tibble with 2,581 rows and 4 columns:
#'
#' \describe{
#'   \item{date}{yyyy-mm-dd of the closing price}
#'   \item{sp500}{S&P500 returns' quantiles}
#'   \item{djia}{DJIA returns' quantiles}
#'   \item{spread_1}{Lagged 10-Year Treasury Constant Maturity Minus 3-Month Treasury Constant Maturity}
#'   \item{returns_sp500}{S&P500 returns}
#'   \item{djia}{DJIA returns}
#' }
#' @source
#' \url{https://fred.stlouisfed.org/series/SP500}
#'
#' \url{https://fred.stlouisfed.org/series/DJIA}
#'
#' \url{https://fred.stlouisfed.org/series/T10Y3M}
#'
"stockreturns"
