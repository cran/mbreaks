#' World Health Organization TB data
#'
#' Data set from the Garcia and Perron study's of ex-post exchange rate.
#'
#' @format ## `real`
#' A data frame with 103 rows and 1 column:
#' \describe{
#'   \item{rate}{Real exchange rate}
#' }
#' @source Garcia, R. and Perron, P., 1996. "An analysis of the real interest rate under regime shifts." Review of Economics and Statistics 111–125.
"real"

#' New Keynesian Phillips curve data
#'
#' Data set from inflation and other macroeconomic variables
#'
#' @format ## `nkpc`
#' A data frame with 151 rows and 12 columns:
#' \describe{
#'   \item{year}{Current period year}
#'   \item{quarter}{Quarter in current period year}
#'   \item{inf}{Inflation rate}
#'   \item{inflag}{Inflation rate in previous period}
#'   \item{inffut}{Expected inflation rate, taken as value of inflation rate of next period}
#'   \item{ygap}{Productivity output gap}
#'   \item{lbs}{}
#'   \item{lbslag}{}
#'   \item{spreadlag}{}
#'   \item{dwlag}{}
#'   \item{dcplag}{}
#' }
#' @source Perron, P. and Yamamoto, Y., 2015. "Using ols to estimate and test for structural changes in models with endogenous regressors." Journal of Applied Econometrics 30, 119–144.
"nkpc"

