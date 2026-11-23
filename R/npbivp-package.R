#' npbivp: Bayesian Nonparametric Integrative Variable Selection
#'
#' Bayesian nonparametric integrative variable selection methods for
#' multi-omics data with continuous and survival outcomes.
#'
#' Main functions:
#'
#' \itemize{
#' \item \code{\link{npbivp}}
#' \item \code{\link{npbivpr}}
#' \item \code{\link{predict.npbivp}}
#' }
#'
#' @references Bhargob Kakoty and Thierry Chekouo. Integrative Nonparametric Bayesian Variable Selection and Prediction, with Application to Multi-Omics Cancer Genomics. Working paper.
#'
#' @keywords internal
#' @useDynLib npbivp, .registration = TRUE
#' @importFrom stats predict rgamma rnorm
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
