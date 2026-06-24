#' Predict from a fitted NPBIVP model
#'
#' Generates posterior predictive conditional means and variances for new
#' observations from a fitted \code{"npbivp"} object.
#'
#' @param object A fitted object returned by \code{\link{npbivp}}.
#' @param newdata A list of test omics data matrices. The list must have the
#'   same number of omics platforms and same number of columns per platform as
#'   the training data used in \code{object}.
#' @param newcovar Optional matrix of test clinical covariates. Required if the
#'   model was fitted with covariates.
#' @param thin Thinning interval for posterior samples used in prediction.
#'   Default is 5.
#' @param start First saved MCMC sample index to use. Default is 1.
#' @param end Last saved MCMC sample index to use. Default is \code{NULL}, which
#'   uses all available saved samples.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"npbivp_pred"}, which is a list containing:
#' \describe{
#'   \item{\code{pred_draws}}{Matrix of posterior predictive conditional means.
#'   Rows correspond to test observations and columns correspond to retained
#'   MCMC samples.}
#'   \item{\code{cond_var_draws}}{Matrix of posterior predictive conditional
#'   variances. Rows correspond to test observations and columns correspond to
#'   retained MCMC samples.}
#'   \item{\code{sample_index}}{Indices of saved MCMC samples used for prediction.}
#'   \item{\code{call}}{The matched function call.}
#' }
#'
#' @details
#' If the original model was fitted with \code{scale_data = TRUE}, the new omics
#' data are min-max scaled using the training-data minima and maxima stored in
#' the fitted object. Similarly, if covariates were used and
#' \code{scale_covar = TRUE}, \code{newcovar} is scaled using the training
#' covariate minima and maxima. The prediction C++ routine is called with
#' \code{scale_data = FALSE} because scaling is handled by the R wrapper.
#'
#' For survival outcomes, censored training log-times are replaced by their
#' posterior mean imputed values before prediction.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 20
#' n_test <- 5
#'
#' X1 <- matrix(rnorm(n * 30), n, 30)
#' X1_test <- matrix(rnorm(n_test * 30), n_test, 30)
#' y <- rnorm(n)
#'
#' fit <- npbivp(
#'   data = list(X1),
#'   y = y,
#'   iter = 100,
#'   burn_in = 50
#' )
#'
#' pred <- predict(fit, newdata = list(X1_test), thin = 2)
#' dim(pred$pred_draws)
#' }
#'
#' @export
predict.npbivp <- function(object,
                           newdata,
                           newcovar = NULL,
                           thin = 5,
                           start = 1,
                           end = NULL,
                           ...) {
  if (!inherits(object, "npbivp")) {
    stop("object must be a fitted npbivp object.")
  }

  if (!is.list(newdata)) {
    stop("newdata must be a list of omics data matrices.")
  }

  if (length(newdata) != object$m) {
    stop("newdata must contain the same number of omics platforms as the training data.")
  }

  newdata <- lapply(newdata, as.matrix)

  if (length(unique(sapply(newdata, nrow))) != 1) {
    stop("All test omics matrices in newdata must have the same number of rows.")
  }

  for (j in seq_len(object$m)) {
    if (ncol(newdata[[j]]) != object$p[j]) {
      stop(
        paste0(
          "Platform ", j, " in newdata has ", ncol(newdata[[j]]),
          " columns, but the training data had ", object$p[j], "."
        )
      )
    }
  }

  n_test <- nrow(newdata[[1]])

  if (object$has_covar) {
    if (is.null(newcovar)) {
      stop("newcovar must be supplied because the model was fitted with covariates.")
    }

    newcovar <- as.matrix(newcovar)

    if (nrow(newcovar) != n_test) {
      stop("newcovar must have the same number of rows as newdata.")
    }

    if (ncol(newcovar) != ncol(object$covar_train)) {
      stop("newcovar must have the same number of columns as the training covariate matrix.")
    }
  }

  if (!object$has_covar && !is.null(newcovar)) {
    warning("newcovar was supplied but will be ignored because the model was fitted without covariates.")
  }

  if (!is.numeric(thin) || length(thin) != 1 || thin <= 0) {
    stop("thin must be a positive integer.")
  }

  if (!is.numeric(start) || length(start) != 1 || start <= 0) {
    stop("start must be a positive integer.")
  }

  thin <- as.integer(thin)
  start <- as.integer(start)

  samples <- object$samples
  n_saved <- nrow(samples)

  if (is.null(end)) {
    end <- n_saved
  }

  if (!is.numeric(end) || length(end) != 1 || end <= 0) {
    stop("end must be a positive integer or NULL.")
  }

  end <- as.integer(end)

  if (start > end) {
    stop("start must be less than or equal to end.")
  }

  sample_index <- seq(start, end, by = thin)

  if (any(sample_index < 1) || any(sample_index > n_saved)) {
    stop("Invalid start/end/thin combination.")
  }

  if (object$scale_data) {
    data_train_cpp <- Map(
      apply_minmax_params,
      object$data_train,
      object$data_scale_params
    )

    newdata_cpp <- Map(
      apply_minmax_params,
      newdata,
      object$data_scale_params
    )
  } else {
    data_train_cpp <- object$data_train
    newdata_cpp <- newdata
  }

  if (object$has_covar) {
    if (object$scale_covar) {
      covar_train_cpp <- apply_minmax_params(
        object$covar_train,
        object$covar_scale_params
      )

      newcovar_cpp <- apply_minmax_params(
        newcovar,
        object$covar_scale_params
      )
    } else {
      covar_train_cpp <- object$covar_train
      newcovar_cpp <- newcovar
    }
  }

  p_fit <- object$p_fit
  m_fit <- object$m_fit
  total_p_fit <- sum(p_fit)

  y_train_pred <- object$y_train

  if (object$outcome_type == "survival") {
    cens_id <- which(object$obs_ind == 0)

    if (length(cens_id) > 0) {
      cens_start <- 2 * total_p_fit + 2 * m_fit + 3
      cens_end <- cens_start + length(cens_id) - 1

      cens_draws <- samples[sample_index, cens_start:cens_end, drop = FALSE]
      cens_mean <- colMeans(cens_draws)

      y_train_pred[cens_id] <- cens_mean
    }
  }

  rho_start <- total_p_fit + 1
  rho_end <- 2 * total_p_fit

  param_start <- 2 * total_p_fit + 1
  param_end <- 2 * total_p_fit + 2 * m_fit + 2

  samples_pred <- samples[
    sample_index,
    c(rho_start:rho_end, param_start:param_end),
    drop = FALSE
  ]

  if (object$has_covar) {
    pred <- pred_cov_cpp(
      data_types_train = data_train_cpp,
      data_types_test = newdata_cpp,
      y = y_train_pred,
      W_train = covar_train_cpp,
      W_test = newcovar_cpp,
      samples = samples_pred,
      scale_data = FALSE
    )
  } else {
    W_train <- matrix(0, object$n, 3)
    W_test <- matrix(0, n_test, 3)

    pred <- pred_nocov_cpp(
      data_types_train = data_train_cpp,
      data_types_test = newdata_cpp,
      y = y_train_pred,
      W_train = W_train,
      W_test = W_test,
      samples = samples_pred,
      scale_data = FALSE
    )
  }

  pred_draws <- pred[[1]]
  cond_var_draws <- pred[[2]]

  out <- list(
    pred_draws = pred_draws,
    cond_var_draws = cond_var_draws,
    sample_index = sample_index,
    call = match.call()
  )

  class(out) <- "npbivp_pred"
  out
}


#' @rdname predict.npbivp
#' @export
npbivp_pred <- function(object, ...) {
  predict(object, ...)
}
