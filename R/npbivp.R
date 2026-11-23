#' Fit the NPBIVP model for multi-omics data
#'
#' Fits a Bayesian nonparametric integrative variable selection model for one or
#' more omics platforms. The model can be used with continuous outcomes or
#' right-censored survival outcomes. Clinical covariates may optionally be
#' included as an additional platform, but variable selection is not performed
#' on the covariates.
#'
#' @param data A list of omics data matrices. Each element of the list should
#'   contain one omics platform, with rows corresponding to individuals and
#'   columns corresponding to features.
#' @param y A numeric outcome vector of length \eqn{n}. For
#'   \code{outcome_type = "continuous"}, this is the quantitative outcome.
#'   For \code{outcome_type = "survival"}, this should be the log survival
#'   or log follow-up time.
#' @param outcome_type Character string specifying the outcome type. Choices are
#'   \code{"continuous"} and \code{"survival"}. Default is
#'   \code{"continuous"}.
#' @param obs_ind Observation indicator for survival outcomes. This should be a
#'   0/1 vector of length \eqn{n}, where 1 denotes an observed event time and
#'   0 denotes a right-censored observation. Required when
#'   \code{outcome_type = "survival"}. Ignored for continuous outcomes.
#' @param covar Optional matrix of clinical covariates, with rows corresponding
#'   to individuals and columns corresponding to covariates. If supplied, the
#'   covariates are treated as an additional platform. A platform weight is
#'   learned for the covariate block, but covariates are always included and no
#'   variable selection is performed on them.
#' @param scale_covar Logical; whether to min-max scale the covariates to
#'   \eqn{[0,1]}. Default is \code{TRUE}.
#' @param scale_data Logical; whether to min-max scale each omics platform to
#'   \eqn{[0,1]}. Default is \code{TRUE}.
#' @param intercept Logical; whether to include an intercept term in the outcome
#'   model. Default is \code{TRUE}.
#' @param iter Total number of MCMC iterations, including burn-in.
#' @param burn_in Number of initial MCMC iterations to discard.
#' @param pip Prior inclusion probability for features in each omics platform.
#'   If \code{NULL}, the default is \code{rep(0.025, m)}, where \eqn{m} is the
#'   number of omics platforms. A single numeric value applies the same prior
#'   inclusion probability to all platforms. A vector of length \eqn{m} allows
#'   platform-specific prior inclusion probabilities.
#' @param lambdas_prior_param Numeric vector of length 4 giving Gamma prior
#'   parameters for the intercept precision and common scale parameter. The form
#'   is \code{c(a1, b1, a2, b2)}, where \code{a1, b1} are the shape and rate for
#'   the intercept precision and \code{a2, b2} are the shape and rate for the
#'   common scale parameter. Default is \code{rep(1, 4)}.
#' @param r_prior_param Numeric vector of length 2 giving the shape and rate of
#'   the Gamma prior for the error precision. Default is \code{c(2, 0.1)}.
#' @param weight_prior Prior variance of the normal prior assigned to the
#'   transformed platform weights. Default is 1. Ignored when there is only one
#'   fitted platform.
#' @param seed Random seed used before initializing and running the MCMC sampler.
#'
#' @return An object of class \code{"npbivp"}, which is a list containing:
#' \describe{
#'   \item{\code{samples}}{Matrix of saved posterior samples after burn-in.}
#'   \item{\code{acceptance}}{Data frame containing only active acceptance rates.}
#'   \item{\code{acceptance_raw}}{Named numeric vector containing raw C++ acceptance output.}
#'   \item{\code{runtime}}{Total computing time for all MCMC iterations.}
#'   \item{\code{nllik}}{Vector of negative log-likelihood values after burn-in.}
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{outcome_type}}{The outcome type used in the model.}
#'   \item{\code{has_covar}}{Logical indicator of whether covariates were included.}
#'   \item{\code{data_train}}{Original training omics data.}
#'   \item{\code{y_train}}{Original training outcome vector.}
#'   \item{\code{obs_ind}}{Observation indicator used in model fitting.}
#'   \item{\code{covar_train}}{Original training covariate matrix, if supplied.}
#'   \item{\code{data_scale_params}}{Training-data min/max values used for omics scaling.}
#'   \item{\code{covar_scale_params}}{Training-data min/max values used for covariate scaling.}
#'   \item{\code{n}}{Number of individuals.}
#'   \item{\code{m}}{Number of omics platforms supplied by the user.}
#'   \item{\code{m_fit}}{Number of fitted platforms: \code{m} without covariates and
#'   \code{m + 1} with covariates.}
#'   \item{\code{p}}{Number of features in each omics platform.}
#'   \item{\code{p_fit}}{Number of features in each fitted platform, including the
#'   covariate block if supplied.}
#'   \item{\code{pip}}{Prior inclusion probabilities used in model fitting.}
#'   \item{\code{iter}}{Total number of MCMC iterations.}
#'   \item{\code{burn_in}}{Number of burn-in iterations.}
#'   \item{\code{scale_data}}{Logical indicating whether omics data were scaled.}
#'   \item{\code{scale_covar}}{Logical indicating whether covariates were scaled.}
#' }
#'
#' @details
#' If \code{scale_data = TRUE}, omics platforms are min-max scaled in R before
#' calling the C++ sampler, and the training min/max values are stored for use in
#' prediction. If \code{scale_covar = TRUE}, covariates are scaled similarly.
#' The C++ sampler is called with \code{scale_data = FALSE} because scaling is
#' handled by the R wrapper.
#'
#' If \code{covar} is supplied, the covariate matrix is appended internally as
#' the last fitted platform. Therefore, if the user supplies \eqn{m} omics
#' platforms, the fitted model has \eqn{m + 1} platforms.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 20
#' X1 <- matrix(rnorm(n * 30), n, 30)
#' X2 <- matrix(rnorm(n * 40), n, 40)
#' y <- rnorm(n)
#'
#' fit <- npbivp(
#'   data = list(X1, X2),
#'   y = y,
#'   iter = 100,
#'   burn_in = 50
#' )
#'
#' fit$acceptance
#' }
#'
#' @references Bhargob Kakoty and Thierry Chekouo. Integrative Nonparametric Bayesian Variable Selection and Prediction, with Application to Multi-Omics Cancer Genomics. Working paper.
#'
#' @export
npbivp <- function(data, y,
                   outcome_type = c("continuous", "survival"),
                   obs_ind = NULL,
                   covar = NULL,
                   scale_covar = TRUE,
                   scale_data = TRUE,
                   intercept = TRUE,
                   iter = 2000,
                   burn_in = 1000,
                   pip = NULL,
                   lambdas_prior_param = rep(1, 4),
                   r_prior_param = c(2, 0.1),
                   weight_prior = 1,
                   seed = 2026) {
  set.seed(seed)

  outcome_type <- match.arg(outcome_type)

  if (!is.list(data)) stop("data must be a list of omics data matrices.")

  m <- length(data)
  if (m < 1) stop("data must contain at least one omics platform.")

  data <- lapply(data, as.matrix)
  n_vec <- sapply(data, nrow)

  if (length(unique(n_vec)) != 1) {
    stop("All omics data matrices must have the same number of rows.")
  }

  n <- n_vec[1]
  y <- as.numeric(y)

  if (length(y) != n) {
    stop("length(y) must equal the number of rows in each omics data matrix.")
  }

  if (!is.logical(scale_data) || length(scale_data) != 1) {
    stop("scale_data must be TRUE or FALSE.")
  }

  if (!is.logical(scale_covar) || length(scale_covar) != 1) {
    stop("scale_covar must be TRUE or FALSE.")
  }

  if (!is.logical(intercept) || length(intercept) != 1) {
    stop("intercept must be TRUE or FALSE.")
  }

  if (!is.numeric(iter) || length(iter) != 1 || iter <= 0) {
    stop("iter must be a positive integer.")
  }

  if (!is.numeric(burn_in) || length(burn_in) != 1 || burn_in < 0) {
    stop("burn_in must be a non-negative integer.")
  }

  iter <- as.integer(iter)
  burn_in <- as.integer(burn_in)

  if (iter <= burn_in) stop("iter must be larger than burn_in.")

  if (length(lambdas_prior_param) != 4) {
    stop("lambdas_prior_param must be a vector of length 4: c(a1,b1,a2,b2).")
  }

  if (length(r_prior_param) != 2) {
    stop("r_prior_param must be a vector of length 2: c(shape, rate).")
  }

  if (any(lambdas_prior_param <= 0)) {
    stop("lambdas_prior_param values must be positive.")
  }

  if (any(r_prior_param <= 0)) {
    stop("r_prior_param values must be positive.")
  }

  if (weight_prior <= 0) stop("weight_prior must be positive.")

  if (outcome_type == "continuous") {
    if (!is.null(obs_ind)) {
      warning("obs_ind was provided but will be ignored for continuous outcome.")
    }
    obs_ind <- rep(1L, n)
  }

  if (outcome_type == "survival") {
    if (is.null(obs_ind)) {
      stop("For outcome_type = 'survival', obs_ind must be provided.")
    }

    obs_ind <- as.integer(obs_ind)

    if (length(obs_ind) != n) {
      stop("obs_ind must have length equal to length(y).")
    }

    if (!all(obs_ind %in% c(0L, 1L))) {
      stop("obs_ind must contain only 0/1 values: 1 = observed, 0 = right censored.")
    }
  }

  if (is.null(pip)) {
    alpha_prior_param <- rep(0.025, m)
  } else {
    alpha_prior_param <- as.numeric(pip)

    if (length(alpha_prior_param) == 1) {
      alpha_prior_param <- rep(alpha_prior_param, m)
    } else if (length(alpha_prior_param) != m) {
      stop("pip must be either a single value or a vector of length ", m, ".")
    }
  }

  if (any(alpha_prior_param <= 0 | alpha_prior_param >= 1)) {
    stop("pip values must be between 0 and 1.")
  }

  data_scale_params <- lapply(data, get_minmax_params)

  if (scale_data) {
    data_types <- Map(apply_minmax_params, data, data_scale_params)
  } else {
    data_types <- data
  }

  obs_indicator <- obs_ind
  M <- iter

  sigma_prior <- weight_prior
  c_prec <- 1
  prop_scheme <- 2
  weight_upd <- TRUE
  icept_update <- intercept
  same_lambda <- TRUE

  covar_train <- NULL
  covar_scale_params <- NULL

  if (!is.null(covar)) {
    W_original <- as.matrix(covar)

    if (nrow(W_original) != n) {
      stop("covar must have the same number of rows as the omics data.")
    }

    covar_train <- W_original
    covar_scale_params <- get_minmax_params(W_original)

    if (scale_covar) {
      W <- apply_minmax_params(W_original, covar_scale_params)
    } else {
      W <- W_original
    }

    alpha_prior_param <- c(alpha_prior_param, 1)

    c_lambdas <- c(1, rep(1, m + 1))
    sigma_prop <- rep(1, m)

    init <- c(
      rgamma(1, shape = lambdas_prior_param[3], rate = lambdas_prior_param[4]),
      rnorm(m),
      rgamma(1, shape = r_prior_param[1], rate = r_prior_param[2])
    )

    samples <- npbivp_cov_cpp(
      data_types, y, W, obs_indicator, M, burn_in,
      alpha_prior_param, lambdas_prior_param, r_prior_param,
      sigma_prior, sigma_prop, c_lambdas, c_prec, init,
      prop_scheme, weight_upd, icept_update, same_lambda,
      FALSE
    )
  } else {
    W <- matrix(0, n, 3)
    c_lambdas <- c(1, rep(1, m))

    if (m == 1) {
      sigma_prop <- 0
      init <- c(
        rgamma(1, shape = lambdas_prior_param[3], rate = lambdas_prior_param[4]),
        rgamma(1, shape = r_prior_param[1], rate = r_prior_param[2])
      )
    } else {
      sigma_prop <- rep(1, m - 1)
      init <- c(
        rgamma(1, shape = lambdas_prior_param[3], rate = lambdas_prior_param[4]),
        rnorm(m - 1),
        rgamma(1, shape = r_prior_param[1], rate = r_prior_param[2])
      )
    }

    samples <- npbivp_nocov_cpp(
      data_types, y, W, obs_indicator, M, burn_in,
      alpha_prior_param, lambdas_prior_param, r_prior_param,
      sigma_prior, sigma_prop, c_lambdas, c_prec, init,
      prop_scheme, weight_upd, icept_update, same_lambda,
      FALSE
    )
  }

  m_fit <- if (!is.null(covar)) m + 1 else m


  acc <- format_acceptance(samples[[2]], m_fit)

  out <- list(
    samples = samples[[1]],
    acceptance = acc$active,
    acceptance_raw = acc$raw,
    runtime = samples[[3]],
    nllik = samples[[4]],
    call = match.call(),
    outcome_type = outcome_type,
    has_covar = !is.null(covar),
    data_train = data,
    y_train = y,
    obs_ind = obs_ind,
    covar_train = covar_train,
    data_scale_params = data_scale_params,
    covar_scale_params = covar_scale_params,
    n = n,
    m = m,
    m_fit = m_fit,
    p = sapply(data, ncol),
    p_fit = if (!is.null(covar)) c(sapply(data, ncol), ncol(covar_train)) else sapply(data, ncol),
    pip = alpha_prior_param,
    iter = iter,
    burn_in = burn_in,
    scale_data = scale_data,
    scale_covar = scale_covar
  )

  class(out) <- "npbivp"
  out
}
