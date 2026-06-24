#' Fit the NPBIVPR model with regulatory network learning
#'
#' Fits an extension of the NPBIVP model for two omics platforms, where features
#' in the first platform may be regulated by features in the second platform.
#' The model was motivated by gene-miRNA regulatory network learning, where gene
#' expression features may be regulated by miRNAs, but the function can be used
#' for any two-platform setting with a directed candidate regulatory structure.
#'
#' @param data A list of exactly two omics data matrices. The first platform is
#'   treated as the regulated platform and the second as the regulatory platform.
#' @param y Numeric outcome vector.
#' @param outcome_type Either `"continuous"` or `"survival"`.
#' @param obs_ind Observation indicator for survival outcomes; 1 = observed,
#'   0 = right censored.
#' @param covar Optional clinical covariate matrix.
#' @param P_way Optional pathway/group membership matrix for features in the
#'   first platform. It should be a p1 by L 0/1 matrix. If `NULL`, all features
#'   in the first platform are assigned to one global group.
#' @param Cand Optional candidate regulatory matrix. It should be a p1 by p2
#'   0/1 matrix, where entry `(j,k)` indicates whether feature `k` in platform 2
#'   is a candidate regulator of feature `j` in platform 1. If `NULL`, all
#'   platform-2 features are treated as candidate regulators for all platform-1
#'   features.
#' @param nu_R Length-2 vector of probit-link intercepts controlling baseline
#'   prior inclusion probabilities for the two omics platforms before regulatory
#'   network adjustment. Default is `c(-1.5, -1.5)`.
#' @param r_muL_prior_param Shape and rate of the Gamma hyperprior for pathway
#'   effect precision. Default is `c(1, 1)`.
#' @param r_ev_prior_param Shape and rate of the Gamma hyperprior for error
#'   precision in the regulatory network model. Default is `c(1, 1)`.
#' @param R_prior_q Prior inclusion probability for candidate regulatory edges.
#'   Default is `0.025`.
#' @param tau_prior_param Optional length-2 vector of variances for half-normal
#'   priors on coupling parameters. If `NULL`, a heuristic based on `Cand` and
#'   `R_prior_q` is used.
#' @param beta_prior_tau Scale parameter controlling the variance of the normal
#'   spike prior for regulatory effects. Default is 10.
#' @param RN Logical; whether to jointly fit the regulatory network model.
#'   Default is `TRUE`. If `FALSE`, the model reduces to the two-platform NPBIVP
#'   outcome model with the same wrapper structure.
#' @inheritParams npbivp
#'
#' @return An object of class `c("npbivpr", "npbivp")`. It contains the usual
#'   NPBIVP output components plus, when `RN = TRUE`, regulatory-network
#'   summaries:
#' \describe{
#'   \item{samples}{Outcome-model posterior samples.}
#'   \item{acceptance}{Active acceptance-rate summary.}
#'   \item{acceptance_raw}{Raw acceptance-rate vector.}
#'   \item{nllik}{Negative log-likelihood values after burn-in.}
#'   \item{mu_samples}{Posterior samples of pathway/group effects.}
#'   \item{tau_samples}{Posterior samples of coupling parameters.}
#'   \item{RN_prob}{Posterior regulatory edge probability matrix.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 20
#' p1 <- 30
#' p2 <- 40
#'
#' X1 <- matrix(rnorm(n * p1), n, p1)
#' X2 <- matrix(rnorm(n * p2), n, p2)
#' y <- rnorm(n)
#'
#' P_way <- matrix(1, p1, 1)
#' Cand <- matrix(1, p1, p2)
#'
#' fit <- npbivpr(
#'   data = list(X1, X2),
#'   y = y,
#'   P_way = P_way,
#'   Cand = Cand,
#'   iter = 100,
#'   burn_in = 50
#' )
#'
#' fit$RN_prob
#' }
#'
#' @export
npbivpr <- function(data, y,
                    outcome_type = c("continuous", "survival"),
                    obs_ind = NULL,
                    covar = NULL,
                    P_way = NULL,
                    Cand = NULL,
                    nu_R = c(-1.5, -1.5),
                    r_muL_prior_param = c(1, 1),
                    r_ev_prior_param = c(1, 1),
                    R_prior_q = 0.025,
                    tau_prior_param = NULL,
                    beta_prior_tau = 10,
                    RN = TRUE,
                    scale_covar = TRUE,
                    scale_data = TRUE,
                    intercept = TRUE,
                    iter = 2000,
                    burn_in = 1000,
                    lambdas_prior_param = rep(1, 4),
                    r_prior_param = c(2, 0.1),
                    weight_prior = 1,
                    seed = 2026) {
  set.seed(seed)

  outcome_type <- match.arg(outcome_type)

  if (!is.list(data)) stop("data must be a list.")
  if (length(data) != 2) {
    stop("npbivpr requires exactly two omics platforms in data.")
  }

  data <- lapply(data, as.matrix)
  n_vec <- sapply(data, nrow)

  if (length(unique(n_vec)) != 1) {
    stop("Both omics matrices must have the same number of rows.")
  }

  n <- n_vec[1]
  y <- as.numeric(y)

  if (length(y) != n) {
    stop("length(y) must equal the number of rows in the omics matrices.")
  }

  p <- sapply(data, ncol)
  p1 <- p[1]
  p2 <- p[2]
  m <- 2

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
      stop("obs_ind must contain only 0/1 values.")
    }
  }

  if (is.null(P_way)) {
    P_way <- matrix(1, nrow = p1, ncol = 1)
  } else {
    P_way <- as.matrix(P_way)

    if (nrow(P_way) != p1) {
      stop("P_way must have nrow equal to the number of features in platform 1.")
    }

    if (!all(P_way %in% c(0, 1))) {
      stop("P_way must contain only 0/1 values.")
    }

    if (ncol(P_way) < 1) {
      stop("P_way must have at least one column.")
    }
  }

  if (is.null(Cand)) {
    Cand <- matrix(1, nrow = p1, ncol = p2)
  } else {
    Cand <- as.matrix(Cand)

    if (!all(dim(Cand) == c(p1, p2))) {
      stop("Cand must be a p1 by p2 matrix, where p1 = ncol(data[[1]]) and p2 = ncol(data[[2]]).")
    }

    if (!all(Cand %in% c(0, 1))) {
      stop("Cand must contain only 0/1 values.")
    }
  }

  if (length(nu_R) != 2) stop("nu_R must be a vector of length 2.")
  if (length(r_muL_prior_param) != 2) stop("r_muL_prior_param must have length 2.")
  if (length(r_ev_prior_param) != 2) stop("r_ev_prior_param must have length 2.")
  if (length(lambdas_prior_param) != 4) stop("lambdas_prior_param must have length 4.")
  if (length(r_prior_param) != 2) stop("r_prior_param must have length 2.")

  if (any(r_muL_prior_param <= 0)) stop("r_muL_prior_param values must be positive.")
  if (any(r_ev_prior_param <= 0)) stop("r_ev_prior_param values must be positive.")
  if (any(lambdas_prior_param <= 0)) stop("lambdas_prior_param values must be positive.")
  if (any(r_prior_param <= 0)) stop("r_prior_param values must be positive.")
  if (R_prior_q <= 0 || R_prior_q >= 1) stop("R_prior_q must be between 0 and 1.")
  if (beta_prior_tau <= 0) stop("beta_prior_tau must be positive.")

  if (is.null(tau_prior_param)) {
    u <- 1 / 3
    v <- 1 / 5

    xmax <- max(rowSums(Cand))
    ymax <- max(colSums(Cand))

    if (xmax == 0 || ymax == 0) {
      stop("Cand must contain at least one candidate edge in both row and column directions.")
    }

    tau_prior_param <- c(
      (u / (R_prior_q * xmax))^2,
      (v / (R_prior_q * ymax))^2
    )
  } else {
    tau_prior_param <- as.numeric(tau_prior_param)
    if (length(tau_prior_param) != 2) stop("tau_prior_param must have length 2.")
    if (any(tau_prior_param <= 0)) stop("tau_prior_param values must be positive.")
  }

  if (!is.logical(RN) || length(RN) != 1) stop("RN must be TRUE or FALSE.")
  if (!is.logical(scale_data) || length(scale_data) != 1) stop("scale_data must be TRUE or FALSE.")
  if (!is.logical(scale_covar) || length(scale_covar) != 1) stop("scale_covar must be TRUE or FALSE.")
  if (!is.logical(intercept) || length(intercept) != 1) stop("intercept must be TRUE or FALSE.")

  iter <- as.integer(iter)
  burn_in <- as.integer(burn_in)

  if (iter <= burn_in) stop("iter must be larger than burn_in.")

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

    m_fit <- 3
    c_lambdas <- c(1, rep(1, m_fit))
    sigma_prop <- rep(1, m_fit - 1)

    init <- c(
      rgamma(1, shape = lambdas_prior_param[3], rate = lambdas_prior_param[4]),
      rnorm(m_fit - 1),
      rgamma(1, shape = r_prior_param[1], rate = r_prior_param[2])
    )

    fit_raw <- npbivpr_cov_cpp(
      data_types, y, W, obs_indicator, M, burn_in,
      P_way, Cand, nu_R, lambdas_prior_param, r_prior_param,
      r_muL_prior_param, r_ev_prior_param, tau_prior_param,
      beta_prior_tau, R_prior_q, sigma_prior, sigma_prop,
      c_lambdas, c_prec, init, prop_scheme,
      weight_upd, icept_update, same_lambda, RN, FALSE
    )
  } else {
    W <- matrix(0, n, 3)
    m_fit <- 2
    c_lambdas <- c(1, rep(1, m_fit))
    sigma_prop <- rep(1, m_fit - 1)

    init <- c(
      rgamma(1, shape = lambdas_prior_param[3], rate = lambdas_prior_param[4]),
      rnorm(m_fit - 1),
      rgamma(1, shape = r_prior_param[1], rate = r_prior_param[2])
    )

    fit_raw <- npbivpr_nocov_cpp(
      data_types, y, W, obs_indicator, M, burn_in,
      P_way, Cand, nu_R, lambdas_prior_param, r_prior_param,
      r_muL_prior_param, r_ev_prior_param, tau_prior_param,
      beta_prior_tau, R_prior_q, sigma_prior, sigma_prop,
      c_lambdas, c_prec, init, prop_scheme,
      weight_upd, icept_update, same_lambda, RN, FALSE
    )
  }

  acc <- format_acceptance(fit_raw[[2]], m_fit)

  out <- list(
    samples = fit_raw[[1]],
    acceptance = acc$active,
    acceptance_raw = acc$raw,
    runtime = fit_raw[[3]],
    nllik = fit_raw[[4]],
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
    p = p,
    p_fit = if (!is.null(covar)) c(p, ncol(covar_train)) else p,
    iter = iter,
    burn_in = burn_in,
    scale_data = scale_data,
    scale_covar = scale_covar,
    RN = RN,
    P_way = P_way,
    Cand = Cand,
    nu_R = nu_R,
    R_prior_q = R_prior_q,
    tau_prior_param = tau_prior_param,
    beta_prior_tau = beta_prior_tau,
    mu_samples = if (RN) fit_raw[[5]] else NULL,
    tau_samples = if (RN) fit_raw[[6]] else NULL,
    RN_prob = if (RN) fit_raw[[7]] else NULL
  )

  class(out) <- c("npbivpr", "npbivp")
  out
}
