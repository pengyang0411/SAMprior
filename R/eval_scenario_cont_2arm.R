#' Evaluate One Scenario for a Two-Arm Comparative Trial with Continuous Endpoint
#'
#' The \code{eval_scenario_cont_2arm} function is designed to evaluate
#' repeated-sampling operating characteristics for a two-arm comparative trial
#' with a continuous endpoint under one borrowing strategy: self-adapting
#' mixture prior (SAM), robust MAP prior with fixed weight (rMAP), or
#' non-informative prior (NP).
#'
#' The treatment effect is defined as \eqn{\tau = \theta_t - \theta},
#' where \eqn{\theta_t} and \eqn{\theta} denote the true means in the
#' treatment and control arms, respectively.
#'
#' For a given true scenario \eqn{(\theta_t, \theta)}, this function computes
#' the repeated-sampling rejection probability, bias, root mean squared error
#' (RMSE), and mean borrowing weight using one-dimensional numerical
#' integration.
#'
#' @param if.prior Informative prior constructed based on historical data for
#' the control arm, represented (approximately) as a normal mixture prior.
#' @param nf.prior Non-informative prior used as the robustifying component
#' for the control arm prior.
#' @param prior.t Prior used for the treatment arm. If missing, the default
#' value is set to be \code{nf.prior}.
#' @param n.t Number of subjects in the treatment arm.
#' @param n Number of subjects in the control arm.
#' @param sigma.t Known sampling standard deviation in the treatment arm.
#' @param sigma Known sampling standard deviation in the control arm.
#' @param theta.t True treatment arm mean.
#' @param theta True control arm mean.
#' @param cutoff Posterior probability cutoff used for decision making.
#' Rejection occurs if the posterior probability exceeds \code{cutoff}.
#' @param delta Clinically significant difference used for the SAM prior.
#' This argument is only used when \code{method = "SAM"}.
#' @param method Borrowing strategy for the control arm. Must be one of
#' \code{"SAM"}, \code{"rMAP"}, or \code{"NP"}.
#' @param alternative Direction of the posterior decision. Must be one of
#' \code{"greater"} (for superiority) or \code{"less"} (for inferiority).
#' @param margin Clinical margin. Must be a non-negative scalar. The default
#' value is \code{0}.
#' @param weight_rMAP Weight assigned to the informative prior component
#' (\eqn{0 \leq} \code{weight_rMAP} \eqn{\leq 1}) for the robust MAP prior.
#' This argument is only used when \code{method = "rMAP"}. The default value is
#' 0.5.
#' @param method.w Methods used to determine the mixture weight for SAM priors.
#' The default method is "LRT" (Likelihood Ratio Test), the alternative option
#' is "PPR" (Posterior Probability Ratio). See \code{\link{SAM_weight}} for
#' more details.
#' @param prior.odds The prior probability of \eqn{H_0} being true compared to
#' the prior probability of \eqn{H_1} being true using PPR method. The default
#' value is 1. See \code{\link{SAM_weight}} for more details.
#' @param rel.tol Relative tolerance passed to numerical integration.
#' @param n_sd_int Half-width of the numerical integration region for each arm,
#' expressed as a multiple of the corresponding standard error.
#'
#' @return A one-row data frame with the following columns:
#' \describe{
#'   \item{theta}{True control arm mean.}
#'   \item{theta.t}{True treatment arm mean.}
#'   \item{delta_true}{True treatment effect, \eqn{\tau = \theta_t - \theta}.}
#'   \item{method}{Borrowing method used.}
#'   \item{alternative}{Direction of the posterior decision.}
#'   \item{cutoff}{Posterior probability cutoff used for decision making.}
#'   \item{margin}{Clinical margin used for inference.}
#'   \item{reject_prob}{Repeated-sampling rejection probability.}
#'   \item{bias}{Bias of the posterior mean estimator of \eqn{\theta}.}
#'   \item{rmse}{Root mean squared error of the posterior mean estimator of \eqn{\theta}.}
#'   \item{mean_weight}{Average borrowing weight under the specified method.}
#' }
#'
#' @details
#' The rejection probability is computed by reducing the repeated-sampling
#' decision rule to a one-dimensional integral over the control-arm sample mean.
#'
#' Bias and RMSE are evaluated for the posterior mean estimator of the control
#' arm mean \eqn{\theta}. Both are computed from one-dimensional first and
#' second moments of the control-arm posterior mean. The mean borrowing weight
#' is computed by one-dimensional integration over the control-arm sample mean.
#'
#' Under null scenarios, \code{reject_prob} corresponds to type I error. Under
#' alternative scenarios, it corresponds to power.
#'
#' @export
eval_scenario_cont_2arm <- function(if.prior, nf.prior, prior.t = nf.prior,
                                    n.t, n,
                                    sigma.t, sigma,
                                    theta.t, theta,
                                    cutoff,
                                    delta,
                                    method = c("SAM", "rMAP", "NP"),
                                    alternative = c("greater", "less"),
                                    margin = 0,
                                    weight_rMAP = 0.5,
                                    method.w = "LRT",
                                    prior.odds = 1,
                                    rel.tol = 1e-6,
                                    n_sd_int = 8) {
  method <- match.arg(method)
  alternative <- match.arg(alternative)

  if (!is.numeric(n.t) || length(n.t) != 1 || !is.finite(n.t) || n.t <= 0) {
    stop("`n.t` must be a positive scalar.")
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    stop("`n` must be a positive scalar.")
  }
  if (!is.numeric(sigma.t) || length(sigma.t) != 1 || !is.finite(sigma.t) || sigma.t <= 0) {
    stop("`sigma.t` must be a positive scalar.")
  }
  if (!is.numeric(sigma) || length(sigma) != 1 || !is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be a positive scalar.")
  }
  if (!is.numeric(theta.t) || length(theta.t) != 1 || !is.finite(theta.t)) {
    stop("`theta.t` must be a finite scalar.")
  }
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta)) {
    stop("`theta` must be a finite scalar.")
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff <= 0 || cutoff >= 1) {
    stop("`cutoff` must be a scalar in (0, 1).")
  }
  if (!is.numeric(margin) || length(margin) != 1 || !is.finite(margin) || margin < 0) {
    stop("`margin` must be a non-negative scalar.")
  }
  if (!method.w %in% c("LRT", "PPR")) {
    stop("`method.w` must be either \"LRT\" or \"PPR\".")
  }
  if (!is.numeric(prior.odds) || length(prior.odds) != 1 ||
      !is.finite(prior.odds) || prior.odds <= 0) {
    stop("`prior.odds` must be a positive scalar.")
  }
  if (!is.numeric(rel.tol) || length(rel.tol) != 1 || !is.finite(rel.tol) || rel.tol <= 0) {
    stop("`rel.tol` must be a positive scalar.")
  }
  if (!is.numeric(n_sd_int) || length(n_sd_int) != 1 || !is.finite(n_sd_int) || n_sd_int <= 0) {
    stop("`n_sd_int` must be a positive scalar.")
  }
  if (method == "rMAP") {
    if (is.null(weight_rMAP) || length(weight_rMAP) != 1 ||
        !is.finite(weight_rMAP) || weight_rMAP < 0 || weight_rMAP > 1) {
      stop("`weight_rMAP` must be a scalar in [0, 1] when `method = \"rMAP\"`.")
    }
  }

  se_t <- sigma.t / sqrt(n.t)
  se   <- sigma / sqrt(n)

  lower_t <- theta.t - n_sd_int * se_t
  upper_t <- theta.t + n_sd_int * se_t
  lower   <- theta   - n_sd_int * se
  upper   <- theta   + n_sd_int * se

  ## ---- control posterior object given ybar ----
  post_c_fun <- function(ybar) {
    if (method == "SAM") {
      .control_sam_update(
        if.prior = if.prior,
        nf.prior = nf.prior,
        ybar = ybar,
        n = n,
        sigma = sigma,
        delta = delta,
        weight_fun = weight_fun_normmix,
        method.w = method.w,
        prior.odds = prior.odds
      )
    } else if (method == "rMAP") {
      .control_fixed_update(
        if.prior = if.prior,
        nf.prior = nf.prior,
        ybar = ybar,
        n = n,
        sigma = sigma,
        weight = weight_rMAP
      )
    } else {
      post_c <- .posterior_mix_update(
        prior = nf.prior,
        ybar = ybar,
        n = n,
        sigma = sigma
      )
      post_c$weight <- 0
      post_c
    }
  }

  ## ---- posterior mean of theta given ybar ----
  post_mean_fun <- function(ybar) {
    post_c <- post_c_fun(ybar)
    sum(post_c$w * post_c$mu)
  }

  ## ---- borrowing weight given ybar ----
  weight_fun_eval <- function(ybar) {
    if (method == "SAM") {
      weight_fun_normmix(
        ybar = ybar,
        if.prior = if.prior,
        nf.prior = nf.prior,
        n = n,
        sigma = sigma,
        delta = delta,
        method.w = method.w,
        prior.odds = prior.odds
      )
    } else if (method == "rMAP") {
      weight_rMAP
    } else {
      0
    }
  }

  ## ---- posterior probability threshold in ybar_t for fixed ybar ----
  find_ybar_t_cutoff_general <- function(ybar) {
    gfun <- function(ybar_t) {
      ps <- post_summary_cont_2arm(
        ybar_t = ybar_t,
        ybar = ybar,
        if.prior = if.prior,
        nf.prior = nf.prior,
        prior.t = prior.t,
        n.t = n.t,
        n = n,
        sigma.t = sigma.t,
        sigma = sigma,
        delta = delta,
        cutoff = cutoff,
        method = method,
        alternative = alternative,
        margin = margin,
        weight_rMAP = weight_rMAP,
        method.w = method.w,
        prior.odds = prior.odds
      )
      ps$post_prob - cutoff
    }

    gl <- gfun(lower_t)
    gu <- gfun(upper_t)

    if (alternative == "greater") {
      if (gl > 0) return(-Inf)
      if (gu < 0) return( Inf)
    } else {
      if (gl < 0) return(-Inf)
      if (gu > 0) return( Inf)
    }

    uniroot(gfun, interval = c(lower_t, upper_t), tol = 1e-8)$root
  }

  ## ---- 1D reject probability over ybar ----
  reject_integrand <- function(ybar) {
    vapply(ybar, function(x) {
      y_cut <- find_ybar_t_cutoff_general(x)

      pr_reject_given_y <-
        if (is.infinite(y_cut) && y_cut < 0) {
          if (alternative == "greater") 1 else 0
        } else if (is.infinite(y_cut) && y_cut > 0) {
          if (alternative == "greater") 0 else 1
        } else if (alternative == "greater") {
          1 - pnorm(y_cut, mean = theta.t, sd = se_t)
        } else {
          pnorm(y_cut, mean = theta.t, sd = se_t)
        }

      pr_reject_given_y * dnorm(x, mean = theta, sd = se)
    }, numeric(1))
  }

  reject_prob <- integrate(
    f = reject_integrand,
    lower = lower,
    upper = upper,
    rel.tol = rel.tol,
    subdivisions = 200L
  )$value

  ## ---- 1D control moments for bias/RMSE ----
  mean_integrand <- function(ybar) {
    vapply(ybar, function(x) {
      post_mean_fun(x) * dnorm(x, mean = theta, sd = se)
    }, numeric(1))
  }

  mean_sq_integrand <- function(ybar) {
    vapply(ybar, function(x) {
      post_mean_fun(x)^2 * dnorm(x, mean = theta, sd = se)
    }, numeric(1))
  }

  weight_integrand <- function(ybar) {
    vapply(ybar, function(x) {
      weight_fun_eval(x) * dnorm(x, mean = theta, sd = se)
    }, numeric(1))
  }

  mean_post <- integrate(
    f = mean_integrand,
    lower = lower,
    upper = upper,
    rel.tol = rel.tol,
    subdivisions = 200L
  )$value

  mean_sq_post <- integrate(
    f = mean_sq_integrand,
    lower = lower,
    upper = upper,
    rel.tol = rel.tol,
    subdivisions = 200L
  )$value

  mean_weight <- integrate(
    f = weight_integrand,
    lower = lower,
    upper = upper,
    rel.tol = rel.tol,
    subdivisions = 200L
  )$value

  ## ---- bias and RMSE for theta ----
  bias <- mean_post - theta
  rmse <- sqrt(pmax(mean_sq_post - 2 * theta * mean_post + theta^2, 0))

  data.frame(
    theta = theta,
    theta.t = theta.t,
    delta_true = theta.t - theta,
    method = method,
    alternative = alternative,
    cutoff = cutoff,
    margin = margin,
    reject_prob = reject_prob,
    bias = bias,
    rmse = rmse,
    mean_weight = mean_weight
  )
}
