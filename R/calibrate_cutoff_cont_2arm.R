#' Calibrate Posterior Probability Cutoff for a Two-Arm Comparative Trial with Continuous Endpoint
#'
#' The \code{calibrate_cutoff_cont_2arm} function is designed to calibrate the
#' posterior probability cutoff for a two-arm comparative trial with a
#' continuous endpoint under one borrowing strategy: self-adapting mixture prior
#' (SAM), robust MAP prior with fixed weight (rMAP), or non-informative prior
#' (NP).
#'
#' The calibrated cutoff is chosen so that the repeated-sampling rejection
#' probability under a specified scenario equals the target value.
#'
#' @param if.prior Informative prior constructed based on historical data for
#' the control arm, represented (approximately) as a normal mixture prior.
#' @param nf.prior Non-informative prior used for the mixture.
#' @param prior.t Prior used for the treatment arm. If missing, the default
#' value is set to be \code{nf.prior}.
#' @param target Target rejection probability under the calibration scenario,
#' typically the desired type I error rate when the calibration scenario
#' corresponds to the null hypothesis.
#' @param n.t Sample size for the treatment arm.
#' @param n Sample size for the control arm.
#' @param sigma.t Known sampling standard deviation for the treatment arm.
#' @param sigma Known sampling standard deviation for the control arm.
#' @param theta.t The mean for the treatment arm under the calibration scenario.
#' @param theta The mean for the control arm under the calibration scenario.
#' @param delta Clinically significant difference used for the SAM prior.
#' This argument is only used when \code{method = "SAM"}.
#' @param method Borrowing strategy for the control arm. Must be one of
#' \code{"SAM"}, \code{"rMAP"}, or \code{"NP"}.
#' @param alternative Direction of the posterior decision. Must be one of
#' \code{"greater"} (for superiority) or \code{"less"} (for inferiority).
#' @param margin Clinical margin. Must be non-negative. Default is \code{0}.
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
#' @param interval Search interval for the posterior probability cutoff. The
#' default is from 0.5 to 0.999.
#' @param rel.tol Tolerance passed to \code{\link[stats]{uniroot}}.
#' @param oc_rel.tol Relative tolerance passed to scenario-level numerical
#' integration.
#' @param n_sd_int Half-width of the numerical integration region for each arm,
#' expressed as a multiple of the corresponding standard error.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{cutoff}{Calibrated posterior probability cutoff.}
#'   \item{objective}{Value of the root-finding objective at the solution.}
#'   \item{target}{Target rejection probability.}
#'   \item{method}{Borrowing method used.}
#'   \item{alternative}{Direction of the posterior decision.}
#'   \item{margin}{Clinical margin used for inference.}
#'   \item{theta}{True control-arm mean under the calibration scenario.}
#'   \item{theta.t}{True treatment-arm mean under the calibration scenario.}
#'   \item{interval}{Search interval used for calibration.}
#' }
#'
#' @details
#' The function solves for the posterior probability cutoff such that the
#' scenario-level rejection probability returned by
#' \code{\link{eval_scenario_cont_2arm}} matches \code{target}. Under null
#' scenarios, this corresponds to calibration of the type I error.
#'
#' @export
calibrate_cutoff_cont_2arm <- function(if.prior, nf.prior, prior.t = nf.prior,
                                       target = 0.05,
                                       n.t, n,
                                       sigma.t, sigma,
                                       theta.t = NULL,
                                       theta = NULL,
                                       delta,
                                       method = c("SAM", "rMAP", "NP"),
                                       alternative = c("greater", "less"),
                                       margin = 0,
                                       weight_rMAP = 0.5,
                                       method.w = "LRT",
                                       prior.odds = 1,
                                       interval = c(0.5, 0.999),
                                       rel.tol = 1e-5,
                                       oc_rel.tol = 1e-6,
                                       n_sd_int = 8) {
  method <- match.arg(method)
  alternative <- match.arg(alternative)

  if (!is.numeric(target) || length(target) != 1 || !is.finite(target) ||
      target <= 0 || target >= 1) {
    stop("`target` must be a scalar in (0, 1).")
  }
  if (!is.numeric(interval) || length(interval) != 2 ||
      !all(is.finite(interval)) || interval[1] >= interval[2]) {
    stop("`interval` must be an increasing length-2 numeric vector.")
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
  if (!is.numeric(oc_rel.tol) || length(oc_rel.tol) != 1 ||
      !is.finite(oc_rel.tol) || oc_rel.tol <= 0) {
    stop("`oc_rel.tol` must be a positive scalar.")
  }
  if (!is.numeric(margin) || length(margin) != 1 || !is.finite(margin) || margin < 0) {
    stop("`margin` must be a non-negative scalar.")
  }
  if (is.null(theta.t)) theta.t <- summary(if.prior)["mean"]
  if (is.null(theta)) theta <- summary(if.prior)["mean"]

  f <- function(cutoff) {
    out <- eval_scenario_cont_2arm(
      if.prior = if.prior,
      nf.prior = nf.prior,
      prior.t = prior.t,
      n.t = n.t,
      n = n,
      sigma.t = sigma.t,
      sigma = sigma,
      theta.t = theta.t,
      theta = theta,
      cutoff = cutoff,
      delta = delta,
      method = method,
      alternative = alternative,
      margin = margin,
      weight_rMAP = weight_rMAP,
      method.w = method.w,
      prior.odds = prior.odds,
      rel.tol = oc_rel.tol,
      n_sd_int = n_sd_int
    )

    out$reject_prob - target
  }

  f_lo <- f(interval[1])
  f_hi <- f(interval[2])

  if (!is.finite(f_lo) || !is.finite(f_hi)) {
    stop("Calibration failed because the objective could not be evaluated at the interval endpoints.")
  }

  if (f_lo * f_hi > 0) {
    stop(
      sprintf(
        "Root not bracketed on [%.3f, %.3f]: f(lower)=%.6f, f(upper)=%.6f",
        interval[1], interval[2], f_lo, f_hi
      )
    )
  }

  out <- uniroot(f, interval = interval, tol = rel.tol)

  list(
    cutoff = out$root,
    objective = out$f.root,
    target = target,
    method = method,
    alternative = alternative,
    margin = margin,
    theta = theta,
    theta.t = theta.t,
    interval = interval
  )
}
