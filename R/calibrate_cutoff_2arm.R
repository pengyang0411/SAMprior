#' Calibrate Posterior Probability Cutoff for a Two-Arm Comparative Trial
#'
#' The \code{calibrate_cutoff_2arm} function is designed to calibrate the
#' posterior probability cutoff for a two-arm comparative trial under the
#' specified endpoint type and borrowing strategy: self-adapting mixture
#' prior (SAM), robust MAP prior with fixed weight (rMAP), or non-informative
#' prior (NP).
#'
#' For binary endpoints, the function dispatches to
#' \code{\link{calibrate_cutoff_bin_2arm}}. For continuous endpoints, it
#' dispatches to \code{\link{calibrate_cutoff_cont_2arm}}.
#'
#' @param if.prior Informative prior constructed from historical data,
#' represented (approximately) as a mixture of conjugate distributions.
#' @param nf.prior Non-informative prior used for the mixture and as the
#' robustifying component for the control arm prior.
#' @param prior.t Prior used for the treatment arm. If missing, the default
#' value is set to be \code{nf.prior}.
#' @param target Target rejection probability under the calibration scenario,
#' typically the desired type I error rate when the calibration scenario
#' corresponds to the null hypothesis.
#' @param n.t Sample size for the treatment arm.
#' @param n Sample size for the control arm.
#' @param theta.t The response rate (binary endpoints) or mean
#' (continuous endpoints) for the treatment arm under the calibration scenario.
#' @param theta The response rate (binary endpoints) or mean
#' (continuous endpoints) for the control arm under the calibration scenario.
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
#' @param oc_rel.tol Relative tolerance passed to scenario-level calculations.
#' For continuous endpoints this controls numerical integration in the
#' operating-characteristic engine; for binary endpoints this controls numerical
#' integration in posterior probability calculations.
#' @param ... Additional arguments required by the endpoint-specific method.
#' For continuous endpoints, this includes \code{sigma.t}, \code{sigma}, and
#' optionally \code{n_sd_int}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{cutoff}{Calibrated posterior probability cutoff.}
#'   \item{objective}{Value of the root-finding objective at the solution.}
#'   \item{target}{Target rejection probability.}
#'   \item{method}{Borrowing method used.}
#'   \item{alternative}{Direction of the posterior decision.}
#'   \item{margin}{Clinical margin used for inference.}
#'   \item{theta}{True control-arm parameter under the calibration scenario.}
#'   \item{theta.t}{True treatment-arm parameter under the calibration scenario.}
#'   \item{interval}{Search interval used for calibration.}
#' }
#'
#' @details
#' The function calibrates the posterior probability cutoff \eqn{C} for a
#' two-arm comparative trial so that the repeated-sampling rejection
#' probability under the specified calibration scenario equals
#' \code{target}. Inference is based on the posterior probability of the
#' treatment effect \eqn{\theta_t - \theta}, where \eqn{\theta_t} and
#' \eqn{\theta} denote the treatment-arm and control-arm parameters,
#' respectively.
#'
#' For \code{alternative = "greater"}, the posterior decision rule is based on
#' \deqn{P(\theta_t - \theta > \Delta \mid D) > C,}
#' where \eqn{\Delta} denotes the clinical margin and \eqn{D} denotes the
#' observed data. Under the default settings, calibration is performed under
#' the null scenario \eqn{\theta_t = \theta = \theta_h}, where
#' \eqn{\theta_h = E(\code{if.prior})} is the posterior mean of the
#' informative prior constructed from historical data. This corresponds to a
#' superiority test when \eqn{\Delta = 0}.
#'
#' For \code{alternative = "less"}, the posterior decision rule is based on
#' \deqn{P(\theta_t - \theta < -\Delta \mid D) > C.}
#' Calibration is again performed under the null scenario
#' \eqn{\theta_t = \theta = \theta_h}, where
#' \eqn{\theta_h = E(\code{if.prior})}. This corresponds to an
#' inferiority test when \eqn{\Delta = 0}.
#'
#' This function acts as a generic front-end and dispatches according to the
#' class of \code{if.prior}. For binary endpoints, it calls
#' \code{\link{calibrate_cutoff_bin_2arm}}. For continuous endpoints, it calls
#' \code{\link{calibrate_cutoff_cont_2arm}}.
#'
#' @seealso \code{\link{calibrate_cutoff_bin_2arm}},
#' \code{\link{calibrate_cutoff_cont_2arm}}
#'
#' @examples
#' ## Example: calibrate the posterior probability cutoff for a two-arm
#' ## binary trial using a SAM prior for the control arm
#'
#' ## Informative prior constructed from historical data
#' if.prior <- mixbeta(c(1, 20, 40))
#'
#' ## Sample sizes
#' n.t <- 100
#' n <- 50
#'
#' ## Calibrate the posterior probability cutoff
#' PPC <- calibrate_cutoff_2arm(
#'   if.prior = if.prior,             ## Informative prior from historical data
#'   nf.prior = mixbeta(c(1, 1, 1)),  ## Non-informative prior for SAM prior
#'   prior.t = mixbeta(c(1, 1, 1)),   ## Prior for treatment arm
#'   target = 0.05,                   ## Targeted type I error rate
#'   n.t = n.t,                       ## Sample size for treatment arm
#'   n = n,                           ## Sample size for control arm
#'   theta.t = summary(if.prior)["mean"],  ## True treatment-arm response rate
#'   theta = summary(if.prior)["mean"],    ## True control-arm response rate
#'   method = "SAM",                  ## Borrowing method
#'   delta = 0.2,                     ## Clinically significant difference
#'   alternative = "greater",         ## Superiority test
#'   margin = 0                       ## Clinical margin
#' )
#'
#' PPC
#'
#'
#' @export
calibrate_cutoff_2arm <- function(if.prior, nf.prior, prior.t = nf.prior,
                                  target = 0.05,
                                  n.t, n,
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
                                  ...) {
  UseMethod("calibrate_cutoff_2arm")
}

#' @export
calibrate_cutoff_2arm.default <- function(if.prior, nf.prior, prior.t = nf.prior,
                                          target = 0.05,
                                          n.t, n,
                                          theta.t = NULL, theta = NULL,
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
                                          ...) {
  stop("Unsupported prior class for `calibrate_cutoff_2arm()`.")
}

#' @describeIn calibrate_cutoff_2arm
#' Calibrate the posterior probability cutoff for a two-arm binary endpoint.
#' @export
calibrate_cutoff_2arm.betaMix <- function(if.prior, nf.prior, prior.t = nf.prior,
                                          target = 0.05,
                                          n.t, n,
                                          theta.t = NULL, theta = NULL,
                                          delta,
                                          method = c("SAM", "rMAP", "NP"),
                                          alternative = c("greater", "less"),
                                          margin = 0,
                                          weight_rMAP = 0.5,
                                          method.w = "LRT",
                                          prior.odds = 1,
                                          interval = c(0.5, 0.999),
                                          rel.tol = 1e-5,
                                          oc_rel.tol = 1e-8,
                                          ...) {
  calibrate_cutoff_bin_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    target = target,
    n.t = n.t,
    n = n,
    theta.t = theta.t,
    theta = theta,
    delta = delta,
    method = method,
    alternative = alternative,
    margin = margin,
    weight_rMAP = weight_rMAP,
    method.w = method.w,
    prior.odds = prior.odds,
    interval = interval,
    rel.tol = rel.tol,
    oc_rel.tol = oc_rel.tol
  )
}

#' @describeIn calibrate_cutoff_2arm
#' Calibrate the posterior probability cutoff for a two-arm continuous endpoint.
#' @param sigma.t Known sampling standard deviation in the treatment arm.
#' @param sigma Known sampling standard deviation in the control arm.
#' @param n_sd_int Half-width of the numerical integration region for each arm,
#' expressed as a multiple of the corresponding standard error.
#' @export
calibrate_cutoff_2arm.normMix <- function(if.prior, nf.prior, prior.t = nf.prior,
                                          target = 0.05,
                                          n.t, n,
                                          theta.t = NULL, theta = NULL,
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
                                          sigma.t,
                                          sigma,
                                          n_sd_int = 8,
                                          ...) {
  if (missing(sigma.t)) {
    stop("`sigma.t` must be provided for continuous endpoints.")
  }
  if (missing(sigma)) {
    stop("`sigma` must be provided for continuous endpoints.")
  }

  calibrate_cutoff_cont_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    target = target,
    n.t = n.t,
    n = n,
    sigma.t = sigma.t,
    sigma = sigma,
    theta.t = theta.t,
    theta = theta,
    delta = delta,
    method = method,
    alternative = alternative,
    margin = margin,
    weight_rMAP = weight_rMAP,
    method.w = method.w,
    prior.odds = prior.odds,
    interval = interval,
    rel.tol = rel.tol,
    oc_rel.tol = oc_rel.tol,
    n_sd_int = n_sd_int
  )
}
