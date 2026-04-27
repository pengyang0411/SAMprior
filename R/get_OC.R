#' Generating Operating Characteristics of SAM Priors Using Analytical Engines
#'
#' The \code{get_OC} function generates operating characteristics of SAM
#' priors using the analytical operating characteristic engines for two-arm
#' trials with binary or continuous endpoints. As an option, the operating
#' characteristics of robust MAP priors can also be generated for comparison.
#'
#' Compared with the original simulation-based implementation, this function
#' does not rely on trial simulation. Instead, it first calibrates the
#' posterior probability cutoff for each borrowing method to achieve the
#' target type I error, and then evaluates operating characteristics
#' analytically across the requested scenarios.
#'
#' @param if.prior Informative prior constructed from historical data,
#' represented (approximately) as a mixture of conjugate distributions.
#' @param theta.h Estimate of the treatment effect based on historical data.
#' Included for interface compatibility. If missing, the default value is set
#' to the posterior mean estimate from \code{if.prior}.
#' @param method.w Methods used to determine the mixture weight for SAM priors.
#' The default method is LRT (Likelihood Ratio Test), the alternative option can
#' be PPR (Posterior Probability Ratio). See \code{\link{SAM_weight}} for more
#' details.
#' @param prior.odds The prior probability of \eqn{H_0} being true compared to
#' the prior probability of \eqn{H_1} being true using PPR method. The default
#' value is 1. See \code{\link{SAM_weight}} for more details.
#' @param nf.prior Non-informative prior used as the robustifying component
#' for the control arm prior.
#' @param prior.t Prior used for the treatment arm. If missing, the default
#' value is set to be \code{nf.prior}.
#' @param delta Clinically significant difference used for the SAM prior.
#' @param n Sample size for the control arm.
#' @param n.t Sample size for the treatment arm.
#' @param target Target type I error used to calibrate the posterior probability
#' cutoff for each method. The default value is typically 0.05.
#' @param if.rMAP Whether to evaluate the operating characteristics of the
#' robust MAP prior for comparison. The default value is \code{FALSE}.
#' @param weight.rMAP Weight assigned to the informative prior component
#' (\eqn{0 \le} \code{weight.rMAP} \eqn{\le 1}) for the robust MAP prior.
#' @param theta A vector of the response rate (binary endpoints) or mean
#' (continuous endpoints) for the control arm.
#' @param theta.t A vector of the response rate (binary endpoints) or mean
#' (continuous endpoints) for the treatment arm.
#' @param alternative Direction of the posterior decision. Must be one of
#' \code{"greater"} (for superiority) or \code{"less"} (for inferiority).
#' @param margin Clinical margin. Must be a non-negative scalar. The default
#' value is \code{0}.
#' @param rel.tol Tolerance passed to numerical root finding.
#' @param oc_rel.tol Relative tolerance passed to operating characteristic
#' evaluation.
#' @param interval Search interval for the posterior probability cutoff.
#' @param n_sd_int Half-width of the numerical integration region for each arm,
#' expressed as a multiple of the corresponding standard error. Used for
#' continuous endpoints only.
#' @param ... Additional parameters. For continuous endpoints, this includes
#' \code{sigma}.
#'
#' @details
#' For each borrowing method, the function first calibrates the posterior
#' probability cutoff so that the repeated-sampling rejection probability under
#' the boundary null scenario equals the target type I error \code{target}.
#' Specifically, calibration is based on the first value of \code{theta}. Let
#' \eqn{\theta = \code{theta}[1]} denote the control-arm parameter under the
#' calibration scenario. Then the treatment-arm parameter is set to
#' \eqn{\theta_t = \theta + margin} when \code{alternative = "greater"}, and
#' to \eqn{\theta_t = \theta - margin} when \code{alternative = "less"}.
#' Thus, when \code{margin = 0}, calibration is performed under the null
#' scenario \eqn{\theta_t = \theta}, corresponding to no treatment effect
#' difference between the treatment and control arms. After the cutoff is
#' calibrated, the function evaluates the operating characteristics across all
#' requested scenarios using the analytical two-arm engines.
#'
#' @return A data frame with one row per scenario-method combination. The
#' columns are:
#' \describe{
#'   \item{Scenarios}{Scenario index.}
#'   \item{theta}{True control-arm parameter value.}
#'   \item{theta.t}{True treatment-arm parameter value.}
#'   \item{Methods}{Borrowing method, one of \code{"NP"}, \code{"rMAP"}, or \code{"SAM"}.}
#'   \item{Cutoffs}{Calibrated posterior probability cutoff for the method.}
#'   \item{Bias of theta}{Bias of the posterior mean estimator of \eqn{\theta}.}
#'   \item{RMSE of theta}{Root mean squared error of the posterior mean estimator of \eqn{\theta}.}
#'   \item{Weight}{Average borrowing weight under the method.}
#'   \item{Probability of Rejection}{Repeated-sampling rejection probability.}
#' }
#'
#' @references
#' Yang P, Zhao Y, Nie L, Vallejo J, Yuan Y.
#' SAM: Self-adapting mixture prior to dynamically borrow information from
#' historical data in clinical trials.
#' \emph{Biometrics} 2023;79(4):2857-2868.
#'
#' Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D,
#' Neuenschwander B.
#' Robust meta-analytic-predictive priors in clinical trials with historical
#' control information.
#' \emph{Biometrics} 2014;70(4):1023-1032.
#'
#' @seealso \code{\link{eval_oc_bin_2arm}},
#' \code{\link{eval_oc_cont_2arm}}, \code{\link{calibrate_cutoff_bin_2arm}},
#' \code{\link{calibrate_cutoff_cont_2arm}}
#'
#' @examples
#' ## Example: operating characteristics for a two-arm binary trial
#' ## using a SAM prior without rMAP comparison
#'
#' ## Informative prior constructed from historical data
#' if.prior <- mixbeta(c(1, 20, 40))
#'
#' ## Evaluate operating characteristics
#' OC <- get_OC(
#'   if.prior = if.prior,
#'   nf.prior = mixbeta(c(1, 1, 1)),   ## Non-informative prior for control mixture
#'   prior.t = mixbeta(c(1, 1, 1)),    ## Prior for treatment arm
#'   delta = 0.2,                      ## Clinically significant difference for SAM
#'   n = 50,                           ## Sample size for control arm
#'   n.t = 100,                        ## Sample size for treatment arm
#'   target = 0.05,                    ## Target type I error for cutoff calibration
#'   if.rMAP = FALSE,                  ## Do not include rMAP comparison
#'   theta = c(summary(if.prior)["mean"], summary(if.prior)["mean"]),
#'   theta.t = c(summary(if.prior)["mean"], 0.50),
#'   alternative = "greater",          ## Superiority test
#'   margin = 0                        ## Clinical margin
#' )
#'
#' OC
#'
#' @export
get_OC <- function(if.prior, theta.h, method.w, prior.odds, nf.prior,
                   prior.t = nf.prior,
                   delta, n, n.t, target = 0.05,
                   if.rMAP = FALSE, weight.rMAP = 0.5,
                   theta, theta.t,
                   alternative = c("greater", "less"),
                   margin = 0,
                   rel.tol = 1e-5,
                   oc_rel.tol = 1e-6,
                   interval = c(0.5, 0.999),
                   n_sd_int = 8, ...) {

  UseMethod("get_OC")
}

#' @export
get_OC.default <- function(if.prior, theta.h, method.w, prior.odds, nf.prior,
                           prior.t = nf.prior,
                           delta, n, n.t, target = 0.05,
                           if.rMAP = FALSE, weight.rMAP = 0.5,
                           theta, theta.t,
                           alternative = c("greater", "less"),
                           margin = 0,
                           rel.tol = 1e-5,
                           oc_rel.tol = 1e-6,
                           interval = c(0.5, 0.999),
                           n_sd_int = 8, ...) {

  "Unknown density"
}

.combine_oc_results_long_format <- function(res_list, cal_list) {
  preferred_order <- c("NP", "rMAP", "SAM")
  method_names <- preferred_order[preferred_order %in% names(res_list)]

  n_scen <- nrow(res_list[[method_names[1]]])
  out_list <- vector("list", n_scen * length(method_names))
  idx <- 1L

  for (s in seq_len(n_scen)) {
    for (m in method_names) {
      res <- res_list[[m]]
      cal <- cal_list[[m]]

      out_list[[idx]] <- data.frame(
        Scenarios = s,
        theta = res$theta[s],
        theta.t = res$theta.t[s],
        Methods = m,
        Cutoffs = cal$cutoff,
        `Bias of theta` = res$bias[s],
        `RMSE of theta` = res$rmse[s],
        Weight = res$mean_weight[s],
        `Probability of Rejection` = res$reject_prob[s]
      )

      idx <- idx + 1L
    }
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL

  num_cols <- sapply(out, is.numeric)
  num_cols["Scenarios"] <- FALSE
  out[num_cols] <- lapply(out[num_cols], function(x) round(x, 4))

  out
}

#' @export
get_OC.betaMix <- function(if.prior, theta.h, method.w, prior.odds, nf.prior,
                           prior.t = nf.prior,
                           delta, n, n.t, target = 0.05,
                           if.rMAP = FALSE, weight.rMAP = 0.5,
                           theta, theta.t,
                           alternative = c("greater", "less"),
                           margin = 0,
                           rel.tol = 1e-5,
                           oc_rel.tol = 1e-6,
                           interval = c(0.5, 0.999),
                           n_sd_int = 8, ...) {

  alternative <- match.arg(alternative)

  if (length(theta) != length(theta.t)) {
    stop("Theta under control and treatment should be the same length!")
  }
  if (missing(if.prior)) {
    stop("Please input the informative prior!")
  }
  if (missing(n)) {
    stop("Please input the sample size for control arm!")
  }
  if (missing(n.t)) {
    stop("Please input the sample size for treatment arm!")
  }
  if (missing(delta)) {
    stop("Please input clinically significant difference!")
  }
  if (missing(theta.h)) {
    theta.h <- summary(if.prior)["mean"]
  }
  if (missing(method.w)) {
    method.w <- "LRT"
  }
  if (missing(prior.odds)) {
    prior.odds <- 1
  }
  if (missing(nf.prior)) {
    nf.prior <- mixbeta(nf.prior = c(1, 1, 1))
  }
  if (!is.numeric(margin) || length(margin) != 1 || !is.finite(margin) || margin < 0) {
    stop("`margin` must be a non-negative scalar.")
  }

  theta_cal <- theta[1]
  theta.t_cal <- if (alternative == "greater") {
    theta_cal + margin
  } else {
    theta_cal - margin
  }

  if (theta.t_cal < 0 || theta.t_cal > 1) {
    stop("Calibration scenario for binary endpoint is outside [0, 1].")
  }

  cal_list <- list()
  res_list <- list()

  message("Calibrating the posterior probability cutoff for NP prior method ...")
  cal_list$NP <- calibrate_cutoff_bin_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    target = target,
    n.t = n.t,
    n = n,
    theta.t = theta.t_cal,
    theta = theta_cal,
    delta = delta,
    method = "NP",
    alternative = alternative,
    margin = margin,
    weight_rMAP = NULL,
    method.w = method.w,
    prior.odds = prior.odds,
    interval = interval,
    rel.tol = rel.tol,
    oc_rel.tol = oc_rel.tol
  )

  if (if.rMAP) {
    message("Calibrating the posterior probability cutoff for rMAP prior method ...")
    cal_list$rMAP <- calibrate_cutoff_bin_2arm(
      if.prior = if.prior,
      nf.prior = nf.prior,
      prior.t = prior.t,
      target = target,
      n.t = n.t,
      n = n,
      theta.t = theta.t_cal,
      theta = theta_cal,
      delta = delta,
      method = "rMAP",
      alternative = alternative,
      margin = margin,
      weight_rMAP = weight.rMAP,
      method.w = method.w,
      prior.odds = prior.odds,
      interval = interval,
      rel.tol = rel.tol,
      oc_rel.tol = oc_rel.tol
    )
  }

  message("Calibrating the posterior probability cutoff for SAM prior method ...")
  cal_list$SAM <- calibrate_cutoff_bin_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    target = target,
    n.t = n.t,
    n = n,
    theta.t = theta.t_cal,
    theta = theta_cal,
    delta = delta,
    method = "SAM",
    alternative = alternative,
    margin = margin,
    weight_rMAP = NULL,
    method.w = method.w,
    prior.odds = prior.odds,
    interval = interval,
    rel.tol = rel.tol,
    oc_rel.tol = oc_rel.tol
  )

  message("Evaluating operating characteristics for NP prior method ...")
  res_list$NP <- eval_oc_bin_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    theta = theta,
    theta.t = theta.t,
    n.t = n.t,
    n = n,
    delta = delta,
    method = "NP",
    cutoff = cal_list$NP$cutoff,
    alternative = alternative,
    margin = margin,
    weight_rMAP = NULL,
    method.w = method.w,
    prior.odds = prior.odds,
    rel.tol = oc_rel.tol
  )

  if (if.rMAP) {
    message("Evaluating operating characteristics for rMAP prior method ...")
    res_list$rMAP <- eval_oc_bin_2arm(
      if.prior = if.prior,
      nf.prior = nf.prior,
      prior.t = prior.t,
      theta = theta,
      theta.t = theta.t,
      n.t = n.t,
      n = n,
      delta = delta,
      method = "rMAP",
      cutoff = cal_list$rMAP$cutoff,
      alternative = alternative,
      margin = margin,
      weight_rMAP = weight.rMAP,
      method.w = method.w,
      prior.odds = prior.odds,
      rel.tol = oc_rel.tol
    )
  }

  message("Evaluating operating characteristics for SAM prior method ...")
  res_list$SAM <- eval_oc_bin_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    theta = theta,
    theta.t = theta.t,
    n.t = n.t,
    n = n,
    delta = delta,
    method = "SAM",
    cutoff = cal_list$SAM$cutoff,
    alternative = alternative,
    margin = margin,
    weight_rMAP = NULL,
    method.w = method.w,
    prior.odds = prior.odds,
    rel.tol = oc_rel.tol
  )

  .combine_oc_results_long_format(
    res_list = res_list,
    cal_list = cal_list
  )
}

#' @export
get_OC.normMix <- function(if.prior, theta.h, method.w, prior.odds, nf.prior,
                           prior.t = nf.prior,
                           delta, n, n.t, target = 0.05,
                           if.rMAP = FALSE, weight.rMAP = 0.5,
                           theta, theta.t,
                           alternative = c("greater", "less"),
                           margin = 0,
                           rel.tol = 1e-5,
                           oc_rel.tol = 1e-6,
                           interval = c(0.5, 0.999),
                           n_sd_int = 8, ..., sigma) {
  alternative <- match.arg(alternative)

  if (length(theta) != length(theta.t)) {
    stop("Theta under control and treatment should be the same length!")
  }
  if (missing(if.prior)) {
    stop("Please input the informative prior!")
  }
  if (missing(n)) {
    stop("Please input the sample size for control arm!")
  }
  if (missing(n.t)) {
    stop("Please input the sample size for treatment arm!")
  }
  if (missing(delta)) {
    stop("Please input clinically significant difference!")
  }
  if (missing(theta.h)) {
    theta.h <- summary(if.prior)["mean"]
  }
  if (missing(method.w)) {
    method.w <- "LRT"
  }
  if (missing(prior.odds)) {
    prior.odds <- 1
  }
  if (missing(sigma)) {
    sigma <- RBesT::sigma(if.prior)
  }
  if (missing(nf.prior)) {
    nf.prior <- mixnorm(
      nf.prior = c(1, summary(if.prior)["mean"], sigma),
      param = "ms"
    )
  }
  if (!is.numeric(margin) || length(margin) != 1 || !is.finite(margin) || margin < 0) {
    stop("`margin` must be a non-negative scalar.")
  }
  theta_cal <- theta[1]
  theta.t_cal <- if (alternative == "greater") {
    theta_cal + margin
  } else {
    theta_cal - margin
  }

  cal_list <- list()
  res_list <- list()

  message("Calibrating the posterior probability cutoff for NP prior method ...")
  cal_list$NP <- calibrate_cutoff_cont_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    target = target,
    n.t = n.t,
    n = n,
    sigma.t = sigma,
    sigma = sigma,
    theta.t = theta.t_cal,
    theta = theta_cal,
    delta = delta,
    method = "NP",
    alternative = alternative,
    margin = margin,
    weight_rMAP = NULL,
    method.w = method.w,
    prior.odds = prior.odds,
    interval = interval,
    rel.tol = rel.tol,
    oc_rel.tol = oc_rel.tol,
    n_sd_int = n_sd_int
  )

  if (if.rMAP) {
    message("Calibrating the posterior probability cutoff for rMAP prior method ...")
    cal_list$rMAP <- calibrate_cutoff_cont_2arm(
      if.prior = if.prior,
      nf.prior = nf.prior,
      prior.t = prior.t,
      target = target,
      n.t = n.t,
      n = n,
      sigma.t = sigma,
      sigma = sigma,
      theta.t = theta.t_cal,
      theta = theta_cal,
      delta = delta,
      method = "rMAP",
      alternative = alternative,
      margin = margin,
      weight_rMAP = weight.rMAP,
      method.w = method.w,
      prior.odds = prior.odds,
      interval = interval,
      rel.tol = rel.tol,
      oc_rel.tol = oc_rel.tol,
      n_sd_int = n_sd_int
    )
  }

  message("Calibrating the posterior probability cutoff for SAM prior method ...")
  cal_list$SAM <- calibrate_cutoff_cont_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    target = target,
    n.t = n.t,
    n = n,
    sigma.t = sigma,
    sigma = sigma,
    theta.t = theta.t_cal,
    theta = theta_cal,
    delta = delta,
    method = "SAM",
    alternative = alternative,
    margin = margin,
    weight_rMAP = NULL,
    method.w = method.w,
    prior.odds = prior.odds,
    interval = interval,
    rel.tol = rel.tol,
    oc_rel.tol = oc_rel.tol,
    n_sd_int = n_sd_int
  )

  message("Evaluating operating characteristics for NP prior method ...")
  res_list$NP <- eval_oc_cont_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    theta = theta,
    theta.t = theta.t,
    n.t = n.t,
    n = n,
    sigma.t = sigma,
    sigma = sigma,
    delta = delta,
    method = "NP",
    cutoff = cal_list$NP$cutoff,
    alternative = alternative,
    margin = margin,
    weight_rMAP = NULL,
    method.w = method.w,
    prior.odds = prior.odds,
    rel.tol = oc_rel.tol,
    n_sd_int = n_sd_int
  )

  if (if.rMAP) {
    message("Evaluating operating characteristics for rMAP prior method ...")
    res_list$rMAP <- eval_oc_cont_2arm(
      if.prior = if.prior,
      nf.prior = nf.prior,
      prior.t = prior.t,
      theta = theta,
      theta.t = theta.t,
      n.t = n.t,
      n = n,
      sigma.t = sigma,
      sigma = sigma,
      delta = delta,
      method = "rMAP",
      cutoff = cal_list$rMAP$cutoff,
      alternative = alternative,
      margin = margin,
      weight_rMAP = weight.rMAP,
      method.w = method.w,
      prior.odds = prior.odds,
      rel.tol = oc_rel.tol,
      n_sd_int = n_sd_int
    )
  }

  message("Evaluating operating characteristics for SAM prior method ...")
  res_list$SAM <- eval_oc_cont_2arm(
    if.prior = if.prior,
    nf.prior = nf.prior,
    prior.t = prior.t,
    theta = theta,
    theta.t = theta.t,
    n.t = n.t,
    n = n,
    sigma.t = sigma,
    sigma = sigma,
    delta = delta,
    method = "SAM",
    cutoff = cal_list$SAM$cutoff,
    alternative = alternative,
    margin = margin,
    weight_rMAP = NULL,
    method.w = method.w,
    prior.odds = prior.odds,
    rel.tol = oc_rel.tol,
    n_sd_int = n_sd_int
  )

  .combine_oc_results_long_format(
    res_list = res_list,
    cal_list = cal_list
  )
}
