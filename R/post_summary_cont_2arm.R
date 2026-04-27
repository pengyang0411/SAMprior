#' Posterior Summary for Two-Arm Comparative Trial with Continuous Endpoint
#'
#' The \code{post_summary_cont_2arm} function is designed to compute the
#' posterior summary for the treatment effect in a two-arm comparative trial
#' with a continuous endpoint under one of three borrowing strategies:
#' self-adapting mixture prior (SAM), robust MAP prior with fixed weight
#' (rMAP), or non-informative prior (NP).
#'
#' The treatment effect is defined as \eqn{\tau = \theta_t - \theta}, where
#' \eqn{\theta_t} and \eqn{\theta} denote the means in the treatment and
#' control arms, respectively. Inference is based on the posterior distribution
#' of \eqn{\tau} given the observed sample means from the two arms.
#'
#' @param ybar_t Observed sample mean in the treatment arm.
#' @param ybar Observed sample mean in the control arm.
#' @param if.prior Informative prior constructed based on historical data for
#' the control arm, represented (approximately) as a normal mixture prior.
#' @param nf.prior Non-informative prior used as the robustifying component
#' for the control arm prior.
#' @param prior.t Prior used for the treatment arm. If missing, the default
#' value is set to be \code{nf.prior}.
#' @param n.t Sample size for the treatment arm.
#' @param n Sample size for the control arm.
#' @param sigma.t Known standard deviation in the treatment arm.
#' @param sigma Known standard deviation in the control arm.
#' @param delta Clinically significant difference used for the SAM prior.
#' This argument is only used when \code{method = "SAM"}.
#' @param cutoff Posterior probability cutoff used for decision making.
#' The null hypothesis is rejected if the posterior tail probability
#' exceeds \code{cutoff}.
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
#'
#' @details
#' The posterior for the treatment arm is obtained by updating
#' \code{prior.t} using the observed sample mean \code{ybar_t}. The posterior
#' for the control arm depends on the selected borrowing strategy:
#'
#' \itemize{
#' \item \code{SAM}: the control arm prior is a mixture of \code{if.prior}
#' and \code{nf.prior}, with adaptive mixture weight determined by the
#' SAM borrowing rule.
#' \item \code{rMAP}: the control arm prior is a fixed mixture of
#' \code{if.prior} and \code{nf.prior}, with fixed weight
#' \code{weight_rMAP}.
#' \item \code{NP}: the control arm prior is \code{nf.prior} only.
#' }
#'
#' When \code{alternative = "greater"}, inference is based on
#' \eqn{P(\theta_t - \theta > margin \mid \bar{Y}_t, \bar{Y})}. When
#' \code{alternative = "less"}, inference is based on
#' \eqn{P(\theta_t - \theta < -margin \mid \bar{Y}_t, \bar{Y})}.
#'
#' The posterior mean and posterior variance of \eqn{\tau} are defined as
#' are defined as
#' \deqn{E(\tau \mid \bar{Y}_t, \bar{Y}) = E(\theta_t \mid \bar{Y}_t) - E(\theta \mid \bar{Y}),}
#' and
#' \deqn{\mathrm{Var}(\tau \mid \bar{Y}_t, \bar{Y}) = \mathrm{Var}(\theta_t \mid \bar{Y}_t) + \mathrm{Var}(\theta \mid \bar{Y}),}
#' where independence between the treatment and control arm posteriors is
#' assumed conditional on the current and historical data.
#'
#' @return A list containing the posterior probability in the requested
#' direction, posterior mean and variance of \eqn{\tau}, decision indicator,
#' borrowing weight used for the control arm prior, and the corresponding
#' trial data and method information.
#'
#' @seealso \code{\link{SAM_weight}}, \code{\link{SAM_prior}}
#'
#' @export
post_summary_cont_2arm <- function(ybar_t, ybar,
                                   if.prior, nf.prior, prior.t = nf.prior,
                                   n.t, n,
                                   sigma.t, sigma,
                                   delta,
                                   cutoff,
                                   method = c("SAM", "rMAP", "NP"),
                                   alternative = c("greater", "less"),
                                   margin = 0,
                                   weight_rMAP = 0.5,
                                   method.w = "LRT",
                                   prior.odds = 1) {
  method <- match.arg(method)
  alternative <- match.arg(alternative)

  if (!is.numeric(ybar_t) || length(ybar_t) != 1 || !is.finite(ybar_t)) {
    stop("`ybar_t` must be a finite scalar.")
  }
  if (!is.numeric(ybar) || length(ybar) != 1 || !is.finite(ybar)) {
    stop("`ybar` must be a finite scalar.")
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff <= 0 || cutoff >= 1) {
    stop("`cutoff` must be a scalar in (0, 1).")
  }
  if (!is.numeric(margin) || length(margin) != 1 || !is.finite(margin) || margin < 0) {
    stop("`margin` must be a non-negative scalar.")
  }
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
  if (!method.w %in% c("LRT", "PPR")) {
    stop("`method.w` must be either \"LRT\" or \"PPR\".")
  }
  if (!is.numeric(prior.odds) || length(prior.odds) != 1 ||
      !is.finite(prior.odds) || prior.odds <= 0) {
    stop("`prior.odds` must be a positive scalar.")
  }

  post_t <- .posterior_mix_update(
    prior = prior.t,
    ybar = ybar_t,
    n = n.t,
    sigma = sigma.t
  )

  if (method == "SAM") {
    post_c <- .control_sam_update(
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
    weight_used <- post_c$wSAM

  } else if (method == "rMAP") {
    if (is.null(weight_rMAP) || length(weight_rMAP) != 1 ||
        !is.finite(weight_rMAP) || weight_rMAP < 0 || weight_rMAP > 1) {
      stop("`weight_rMAP` must be a scalar in [0, 1] when `method = \"rMAP\"`.")
    }

    post_c <- .control_fixed_update(
      if.prior = if.prior,
      nf.prior = nf.prior,
      ybar = ybar,
      n = n,
      sigma = sigma,
      weight = weight_rMAP
    )
    weight_used <- weight_rMAP

  } else {
    post_c <- .posterior_mix_update(
      prior = nf.prior,
      ybar = ybar,
      n = n,
      sigma = sigma
    )
    weight_used <- 0
  }

  post_prob_from_margin <- function(m) {
    out <- 0
    for (j in seq_along(post_t$w)) {
      for (k in seq_along(post_c$w)) {
        z_jk <- (post_t$mu[j] - post_c$mu[k] - m) /
          sqrt(post_t$var[j] + post_c$var[k])
        out <- out + post_t$w[j] * post_c$w[k] * pnorm(z_jk)
      }
    }
    out
  }

  post_prob_greater <- post_prob_from_margin(margin)

  post_prob <- if (alternative == "greater") {
    post_prob_greater
  } else {
    1 - post_prob_from_margin(-margin)
  }

  mean_t <- sum(post_t$w * post_t$mu)
  var_t  <- sum(post_t$w * (post_t$var + post_t$mu^2)) - mean_t^2

  mean_c <- sum(post_c$w * post_c$mu)
  var_c  <- sum(post_c$w * (post_c$var + post_c$mu^2)) - mean_c^2

  post_mean <- mean_t - mean_c
  post_var  <- var_t + var_c

  list(
    post_prob = post_prob,
    post_prob_greater = post_prob_greater,
    post_mean = post_mean,
    post_var = post_var,
    decision = as.numeric(post_prob > cutoff),
    weight = weight_used,
    method = method,
    alternative = alternative,
    margin = margin,
    ybar_t = ybar_t,
    ybar = ybar
  )
}
