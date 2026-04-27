#' Update Control Posterior Under SAM Prior
#'
#' Internal helper function to update the control arm posterior distribution
#' under a self-adapting mixture (SAM) prior for a continuous endpoint.
#'
#' The control prior is formed by mixing an informative prior and a
#' non-informative prior using an adaptive borrowing weight determined by a
#' user-supplied weight function. Each informative component and the
#' non-informative component are updated separately under conjugate
#' normal-normal updating, and the posterior mixture weights are obtained by
#' Bayes theorem.
#'
#' @param if.prior Informative normal mixture prior for the control arm,
#' typically a MAP prior constructed from historical data.
#' @param nf.prior Non-informative normal prior used as the robustifying
#' component in the SAM prior.
#' @param ybar Observed sample mean in the control arm.
#' @param n Number of subjects in the control arm.
#' @param sigma Known sampling standard deviation in the control arm.
#' @param delta Clinically significant difference used in the SAM borrowing rule.
#' @param weight_fun Function used to compute the adaptive borrowing weight for
#' the SAM prior. The function must return a value in \eqn{[0,1]}.
#' @param method.w Methods used to determine the mixture weight for SAM priors.
#' The default method is LRT (Likelihood Ratio Test), the alternative option can
#' be PPR (Posterior Probability Ratio). See \code{\link{SAM_weight}} for more
#' details.
#' @param prior.odds The prior probability of \eqn{H_0} being true compared to
#' the prior probability of \eqn{H_1} being true using PPR method. The default
#' value is 1. See \code{\link{SAM_weight}} for more details.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{w}{Posterior mixture weights for the control arm posterior.}
#'   \item{mu}{Posterior means of the updated mixture components.}
#'   \item{var}{Posterior variances of the updated mixture components.}
#'   \item{wSAM}{Adaptive borrowing weight used in the SAM prior.}
#' }
#'
#' @details
#' Let \code{if.prior} denote the informative prior \eqn{\pi_1(\theta)} and
#' \code{nf.prior} denote the non-informative prior \eqn{\pi_0(\theta)}. The
#' SAM prior for the control arm is formed as
#' \deqn{\pi_{\mathrm{SAM}}(\theta) = w \pi_1(\theta) + (1-w)\pi_0(\theta),}
#' where \eqn{w} is the adaptive borrowing weight returned by
#' \code{weight_fun}.
#'
#' Each informative mixture component and the non-informative normal component
#' are updated separately using conjugate normal-normal updating. The posterior
#' mixture weights are proportional to the prior component weights multiplied by
#' the corresponding marginal density of the observed control arm sample mean.
#'
#' This is an internal utility and is not intended for direct user access.
#'
#' @noRd
.control_sam_update <- function(if.prior, nf.prior, ybar, n, sigma, delta,
                                weight_fun, method.w = "LRT", prior.odds = 1) {
  if (!is.function(weight_fun)) {
    stop("`weight_fun` must be a function returning the SAM weight w(ybar).")
  }
  if (!is.numeric(ybar) || length(ybar) != 1 || !is.finite(ybar)) {
    stop("`ybar` must be a finite scalar.")
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    stop("`n` must be a positive scalar.")
  }
  if (!is.numeric(sigma) || length(sigma) != 1 || !is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be a positive scalar.")
  }
  if (!method.w %in% c("LRT", "PPR")) {
    stop("`method.w` must be either \"LRT\" or \"PPR\".")
  }
  if (!is.numeric(prior.odds) || length(prior.odds) != 1 || !is.finite(prior.odds) || prior.odds <= 0) {
    stop("`prior.odds` must be a positive scalar.")
  }
  
  w_sam <- weight_fun(
    ybar = ybar,
    if.prior = if.prior,
    nf.prior = nf.prior,
    n = n,
    sigma = sigma,
    delta = delta,
    method.w = method.w,
    prior.odds = prior.odds
  )
  
  if (!is.finite(w_sam) || w_sam < 0 || w_sam > 1) {
    stop("`weight_fun` must return a finite value in [0, 1].")
  }
  
  inf_mix <- .extract_mixnorm(if.prior)
  nf_mix  <- .extract_mixnorm(nf.prior)
  
  inf_upd <- lapply(seq_len(inf_mix$K), function(k) {
    .normal_update(inf_mix$mu[k], inf_mix$sd[k], ybar = ybar, n = n, sigma = sigma)
  })
  
  nf_upd <- lapply(seq_len(nf_mix$K), function(k) {
    .normal_update(nf_mix$mu[k], nf_mix$sd[k], ybar = ybar, n = n, sigma = sigma)
  })
  
  inf_marg <- vapply(inf_upd, `[[`, numeric(1), "marg")
  nf_marg  <- vapply(nf_upd,  `[[`, numeric(1), "marg")
  
  prior_w_inf <- w_sam * inf_mix$w
  prior_w_nf  <- (1 - w_sam) * nf_mix$w
  
  post_w_unnorm <- c(prior_w_inf * inf_marg, prior_w_nf * nf_marg)
  den <- sum(post_w_unnorm)
  
  if (!is.finite(den) || den <= 0) {
    stop("Posterior weight normalization failed in `.control_sam_update()`.")
  }
  
  post_w <- post_w_unnorm / den
  
  post_mu <- c(
    vapply(inf_upd, `[[`, numeric(1), "mean"),
    vapply(nf_upd,  `[[`, numeric(1), "mean")
  )
  
  post_var <- c(
    vapply(inf_upd, `[[`, numeric(1), "var"),
    vapply(nf_upd,  `[[`, numeric(1), "var")
  )
  
  list(
    w = post_w,
    mu = post_mu,
    var = post_var,
    wSAM = w_sam
  )
}
#' Update Control Posterior Under Fixed-Weight Mixture Prior
#'
#' Internal helper function to update the control arm posterior distribution
#' under a fixed-weight robust MAP prior for a continuous endpoint.
#'
#' The control prior is formed by mixing an informative prior and a
#' non-informative prior using a prespecified fixed borrowing weight. Each
#' informative component and the non-informative component are updated
#' separately under conjugate normal-normal updating, and the posterior mixture
#' weights are obtained by Bayes theorem.
#'
#' @param if.prior Informative normal mixture prior for the control arm,
#' typically a MAP prior constructed from historical data.
#' @param nf.prior Non-informative normal prior used as the robustifying
#' component in the fixed-weight prior.
#' @param ybar Observed sample mean in the control arm.
#' @param n Number of subjects in the control arm.
#' @param sigma Known sampling standard deviation in the control arm.
#' @param weight Fixed borrowing weight assigned to the informative prior.
#' Must be a scalar in \eqn{[0,1]}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{w}{Posterior mixture weights for the control arm posterior.}
#'   \item{mu}{Posterior means of the updated mixture components.}
#'   \item{var}{Posterior variances of the updated mixture components.}
#'   \item{weight}{Fixed borrowing weight used in the prior.}
#' }
#'
#' @details
#' Let \code{if.prior} denote the informative prior \eqn{\pi_1(\theta)} and
#' \code{nf.prior} denote the non-informative prior \eqn{\pi_0(\theta)}. The
#' fixed-weight prior for the control arm is formed as
#' \deqn{\pi(\theta) = w \pi_1(\theta) + (1-w)\pi_0(\theta),}
#' where \eqn{w} is the prespecified fixed borrowing weight.
#'
#' Each informative mixture component and the non-informative normal component
#' are updated separately using conjugate normal-normal updating. The posterior
#' mixture weights are proportional to the prior component weights multiplied by
#' the corresponding marginal density of the observed control arm sample mean.
#'
#' This is an internal utility and is not intended for direct user access.
#'
#' @noRd
.control_fixed_update <- function(if.prior, nf.prior, ybar, n, sigma, weight) {
  if (!is.numeric(weight) || length(weight) != 1 || !is.finite(weight) ||
      weight < 0 || weight > 1) {
    stop("`weight` must be a finite scalar in [0, 1].")
  }
  if (!is.numeric(ybar) || length(ybar) != 1 || !is.finite(ybar)) {
    stop("`ybar` must be a finite scalar.")
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    stop("`n` must be a positive scalar.")
  }
  if (!is.numeric(sigma) || length(sigma) != 1 || !is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be a positive scalar.")
  }
  
  inf_mix <- .extract_mixnorm(if.prior)
  nf_mix  <- .extract_mixnorm(nf.prior)
  
  inf_upd <- lapply(seq_len(inf_mix$K), function(k) {
    .normal_update(
      mu = inf_mix$mu[k],
      sd = inf_mix$sd[k],
      ybar = ybar,
      n = n,
      sigma = sigma
    )
  })
  
  nf_upd <- lapply(seq_len(nf_mix$K), function(k) {
    .normal_update(
      mu = nf_mix$mu[k],
      sd = nf_mix$sd[k],
      ybar = ybar,
      n = n,
      sigma = sigma
    )
  })
  
  inf_marg <- vapply(inf_upd, `[[`, numeric(1), "marg")
  nf_marg  <- vapply(nf_upd,  `[[`, numeric(1), "marg")
  
  prior_w_inf <- weight * inf_mix$w
  prior_w_nf  <- (1 - weight) * nf_mix$w
  
  post_w_unnorm <- c(prior_w_inf * inf_marg, prior_w_nf * nf_marg)
  den <- sum(post_w_unnorm)
  
  if (!is.finite(den) || den <= 0) {
    stop("Posterior weight normalization failed in `.control_fixed_update()`.")
  }
  
  post_w <- post_w_unnorm / den
  
  post_mu <- c(
    vapply(inf_upd, `[[`, numeric(1), "mean"),
    vapply(nf_upd,  `[[`, numeric(1), "mean")
  )
  
  post_var <- c(
    vapply(inf_upd, `[[`, numeric(1), "var"),
    vapply(nf_upd,  `[[`, numeric(1), "var")
  )
  
  list(
    w = post_w,
    mu = post_mu,
    var = post_var,
    weight = weight
  )
}
#' Compute SAM Borrowing Weight for Normal Mixture Prior
#'
#' Internal helper function to compute the self-adapting mixture (SAM)
#' borrowing weight for a continuous endpoint using the summary-statistic
#' version of the likelihood ratio test rule.
#'
#' The weight is computed based on the observed control-arm sample mean and the
#' posterior mean of the informative prior. The resulting weight is used to mix
#' the informative prior and the non-informative prior in the SAM prior
#' construction.
#'
#' @param ybar Observed sample mean in the control arm.
#' @param if.prior Informative normal mixture prior for the control arm.
#' @param nf.prior Non-informative prior. Included for interface consistency
#' but not used in the calculation.
#' @param n Number of subjects in the control arm.
#' @param sigma Known sampling standard deviation in the control arm.
#' @param delta Clinically significant difference used in the SAM borrowing rule.
#' @param method.w Methods used to determine the mixture weight for SAM priors.
#' The default method is LRT (Likelihood Ratio Test), the alternative option can
#' be PPR (Posterior Probability Ratio). See \code{\link{SAM_weight}} for more
#' details.
#' @param prior.odds The prior probability of \eqn{H_0} being true compared to
#' the prior probability of \eqn{H_1} being true using PPR method. The default
#' value is 1. See \code{\link{SAM_weight}} for more details.
#'
#' @return A scalar borrowing weight in \eqn{[0,1]}.
#'
#' @details
#' Let \eqn{\theta_h} denote the posterior mean estimate from the informative
#' prior. Under the default LRT rule, the SAM borrowing weight is defined as
#' \deqn{w = \frac{1}{1 + R},}
#' where
#' \deqn{R = \exp\{\max(R_1, R_2)\},}
#' with
#' \deqn{R_1 = -\frac{1}{2}\frac{n \delta (\delta - 2 \bar{Y} + 2 \theta_h)}{\sigma^2},}
#' and
#' \deqn{R_2 = -\frac{1}{2}\frac{n \delta (\delta + 2 \bar{Y} - 2 \theta_h)}{\sigma^2}.}
#'
#' Under the PPR rule, the likelihood ratio is additionally divided by
#' \code{prior.odds}. See \code{\link{SAM_weight}} for more details.
#'
#' This is an internal utility and is not intended for direct user access.
#'
#' @noRd
weight_fun_normmix <- function(ybar, if.prior, nf.prior, n, sigma, delta,
                               method.w = "LRT", prior.odds = 1) {
  if (!is.numeric(ybar) || length(ybar) != 1 || !is.finite(ybar)) {
    stop("`ybar` must be a finite scalar.")
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    stop("`n` must be a positive scalar.")
  }
  if (!is.numeric(sigma) || length(sigma) != 1 || !is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be a positive scalar.")
  }
  if (!is.numeric(delta) || length(delta) != 1 || !is.finite(delta) || delta < 0) {
    stop("`delta` must be a non-negative scalar.")
  }
  if (!method.w %in% c("LRT", "PPR")) {
    stop("`method.w` must be either \"LRT\" or \"PPR\".")
  }
  if (!is.numeric(prior.odds) || length(prior.odds) != 1 || !is.finite(prior.odds) || prior.odds <= 0) {
    stop("`prior.odds` must be a positive scalar.")
  }
  
  theta_h <- summary(if.prior)["mean"]
  
  R1 <- -0.5 * (n * delta * (delta - 2 * ybar + 2 * theta_h)) / sigma^2
  R2 <- -0.5 * (n * delta * (delta + 2 * ybar - 2 * theta_h)) / sigma^2
  
  a <- max(R1, R2)
  
  if (method.w == "PPR") {
    a <- a - log(prior.odds)
  }
  
  if (a > 700) {
    return(0)
  }
  
  1 / (1 + exp(a))
}