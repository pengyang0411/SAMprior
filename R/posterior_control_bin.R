#' Update Control Posterior Under SAM Prior for a Binary Endpoint
#'
#' Internal helper function to update the control arm posterior distribution
#' under a self-adapting mixture (SAM) prior for a binary endpoint.
#'
#' The control prior is formed by mixing an informative prior and a
#' non-informative prior using an adaptive borrowing weight determined by a
#' user-supplied weight function. Each informative component and the
#' non-informative component are updated separately under conjugate
#' beta-binomial updating, and the posterior mixture weights are obtained by
#' Bayes theorem.
#'
#' @param if.prior Informative beta mixture prior for the control arm,
#' typically a MAP prior constructed from historical data.
#' @param nf.prior Non-informative beta prior used as the robustifying
#' component in the SAM prior.
#' @param x Observed number of responses in the control arm.
#' @param n Number of subjects in the control arm.
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
#'   \item{a}{Posterior first beta shape parameters.}
#'   \item{b}{Posterior second beta shape parameters.}
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
#' This is an internal utility and is not intended for direct user access.
#'
#' @noRd
.control_sam_update_bin <- function(if.prior, nf.prior, x, n, delta, weight_fun,
                                    method.w = "LRT", prior.odds = 1) {
  if (!is.function(weight_fun)) {
    stop("`weight_fun` must be a function returning the SAM weight w(x).")
  }
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x > n) {
    stop("`x` must be a scalar between 0 and `n`.")
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    stop("`n` must be a positive scalar.")
  }
  if (!method.w %in% c("LRT", "PPR")) {
    stop("`method.w` must be either \"LRT\" or \"PPR\".")
  }
  if (!is.numeric(prior.odds) || length(prior.odds) != 1 || !is.finite(prior.odds) || prior.odds <= 0) {
    stop("`prior.odds` must be a positive scalar.")
  }
  
  w_sam <- weight_fun(
    x = x,
    if.prior = if.prior,
    nf.prior = nf.prior,
    n = n,
    delta = delta,
    method.w = method.w,
    prior.odds = prior.odds
  )
  
  if (!is.finite(w_sam) || w_sam < 0 || w_sam > 1) {
    stop("`weight_fun` must return a finite value in [0, 1].")
  }
  
  inf_mix <- .extract_mixbeta(if.prior)
  nf_mix  <- .extract_mixbeta(nf.prior)
  
  inf_upd <- lapply(seq_len(inf_mix$K), function(k) {
    .beta_update(inf_mix$a[k], inf_mix$b[k], x = x, n = n)
  })
  
  nf_upd <- lapply(seq_len(nf_mix$K), function(k) {
    .beta_update(nf_mix$a[k], nf_mix$b[k], x = x, n = n)
  })
  
  inf_marg <- vapply(inf_upd, `[[`, numeric(1), "marg")
  nf_marg  <- vapply(nf_upd,  `[[`, numeric(1), "marg")
  
  prior_w_inf <- w_sam * inf_mix$w
  prior_w_nf  <- (1 - w_sam) * nf_mix$w
  
  post_w_unnorm <- c(prior_w_inf * inf_marg, prior_w_nf * nf_marg)
  den <- sum(post_w_unnorm)
  
  if (!is.finite(den) || den <= 0) {
    stop("Posterior weight normalization failed in `.control_sam_update_bin()`.")
  }
  
  post_w <- post_w_unnorm / den
  
  post_a <- c(
    vapply(inf_upd, `[[`, numeric(1), "a"),
    vapply(nf_upd,  `[[`, numeric(1), "a")
  )
  
  post_b <- c(
    vapply(inf_upd, `[[`, numeric(1), "b"),
    vapply(nf_upd,  `[[`, numeric(1), "b")
  )
  
  list(
    w = post_w,
    a = post_a,
    b = post_b,
    wSAM = w_sam
  )
}

#' Update Control Posterior Under Fixed-Weight Mixture Prior for a Binary Endpoint
#'
#' Internal helper function to update the control arm posterior distribution
#' under a fixed-weight robust MAP prior for a binary endpoint.
#'
#' The control prior is formed by mixing an informative prior and a
#' non-informative prior using a prespecified fixed borrowing weight. Each
#' informative component and the non-informative component are updated
#' separately under conjugate beta-binomial updating, and the posterior
#' mixture weights are obtained by Bayes theorem.
#'
#' @param if.prior Informative beta mixture prior for the control arm.
#' @param nf.prior Non-informative beta prior used as the robustifying
#' component in the fixed-weight prior.
#' @param x Observed number of responses in the control arm.
#' @param n Number of subjects in the control arm.
#' @param weight Fixed borrowing weight assigned to the informative prior.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{w}{Posterior mixture weights for the control arm posterior.}
#'   \item{a}{Posterior first beta shape parameters.}
#'   \item{b}{Posterior second beta shape parameters.}
#'   \item{weight}{Fixed borrowing weight used in the prior.}
#' }
#'
#' @noRd
.control_fixed_update_bin <- function(if.prior, nf.prior, x, n, weight) {
  if (!is.numeric(weight) || length(weight) != 1 || !is.finite(weight) ||
      weight < 0 || weight > 1) {
    stop("`weight` must be a finite scalar in [0, 1].")
  }
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x > n) {
    stop("`x` must be a scalar between 0 and `n`.")
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    stop("`n` must be a positive scalar.")
  }
  
  inf_mix <- .extract_mixbeta(if.prior)
  nf_mix  <- .extract_mixbeta(nf.prior)
  
  inf_upd <- lapply(seq_len(inf_mix$K), function(k) {
    .beta_update(inf_mix$a[k], inf_mix$b[k], x = x, n = n)
  })
  
  nf_upd <- lapply(seq_len(nf_mix$K), function(k) {
    .beta_update(nf_mix$a[k], nf_mix$b[k], x = x, n = n)
  })
  
  inf_marg <- vapply(inf_upd, `[[`, numeric(1), "marg")
  nf_marg  <- vapply(nf_upd,  `[[`, numeric(1), "marg")
  
  prior_w_inf <- weight * inf_mix$w
  prior_w_nf  <- (1 - weight) * nf_mix$w
  
  post_w_unnorm <- c(prior_w_inf * inf_marg, prior_w_nf * nf_marg)
  den <- sum(post_w_unnorm)
  
  if (!is.finite(den) || den <= 0) {
    stop("Posterior weight normalization failed in `.control_fixed_update_bin()`.")
  }
  
  post_w <- post_w_unnorm / den
  
  post_a <- c(
    vapply(inf_upd, `[[`, numeric(1), "a"),
    vapply(nf_upd,  `[[`, numeric(1), "a")
  )
  
  post_b <- c(
    vapply(inf_upd, `[[`, numeric(1), "b"),
    vapply(nf_upd,  `[[`, numeric(1), "b")
  )
  
  list(
    w = post_w,
    a = post_a,
    b = post_b,
    weight = weight
  )
}

#' Compute SAM Borrowing Weight for a Beta Mixture Prior
#'
#' Internal helper function to compute the self-adapting mixture (SAM)
#' borrowing weight for a binary endpoint using the summary-statistic version
#' of the likelihood ratio test rule.
#'
#' @param x Observed number of responses in the control arm.
#' @param if.prior Informative beta mixture prior for the control arm.
#' @param nf.prior Non-informative prior. Included for interface consistency
#' but not used in the calculation.
#' @param n Number of subjects in the control arm.
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
#' prior. The SAM borrowing weight is defined as
#' \deqn{w = \frac{1}{1 + R},}
#' where
#' \deqn{
#' R =
#' \frac{
#' \max\{\Pr(X=x \mid \theta = \min(\theta_h + \delta, 0.99)),
#'      \Pr(X=x \mid \theta = \max(\theta_h - \delta, 0.01))\}
#' }{
#' \Pr(X=x \mid \theta = \theta_h)
#' }.
#' }
#'
#' Under the PPR rule, the likelihood ratio is additionally divided by
#' \code{prior.odds}. See \code{\link{SAM_weight}} for more details.
#'
#' @noRd
weight_fun_betamix <- function(x, if.prior, nf.prior, n, delta,
                               method.w = "LRT", prior.odds = 1) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x > n) {
    stop("`x` must be a scalar between 0 and `n`.")
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    stop("`n` must be a positive scalar.")
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
  
  p1 <- min(theta_h + delta, 0.99)
  p2 <- max(theta_h - delta, 0.01)
  
  num <- max(
    dbinom(x = x, size = n, prob = p1),
    dbinom(x = x, size = n, prob = p2)
  )
  den <- dbinom(x = x, size = n, prob = theta_h)
  
  if (!is.finite(den) || den <= 0) {
    stop("SAM weight calculation failed in `weight_fun_betamix()`.")
  }
  
  R <- num / den
  
  if (method.w == "PPR") {
    R <- R / prior.odds
  }
  
  1 / (1 + R)
}