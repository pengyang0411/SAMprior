#' Evaluate One Scenario for a Two-Arm Comparative Trial with Binary Endpoint
#'
#' The \code{eval_scenario_bin_2arm} function is designed to evaluate
#' repeated-sampling operating characteristics for a two-arm comparative trial
#' with a binary endpoint under one borrowing strategy: self-adapting mixture
#' prior (SAM), robust MAP prior with fixed weight (rMAP), or non-informative
#' prior (NP).
#'
#' The treatment effect is defined as \eqn{\tau = \theta_t - \theta},
#' where \eqn{\theta_t} and \eqn{\theta} denote the true response rates in the
#' treatment and control arms, respectively.
#'
#' For a given true scenario \eqn{(\theta_t, \theta)}, this function computes
#' the repeated-sampling rejection probability, bias, root mean squared error
#' (RMSE), and mean borrowing weight. The rejection probability is accelerated
#' by exploiting monotonicity of the posterior decision in the treatment-arm
#' response count for each fixed control-arm response count.
#'
#' @param if.prior Informative prior constructed based on historical data for
#' the control arm, represented (approximately) as a beta mixture prior.
#' @param nf.prior Non-informative prior used as the robustifying component
#' for the control arm prior.
#' @param prior.t Prior used for the treatment arm. If missing, the default
#' value is set to be \code{nf.prior}.
#' @param n.t Number of subjects in the treatment arm.
#' @param n Number of subjects in the control arm.
#' @param theta.t True treatment arm response rate.
#' @param theta True control arm response rate.
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
#' @param rel.tol Relative tolerance for numerical integration used in
#' posterior probability calculations.
#'
#' @return A one-row data frame with the following columns:
#' \describe{
#'   \item{theta}{True control arm response rate.}
#'   \item{theta.t}{True treatment arm response rate.}
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
#' @export
eval_scenario_bin_2arm <- function(if.prior, nf.prior, prior.t = nf.prior,
                                   n.t, n,
                                   theta.t, theta,
                                   cutoff,
                                   delta,
                                   method = c("SAM", "rMAP", "NP"),
                                   alternative = c("greater", "less"),
                                   margin = 0,
                                   weight_rMAP = 0.5,
                                   method.w = "LRT",
                                   prior.odds = 1,
                                   rel.tol = 1e-8) {
  method <- match.arg(method)
  alternative <- match.arg(alternative)

  if (!is.numeric(n.t) || length(n.t) != 1 || !is.finite(n.t) || n.t <= 0) {
    stop("`n.t` must be a positive scalar.")
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) {
    stop("`n` must be a positive scalar.")
  }
  if (!is.numeric(theta.t) || length(theta.t) != 1 || !is.finite(theta.t) ||
      theta.t < 0 || theta.t > 1) {
    stop("`theta.t` must be a scalar in [0, 1].")
  }
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta) ||
      theta < 0 || theta > 1) {
    stop("`theta` must be a scalar in [0, 1].")
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
  if (!is.numeric(rel.tol) || length(rel.tol) != 1 ||
      !is.finite(rel.tol) || rel.tol <= 0) {
    stop("`rel.tol` must be a positive scalar.")
  }
  if (method == "rMAP") {
    if (is.null(weight_rMAP) || length(weight_rMAP) != 1 ||
        !is.finite(weight_rMAP) || weight_rMAP < 0 || weight_rMAP > 1) {
      stop("`weight_rMAP` must be a scalar in [0, 1] when `method = \"rMAP\"`.")
    }
  }

  if (!exists(".post_prob_bin_from_posts", mode = "function")) {
    stop("Internal helper `.post_prob_bin_from_posts()` not found. Please source `post_summary_bin_2arm.R` first.")
  }

  ## ---- helper: posterior mean from control posterior object ----
  post_mean_from_control <- function(post_c) {
    mean_comp <- post_c$a / (post_c$a + post_c$b)
    sum(post_c$w * mean_comp)
  }

  ## ---- helper: borrowing weight at x ----
  weight_fun_eval <- function(x) {
    if (method == "SAM") {
      weight_fun_betamix(
        x = x,
        if.prior = if.prior,
        nf.prior = nf.prior,
        n = n,
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

  ## ---- cache treatment posteriors for all x.t ----
  post_t_list <- lapply(0:n.t, function(xt) {
    .posterior_mixbeta_update(
      prior = prior.t,
      x = xt,
      n = n.t
    )
  })

  ## ---- cache control posteriors for all x ----
  post_c_list <- lapply(0:n, function(xc) {
    if (method == "SAM") {
      .control_sam_update_bin(
        if.prior = if.prior,
        nf.prior = nf.prior,
        x = xc,
        n = n,
        delta = delta,
        weight_fun = weight_fun_betamix,
        method.w = method.w,
        prior.odds = prior.odds
      )
    } else if (method == "rMAP") {
      .control_fixed_update_bin(
        if.prior = if.prior,
        nf.prior = nf.prior,
        x = xc,
        n = n,
        weight = weight_rMAP
      )
    } else {
      .posterior_mixbeta_update(
        prior = nf.prior,
        x = xc,
        n = n
      )
    }
  })

  ## ---- posterior probability from cached posteriors ----
  post_prob_from_posts <- function(post_t, post_c) {
    if (alternative == "greater") {
      .post_prob_bin_from_posts(
        post_t = post_t,
        post_c = post_c,
        margin = margin,
        rel.tol = rel.tol
      )
    } else {
      1 - .post_prob_bin_from_posts(
        post_t = post_t,
        post_c = post_c,
        margin = -margin,
        rel.tol = rel.tol
      )
    }
  }

  ## ---- find treatment-count rejection boundary for each fixed x ----
  find_xt_boundary <- function(post_c) {
    decision_at <- function(xt) {
      post_prob_from_posts(post_t_list[[xt + 1L]], post_c) > cutoff
    }

    dec_low <- decision_at(0)
    dec_high <- decision_at(n.t)

    if (alternative == "greater") {
      if (dec_low) return(-1L)
      if (!dec_high) return(n.t + 1L)

      lo <- 0L
      hi <- n.t
      while (lo + 1L < hi) {
        mid <- (lo + hi) %/% 2L
        if (decision_at(mid)) {
          hi <- mid
        } else {
          lo <- mid
        }
      }
      hi
    } else {
      if (!dec_low) return(-1L)
      if (dec_high) return(n.t + 1L)

      lo <- 0L
      hi <- n.t
      while (lo + 1L < hi) {
        mid <- (lo + hi) %/% 2L
        if (decision_at(mid)) {
          lo <- mid
        } else {
          hi <- mid
        }
      }
      lo
    }
  }

  ## ---- rejection probability using thresholded binomial tails ----
  reject_prob <- 0

  for (x in 0:n) {
    p_x <- dbinom(x, size = n, prob = theta)
    post_c <- post_c_list[[x + 1L]]
    xt_boundary <- find_xt_boundary(post_c)

    pr_reject_given_x <- if (alternative == "greater") {
      if (xt_boundary < 0L) {
        1
      } else if (xt_boundary > n.t) {
        0
      } else {
        pbinom(xt_boundary - 1L, size = n.t, prob = theta.t, lower.tail = FALSE)
      }
    } else {
      if (xt_boundary < 0L) {
        0
      } else if (xt_boundary > n.t) {
        1
      } else {
        pbinom(xt_boundary, size = n.t, prob = theta.t, lower.tail = TRUE)
      }
    }

    reject_prob <- reject_prob + pr_reject_given_x * p_x
  }

  ## ---- exact control-arm moments for bias / RMSE / mean weight ----
  mean_post <- 0
  mean_sq_post <- 0
  mean_weight <- 0

  for (x in 0:n) {
    p_x <- dbinom(x, size = n, prob = theta)
    post_c <- post_c_list[[x + 1L]]
    post_mean_x <- post_mean_from_control(post_c)
    weight_x <- weight_fun_eval(x)

    mean_post <- mean_post + post_mean_x * p_x
    mean_sq_post <- mean_sq_post + post_mean_x^2 * p_x
    mean_weight <- mean_weight + weight_x * p_x
  }

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
