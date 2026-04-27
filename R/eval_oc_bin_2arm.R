#' Evaluate Multiple Scenarios for a Two-Arm Comparative Trial with Binary Endpoint
#'
#' The \code{eval_oc_bin_2arm} function is designed to evaluate
#' repeated-sampling operating characteristics for multiple scenarios in a
#' two-arm comparative trial with a binary endpoint under one or more
#' borrowing strategies: self-adapting mixture prior (SAM), robust MAP prior
#' with fixed weight (rMAP), or non-informative prior (NP).
#'
#' For each scenario, the function computes the repeated-sampling rejection
#' probability, bias, RMSE, and mean borrowing weight using
#' \code{\link{eval_scenario_bin_2arm}}.
#'
#' @param if.prior Informative prior constructed based on historical data for
#' the control arm, represented (approximately) as a beta mixture prior.
#' @param nf.prior Non-informative prior used as the robustifying component
#' for the control arm prior.
#' @param prior.t Prior used for the treatment arm. If missing, the default
#' value is set to be \code{nf.prior}.
#' @param theta A vector of true control arm response rates.
#' @param theta.t A vector of true treatment arm response rates.
#' @param n.t Sample size for the treatment arm.
#' @param n Sample size for the control arm.
#' @param delta Clinically significant difference used for the SAM prior.
#' This argument is only used when \code{method = "SAM"}.
#' @param method Borrowing methods to evaluate. Any subset of
#' \code{c("SAM", "rMAP", "NP")}.
#' @param cutoff Posterior probability cutoff specification. Either a single
#' numeric value applied to all methods, or a named numeric vector/list with
#' method-specific cutoffs, for example
#' \code{c(SAM = 0.94, rMAP = 0.96, NP = 0.95)}.
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
#' @param rel.tol Relative tolerance passed to scenario-level numerical
#' integration used in \code{\link{post_summary_bin_2arm}}.
#'
#' @return A data frame with one row per scenario-method combination and columns:
#' \describe{
#'   \item{scenario}{Scenario index.}
#'   \item{theta}{True control arm response rate.}
#'   \item{theta.t}{True treatment arm response rate.}
#'   \item{delta_true}{True treatment effect, \eqn{\tau = \theta_t - \theta}.}
#'   \item{method}{Borrowing method used.}
#'   \item{alternative}{Direction of the posterior decision.}
#'   \item{cutoff}{Posterior probability cutoff used.}
#'   \item{margin}{Clinical margin used for inference.}
#'   \item{reject_prob}{Repeated-sampling rejection probability.}
#'   \item{bias}{Bias of the posterior mean estimator of \eqn{\theta}.}
#'   \item{rmse}{Root mean squared error of the posterior mean estimator of \eqn{\theta}.}
#'   \item{mean_weight}{Average borrowing weight under the specified method.}
#' }
#'
#' @details
#' The vectors \code{theta} and \code{theta.t} must have the same length.
#' Each pair \code{(theta[i], theta.t[i])} defines one scenario.
#'
#' The \code{cutoff} argument may be common across methods or method-specific.
#' This is useful when each borrowing method is calibrated separately before
#' operating characteristics are evaluated.
#'
#' @export
eval_oc_bin_2arm <- function(if.prior, nf.prior, prior.t = nf.prior,
                             theta, theta.t,
                             n.t, n,
                             delta,
                             method = c("SAM", "rMAP", "NP"),
                             cutoff,
                             alternative = c("greater", "less"),
                             margin = 0,
                             weight_rMAP = 0.5,
                             method.w = "LRT",
                             prior.odds = 1,
                             rel.tol = 1e-8) {
  alternative <- match.arg(alternative)
  method <- unique(method)

  valid_methods <- c("SAM", "rMAP", "NP")
  if (!all(method %in% valid_methods)) {
    stop("`method` must be a subset of c(\"SAM\", \"rMAP\", \"NP\").")
  }

  if (!is.numeric(theta) || !is.numeric(theta.t)) {
    stop("`theta` and `theta.t` must be numeric vectors.")
  }
  if (length(theta) != length(theta.t)) {
    stop("`theta` and `theta.t` must have the same length.")
  }
  if (length(theta) == 0) {
    stop("At least one scenario must be provided.")
  }
  if (any(!is.finite(theta)) || any(!is.finite(theta.t))) {
    stop("`theta` and `theta.t` must contain only finite values.")
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

  resolve_cutoff <- function(method_name, cutoff) {
    if (is.numeric(cutoff) && length(cutoff) == 1) {
      return(cutoff)
    }

    if (is.list(cutoff)) {
      cutoff <- unlist(cutoff, recursive = TRUE, use.names = TRUE)
    }

    if (!is.numeric(cutoff) || is.null(names(cutoff))) {
      stop("`cutoff` must be either a scalar or a named numeric vector/list by method.")
    }

    if (!(method_name %in% names(cutoff))) {
      stop(sprintf("Missing cutoff for method '%s'.", method_name))
    }

    out <- cutoff[[method_name]]

    if (!is.numeric(out) || length(out) != 1 || !is.finite(out)) {
      stop(sprintf("Cutoff for method '%s' must be a finite scalar.", method_name))
    }

    out
  }

  out_list <- vector("list", length(theta) * length(method))
  idx <- 1L

  for (i in seq_along(theta)) {
    for (m in method) {
      cutoff_m <- resolve_cutoff(m, cutoff)

      res <- eval_scenario_bin_2arm(
        if.prior = if.prior,
        nf.prior = nf.prior,
        prior.t = prior.t,
        n.t = n.t,
        n = n,
        theta.t = theta.t[i],
        theta = theta[i],
        cutoff = cutoff_m,
        delta = delta,
        method = m,
        alternative = alternative,
        margin = margin,
        weight_rMAP = weight_rMAP,
        method.w = method.w,
        prior.odds = prior.odds,
        rel.tol = rel.tol
      )

      res$scenario <- i
      res <- res[, c("scenario", "theta", "theta.t", "delta_true",
                     "method", "alternative", "cutoff", "margin",
                     "reject_prob", "bias", "rmse", "mean_weight")]

      out_list[[idx]] <- res
      idx <- idx + 1L
    }
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}
