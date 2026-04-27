#' Extract Components from a Normal Mixture Prior
#'
#' Internal helper function to extract mixture weights, component means,
#' and component standard deviations from a normal mixture prior object.
#'
#' The function assumes the mixture prior stores its component matrix in
#' \code{prior[[1]]}, with row names \code{"w"}, \code{"m"}, and \code{"s"}
#' corresponding to the mixture weights, component means, and component
#' standard deviations, respectively.
#'
#' @param prior A normal mixture prior object, such as one created by
#' \code{\link[SAMprior]{mixnorm}} or returned by related functions in the
#' \pkg{SAMprior} workflow.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{w}{Numeric vector of normalized mixture weights.}
#'   \item{mu}{Numeric vector of component means.}
#'   \item{sd}{Numeric vector of component standard deviations.}
#'   \item{K}{Number of mixture components.}
#' }
#'
#' @details
#' Missing component entries are removed before returning the extracted
#' quantities. The mixture weights are normalized to sum to 1.
#'
#' This is an internal utility and is not intended for direct user access.
#'
#' @noRd
.extract_mixnorm <- function(prior) {
  comp <- prior[[1]]
  
  if (is.null(comp)) {
    stop("`prior[[1]]` is NULL.")
  }
  
  comp <- as.matrix(comp)
  
  rn <- rownames(comp)
  if (is.null(rn)) {
    stop("`prior[[1]]` must have row names 'w', 'm', and 's'.")
  }
  
  w_idx <- match("w", rn)
  m_idx <- match("m", rn)
  s_idx <- match("s", rn)
  
  if (any(is.na(c(w_idx, m_idx, s_idx)))) {
    stop("Could not find rows 'w', 'm', 's' in `prior[[1]]`.")
  }
  
  w  <- as.numeric(comp[w_idx, , drop = TRUE])
  mu <- as.numeric(comp[m_idx, , drop = TRUE])
  sd <- as.numeric(comp[s_idx, , drop = TRUE])
  
  keep <- !(is.na(w) | is.na(mu) | is.na(sd))
  w  <- w[keep]
  mu <- mu[keep]
  sd <- sd[keep]
  
  if (length(w) == 0) stop("No mixture components found.")
  if (any(sd <= 0)) stop("All component SDs must be positive.")
  
  w <- w / sum(w)
  
  list(
    w  = w,
    mu = mu,
    sd = sd,
    K  = length(w)
  )
}

#' Update a Normal Prior Component
#'
#' Internal helper function to update a single normal prior component using
#' conjugate normal-normal updating based on an observed sample mean.
#'
#' Suppose
#' \deqn{\theta \sim N(\mu, sd^2)}
#' and
#' \deqn{\bar{Y} \mid \theta \sim N(\theta, \sigma^2 / n).}
#' Then this function returns the posterior mean, posterior variance, and the
#' marginal density of \eqn{\bar{Y}} under the prior component.
#'
#' @param mu Prior mean of the normal component.
#' @param sd Prior standard deviation of the normal component.
#' @param ybar Observed sample mean.
#' @param n Sample size corresponding to \code{ybar}.
#' @param sigma Known sampling standard deviation.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{mean}{Posterior mean under the component.}
#'   \item{var}{Posterior variance under the component.}
#'   \item{marg}{Marginal density of \code{ybar} under the component.}
#' }
#'
#' @details
#' The marginal density of \code{ybar} is computed under the implied sampling
#' model
#' \deqn{\bar{Y} \sim N(\mu, sd^2 + \sigma^2 / n).}
#'
#' This is an internal utility and is not intended for direct user access.
#'
#' @noRd
.normal_update <- function(mu, sd, ybar, n, sigma) {
  v_prior <- sd^2
  v_lik   <- sigma^2 / n
  
  v_post <- 1 / (1 / v_prior + 1 / v_lik)
  m_post <- v_post * (mu / v_prior + ybar / v_lik)
  
  m_y <- dnorm(ybar, mean = mu, sd = sqrt(v_prior + v_lik))
  
  list(mean = m_post, var = v_post, marg = m_y)
}

#' Update a Normal Mixture Prior
#'
#' Internal helper function to update all components of a normal mixture prior
#' using conjugate normal-normal updating based on an observed sample mean.
#'
#' Each mixture component is updated separately, and the posterior mixture
#' weights are obtained by Bayes theorem using the component-specific marginal
#' densities of the observed sample mean.
#'
#' @param prior A normal mixture prior object.
#' @param ybar Observed sample mean.
#' @param n Sample size corresponding to \code{ybar}.
#' @param sigma Known sampling standard deviation.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{w}{Posterior mixture weights.}
#'   \item{mu}{Posterior means of the updated mixture components.}
#'   \item{var}{Posterior variances of the updated mixture components.}
#'   \item{marg}{Marginal densities of \code{ybar} under the prior components.}
#' }
#'
#' @details
#' Let the prior be represented as a finite mixture of normal components. Each
#' component is updated using \code{.normal_update()}, and the posterior mixture
#' weights are proportional to the prior mixture weights times the corresponding
#' marginal densities of \code{ybar}.
#'
#' This is an internal utility and is not intended for direct user access.
#'
#' @noRd
.posterior_mix_update <- function(prior, ybar, n, sigma) {
  mix <- .extract_mixnorm(prior)
  
  upd <- lapply(seq_len(mix$K), function(j) {
    .normal_update(mix$mu[j], mix$sd[j], ybar = ybar, n = n, sigma = sigma)
  })
  
  marg <- vapply(upd, `[[`, numeric(1), "marg")
  post_w_unnorm <- mix$w * marg
  den <- sum(post_w_unnorm)
  
  if (!is.finite(den) || den <= 0) {
    stop("Posterior weight normalization failed in `.posterior_mix_update()`.")
  }
  
  post_w <- post_w_unnorm / den
  
  post_mean <- vapply(upd, `[[`, numeric(1), "mean")
  post_var  <- vapply(upd, `[[`, numeric(1), "var")
  
  list(
    w = post_w,
    mu = post_mean,
    var = post_var,
    marg = marg
  )
}