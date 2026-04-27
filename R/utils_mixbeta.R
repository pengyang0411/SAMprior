#' Extract Components from a Beta Mixture Prior
#'
#' Internal helper function to extract mixture weights and beta shape
#' parameters from a beta mixture prior object.
#'
#' The function assumes the mixture prior stores its component matrix in
#' \code{prior[[1]]}, with row names \code{"w"}, \code{"a"}, and \code{"b"}
#' corresponding to the mixture weights and beta shape parameters.
#'
#' @param prior A beta mixture prior object, such as one created by
#' \code{\link[SAMprior]{mixbeta}} or returned by related functions in the
#' \pkg{SAMprior} workflow.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{w}{Numeric vector of normalized mixture weights.}
#'   \item{a}{Numeric vector of first beta shape parameters.}
#'   \item{b}{Numeric vector of second beta shape parameters.}
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
.extract_mixbeta <- function(prior) {
  comp <- prior[[1]]
  
  if (is.null(comp)) {
    stop("`prior[[1]]` is NULL.")
  }
  
  comp <- as.matrix(comp)
  
  rn <- rownames(comp)
  if (is.null(rn)) {
    stop("`prior[[1]]` must have row names 'w', 'a', and 'b'.")
  }
  
  w_idx <- match("w", rn)
  a_idx <- match("a", rn)
  b_idx <- match("b", rn)
  
  if (any(is.na(c(w_idx, a_idx, b_idx)))) {
    stop("Could not find rows 'w', 'a', 'b' in `prior[[1]]`.")
  }
  
  w <- as.numeric(comp[w_idx, , drop = TRUE])
  a <- as.numeric(comp[a_idx, , drop = TRUE])
  b <- as.numeric(comp[b_idx, , drop = TRUE])
  
  keep <- !(is.na(w) | is.na(a) | is.na(b))
  w <- w[keep]
  a <- a[keep]
  b <- b[keep]
  
  if (length(w) == 0) stop("No mixture components found.")
  if (any(a <= 0) || any(b <= 0)) {
    stop("All beta shape parameters must be positive.")
  }
  
  w <- w / sum(w)
  
  list(
    w = w,
    a = a,
    b = b,
    K = length(w)
  )
}

#' Update a Beta Prior Component
#'
#' Internal helper function to update a single beta prior component using
#' conjugate beta-binomial updating based on an observed number of responses.
#'
#' Suppose
#' \deqn{\theta \sim \mathrm{Beta}(a, b)}
#' and
#' \deqn{X \mid \theta \sim \mathrm{Binomial}(n, \theta).}
#' Then this function returns the posterior beta parameters and the
#' beta-binomial marginal probability of \eqn{X}.
#'
#' @param a First beta shape parameter.
#' @param b Second beta shape parameter.
#' @param x Observed number of responses.
#' @param n Number of subjects.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{a}{Posterior first beta shape parameter.}
#'   \item{b}{Posterior second beta shape parameter.}
#'   \item{marg}{Beta-binomial marginal probability of \code{x}.}
#' }
#'
#' @details
#' The marginal probability of \code{x} is
#' \deqn{
#' \Pr(X=x)=\binom{n}{x}\frac{B(a+x,b+n-x)}{B(a,b)},
#' }
#' where \eqn{B(\cdot,\cdot)} denotes the beta function.
#'
#' This is an internal utility and is not intended for direct user access.
#'
#' @noRd
.beta_update <- function(a, b, x, n) {
  a_post <- a + x
  b_post <- b + n - x
  
  marg <- choose(n, x) * beta(a_post, b_post) / beta(a, b)
  
  list(
    a = a_post,
    b = b_post,
    marg = marg
  )
}

#' Update a Beta Mixture Prior
#'
#' Internal helper function to update all components of a beta mixture prior
#' using conjugate beta-binomial updating based on an observed number of
#' responses.
#'
#' Each mixture component is updated separately, and the posterior mixture
#' weights are obtained by Bayes theorem using the component-specific
#' beta-binomial marginal probabilities.
#'
#' @param prior A beta mixture prior object.
#' @param x Observed number of responses.
#' @param n Number of subjects.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{w}{Posterior mixture weights.}
#'   \item{a}{Posterior first beta shape parameters.}
#'   \item{b}{Posterior second beta shape parameters.}
#'   \item{marg}{Beta-binomial marginal probabilities under the prior components.}
#' }
#'
#' @details
#' Let the prior be represented as a finite mixture of beta components. Each
#' component is updated using \code{.beta_update()}, and the posterior mixture
#' weights are proportional to the prior mixture weights times the
#' corresponding beta-binomial marginal probabilities.
#'
#' This is an internal utility and is not intended for direct user access.
#'
#' @noRd
.posterior_mixbeta_update <- function(prior, x, n) {
  mix <- .extract_mixbeta(prior)
  
  upd <- lapply(seq_len(mix$K), function(j) {
    .beta_update(mix$a[j], mix$b[j], x = x, n = n)
  })
  
  marg <- vapply(upd, `[[`, numeric(1), "marg")
  post_w_unnorm <- mix$w * marg
  den <- sum(post_w_unnorm)
  
  if (!is.finite(den) || den <= 0) {
    stop("Posterior weight normalization failed in `.posterior_mixbeta_update()`.")
  }
  
  post_w <- post_w_unnorm / den
  
  post_a <- vapply(upd, `[[`, numeric(1), "a")
  post_b <- vapply(upd, `[[`, numeric(1), "b")
  
  list(
    w = post_w,
    a = post_a,
    b = post_b,
    marg = marg
  )
}