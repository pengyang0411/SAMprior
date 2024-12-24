#' Calculating the Propensity Score-Integrated Informative Priors
#'
#' The \code{PS_prior} function is designed to calculate the Propensity Score-Integrated
#' (PS) informative prior constructed based on historical data.
#'
#' @param formula A two-sided \code{\link{formula}} object containing the study
#' indicator and covariates to be used in creating the distance measure used
#' in the matching. This formula will be supplied to the functions that estimate
#' the distance measure. For example, the formula should be specified as
#' \code{G ~ X1 + X2 + \dots} where \code{G} represents the name of study indicator
#'  and \code{X1} and \code{X2} are covariates.
#' @param data A data frame containing the variables named in \code{formula}
#' and possible other arguments.
#' @param outcome The variable name of the outcome.
#' @param study The variable name of the study indicator.
#' @param treat The variable name of the treatment indicator.
#' @param method The matching method to be used. The allowed methods are
#' \code{"nearest"} for nearest neighbor matching (on
#' the propensity score by default), \code{"optimal"} [\code{\link{method_optimal}}]
#' for optimal pair matching, \code{"full"} [\code{\link{method_full}}] for optimal
#' full matching, \code{"genetic"} [\code{\link{method_genetic}}] for genetic
#' matching, \code{"cem"} [\code{\link{method_cem}}] for coarsened exact matching,
#' \code{"exact"} [\code{\link{method_exact}}] for exact matching,
#' \code{"cardinality"} [\code{\link{method_cardinality}}] for cardinality and
#' template matching, and \code{"subclass"} [\code{\link{method_subclass}}] for
#' subclassification. When set to \code{"NULL"}, no matching will occur, but
#' propensity score estimation and common support restrictions will still occur
#' if requested. See the linked pages for each method for more details on what
#' these methods do, how the arguments below are used by each on, and what
#' additional arguments are allowed.
#' @param distance The distance measure to be used. Can be either the name of a
#' method of estimating propensity scores (e.g., \code{"glm"}), the name of a
#' method of computing a distance matrix from the covariates (e.g.,
#' \code{"mahalanobis"}), a vector of already-computed distance measures, or a
#' matrix of pairwise distances. See [\code{\link{distance}}] for allowable
#' options. The default is \code{"glm"} for propensity scores estimated with
#' logistic regression using \code{glm()}. Ignored for some methods; see individual
#' methods pages for information on whether and how the distance measure is
#' used.
#' @param ratio For methods that allow it, how many historical control units should be
#' matched to each current control unit in \eqn{k:1} matching. Should be a single integer
#' value. See the individual methods pages for information on whether and how
#' this argument is used. The default is 1 for 1:1 matching.
#' @param ps.method PS method utilize to calculate an informative prior
#' based on historical data. The allowed methods are \code{"Weighting"}
#' or \code{"Matching"}. The default method is \code{"Weighting"}.
#' @param trim Lower and upper bound of trimming used in \code{"Weighting"}.
#' The default is [0.1,0.9].
#'
#'
#'
#' @details This function aims to calculate informative priors using historical
#'  data by incorporating covariate information to enhance borrowing strength
#'  and address prior-data conflicts.
#'
#'  Let \eqn{G} be the study indicator, where \eqn{G = 1} indicate patient is
#'  from current control study, and \eqn{G = 0} indicate patient is from
#'  historical control study. Given the covariates data \eqn{X}, the propensity
#'  score is defined as follows,
#'  \deqn{e(X) = \Pr(G = 1 | X),}
#'  where \code{\link{distance}} allows different methods to estimate the
#'  propensity scores.
#'
#'  Calculate informative prior through PS matching is to identify a subset of
#'  historical data (\eqn{D_h^*}) that have similar PS as current control data
#'  (\eqn{D}). Various algorithms are available for PS matching, please refer
#'  to \code{method}. The informative prior can then be calculated based on
#'  the matched historical dataset.
#'
#'  Alternative, we can utilize the inverse probability of treatment weighting
#'  (IPTW) to adjust the distribution of \eqn{X} in historical data \eqn{D_h},
#'  making it similar to that in \eqn{D}. Specifically, for the \eqn{i}th
#'  subject, we assign a weight \eqn{\alpha_i} to the outcome \eqn{y_i} in
#'  \eqn{D_h} based on its PS \eqn{e(X_i)} and a fixed weight \eqn{\alpha_i = 1}
#'  to \eqn{X_i} in \eqn{D}, as follows:
#'  \deqn{\alpha_i = G_1 + (1 - G_i) \frac{e(X_i)}{1 - e(X_i)}.}
#'  To avoid extremely large weights that may compromise IPTW, symmetric
#'  trimming rule can be used to trim the tails of the PS distribution by
#'  input \code{trim} with default [0.1,0.9], that is to trim observations
#'  whose estimated PS is outside of this range.
#'
#'  To standardized \eqn{\alpha}, we compute the effective sample size (ESS),
#'  which approximately reflects the level of precision or equivalently its
#'  sample size, retained in the sample after weight as
#'  \eqn{n^{*}_h = (\sum \alpha_i)^2 / \sum{\alpha_i^2}}. The standardized weight
#'  is given by
#'  \deqn{\alpha_i^{*} = G_i + (1 - G_i)\frac{G_i}{\sum{\alpha_i} / n_h^{*}}.}
#'
#'  For binary endpoint \eqn{Y \sim Ber(\theta)}, the informative
#'  prior \eqn{\pi_1(\theta)} can be constructed as follows,
#'  \deqn{\pi_1(\theta) \propto L(\theta | D_h, \alpha^{*}) \pi_0(\theta)
#'   = Beta(a + \sum \alpha_i^{*}y_i, b + n_h^* - \sum \alpha_i^{*}y_i )\},}
#'  where \eqn{\pi_0(\theta)} is a non-informative prior, a natural choice is
#'  \eqn{Beta(a, b)}, with \eqn{a = b = 1}.
#'
#'  For continuous endpoint \eqn{Y \sim N(0, \sigma^2)}, suppose \eqn{\sigma^2}
#'  is unknown, with non-informative prior \eqn{p(\theta, \sigma^2) \propto 1/\sigma^2},
#'  \eqn{\pi_1(\theta)} follows a student-\eqn{t} distribution with degree of
#'  freedom \eqn{n_h^{*} - 1}. Given that \eqn{n_h^{*}} is moderate and large,
#'  it can be approximated by a normal distribution
#'  \eqn{N(\bar{y}^{*}, {s^{*}}^2 / n_h^{*})} with
#'  \deqn{\bar{y}^{*} = \sum \alpha_i^* y_i / \alpha_i^*, ~~ {s^{*}}^2 =
#'  \sum \alpha_i^* (y_i - \bar{y}^{*})^2 / (n_h^{*}  - 1).}
#'
#' @return Displays the informative prior calculated from historical data
#' based on the selected PS method.
#'
#' @references Zhao Y, Laird G, Chen J, Yuan Y. PS-SAM: doubly robust
#' propensity-score-integrated self-adapting mixture prior to dynamically
#' borrow information from historical data.
#'
#' @seealso \code{\link{matchit}}
#'
#' @examples
#' ## Load example data
#' data('PS_SAM_data')
#' ## Subset the data to contain historical data and current control
#' dat <- PS_SAM_data[PS_SAM_data$A == 0, ]
#' str(dat)
#'
#' ## Examples for binary endpoints
#' ## Generate the informative prior based on historical data using PS Matching
#' summary(PS_prior(formula = 'G ~ X_1 + X_2 + X_3',
#'                  data = dat, ps.method = 'Matching', method = 'nearest',
#'                  outcome = 'Y_binary', study = 'G', treat = 'A'))
#'
#' ## Generate the informative prior based on historical data using PS Weighting
#' summary(PS_prior(formula = 'G ~ X_1 + X_2 + X_3',
#'                  data = dat, ps.method = 'Weighting',
#'                  outcome = 'Y_binary', study = 'G', treat = 'A'))
#'
#' ## Examples for continuous endpoints
#' ## Generate the informative prior based on historical data using PS Matching
#' summary(PS_prior(formula = 'G ~ X_1 + X_2 + X_3',
#'                  data = dat, ps.method = 'Matching', method = 'nearest',
#'                  outcome = 'Y_continuous', study = 'G', treat = 'A'))
#'
#' ## Generate the informative prior based on historical data using PS Weighting
#' summary(PS_prior(formula = 'G ~ X_1 + X_2 + X_3',
#'                  data = dat, ps.method = 'Weighting',
#'                  outcome = 'Y_continuous', study = 'G', treat = 'A'))
#'
#'
#' @import Metrics
#' @import RBesT
#' @import assertthat
#' @import checkmate
#' @import ggplot2
#' @import stats
#' @import MatchIt
#' @export
PS_prior <- function(formula, data, outcome, study, treat, method, distance, ratio, ps.method, trim){


  if(missing(outcome)){
    stop("Plase input the name of outcome variable.")
  }else if(!outcome %in% colnames(data)){
    stop("The name of outcome variable does not match any column names in the dataset.")
  }

  if(missing(study)){
    stop("Plase input the name of study indicator.")
  }else if(!study %in% colnames(data)){
    stop("The name of study indicator does not match any column names in the dataset.")
  }

  if(missing(treat)){
    stop("Plase input the name of treatment indicator.")
  }else if(!treat %in% colnames(data)){
    stop("The name of treatment indicator does not match any column names in the dataset.")
  }

  if(missing(distance)){
    message("Using default glm method to estimate the propensity score.")
    distance <- 'glm'
  }

  if(missing(ps.method)){
    message("Using default mehtod Weighting to construct PS prior.")
    ps.method <- 'Weighting'
  }else{
    if(sum(ps.method %in% c('Matching', 'Weighting')) == 0){
      stop("Plase input the ps.method as either Weighting or Matching.")
    }
  }

  if(missing(method) & ps.method == 'Matching'){
    message("Using default nearest neighbor matching method for propensity score matching.")
    method <- 'nearest'
  }else if(missing(method)){
    method <- 'nearest'
  }

  if(missing(ratio) & ps.method == 'Matching') {
    message("Using default ratio as 1 for propensity score matching.")
    ratio <- 1
  }else if(missing(ratio)){
    ratio <- NULL
  }

  if(missing(trim) & ps.method == 'Weighting'){
    message("Using default triming range as [0.1, 0.9].")
    trim <- c(0.1, 0.9)
  }else if(missing(trim)){
    trim <- c(0.1, 0.9)
  }
  checkmate::assert_number(trim[1], lower=0, upper=1)
  checkmate::assert_number(trim[2], lower=0, upper=1)

  ## Determine the outcome data
  outcome_data <- data[ , outcome]

  # Determine outcome type
  is_binary <- all(outcome_data %in% c(0, 1))
  is_continuous <- is.numeric(outcome_data) && !is_binary

  if (is_binary) {
    # message("Binary outcome detected.")
    out <- PS_prior.beta(formula   = formula,
                         data      = data,
                         outcome   = outcome,
                         study     = study,
                         treat     = treat,
                         method    = method,
                         distance  = distance,
                         ratio     = ratio,
                         ps.method = ps.method,
                         trim      = trim)
  } else if (is_continuous) {
    # message("Continuous outcome detected.")
    out <- PS_prior.norm(formula   = formula,
                         data      = data,
                         outcome   = outcome,
                         study     = study,
                         treat     = treat,
                         method    = method,
                         distance  = distance,
                         ratio     = ratio,
                         ps.method = ps.method,
                         trim      = trim)
  } else {
    stop("Unsupported outcome type. Outcome must be either binary or continuous.")
  }

  return(out)

}

#' @describeIn PS_prior The function calculates the Propensity Score-Integrated
#' informative prior based on historical data for binary and continuous endpoint.
#' @export
PS_prior.default <- function(formula, data, outcome, study, treat, method, distance, ratio, ps.method, trim) "Unknown density"

#' @describeIn PS_prior The function calculates the Propensity Score-Integrated
#' informative prior based on historical data for binary endpoint.
#' @export
PS_prior.beta <- function(formula, data, outcome, study, treat, method, distance, ratio, ps.method, trim) {


  ## PS Matching
  if(ps.method == 'Matching'){

    ## Estimate the propensity score via MatchIt
    m.out <- MatchIt::matchit(formula = as.formula(formula), data = data, method = method,
                              distance = distance, ratio = ratio, caliper = 0.2)
    out   <- match.data(m.out, data = data)

    ## Subset the external data from matching
    dat_h <- subset(out, out[ , treat] == 0 & out[ , study] == 0)
    dat_h <- dat_h[ , outcome]

    ## Summary statsitics of subsetted historical data
    n_h <- length(dat_h)
    x_h <- sum(dat_h)

    ## Construct the informative prior
    if.prior <- RBesT::mixbeta(c(1, x_h, n_h - x_h))

  }

  ## IPW
  if(ps.method == 'Weighting'){

    ## Estimate the propensity score via MatchIt
    m.out <- MatchIt::matchit(formula = as.formula(formula), data = data,
                              distance = distance)
    out   <- match.data(m.out, data = data)

    ## Triming
    data$distance <- m.out$distance
    data_trimed <- subset(data, data$distance >= trim[1] & data$distance <= trim[2])
    ## Subset to historical data
    data_trimed <- subset(data_trimed, data_trimed[ , study] == 0 & data_trimed[ , treat] == 0)

    ## IPTW-ATT weight
    w_alpha <- data_trimed$distance / (1 - data_trimed$distance)

    # ESS
    n_h0 <- ((sum(w_alpha))^2) / (sum((w_alpha)^2))

    ## Normalized weight
    w_c  <- w_alpha / sum(w_alpha) *  n_h0

    ## Sum weighted of historical data
    x_h0 <- sum(data_trimed[ , outcome] * w_c)

    ## Construct the informative prior
    if.prior <- RBesT::mixbeta(c(1, x_h0, n_h0 - x_h0))

  }

  return(if.prior)

}

#' @describeIn PS_prior The function calculates the Propensity Score-Integrated
#' informative prior based on historical data for continuous endpoint.
#' @export
PS_prior.norm <- function(formula, data, outcome, study, treat, method, distance, ratio, ps.method, trim) {

  ## PS Matching
  if(ps.method == 'Matching'){

    ## Estimate the propensity score via MatchIt
    m.out <- MatchIt::matchit(formula = as.formula(formula), data = data, method = method,
                              distance = distance, ratio = ratio, caliper = 0.2)
    out   <- match.data(m.out, data = data)

    ## Subset the external data from matching
    dat_h <- subset(out, out[ , treat] == 0 & out[ , study] == 0)
    dat_h <- dat_h[ , outcome]

    ## Summary statsitics of subsetted historical data
    n_h <- length(dat_h)
    x_h <- mean(dat_h)

    ## Construct the informative prior
    if.prior <- RBesT::mixnorm(c(1,  x_h , sqrt(var(dat_h)/n_h)), sigma = sqrt(var(dat_h)))

  }

  ## IPW
  if(ps.method == 'Weighting'){


    ## Estimate the propensity score via MatchIt
    m.out <- MatchIt::matchit(formula = as.formula(formula), data = data,
                              distance = distance)
    out   <- match.data(m.out, data = data)

    # Triming
    data$distance <- m.out$distance
    data_trimed <- subset(data, data$distance >= trim[1] & data$distance <= trim[2])
    ## Subset to historical data
    data_trimed <- subset(data_trimed, data_trimed[ , study] == 0 & data_trimed[ , treat] == 0)

    ## Subset to historical data
    data_trimed <- subset(data_trimed, data_trimed[ , study] == 0)

    ## IPTW-ATT weight
    w_alpha <- data_trimed$distance / (1 - data_trimed$distance)

    # ESS
    n_h0 <- ((sum(w_alpha))^2) / (sum((w_alpha)^2))

    ## Normalized weight
    w_c  <- w_alpha / sum(w_alpha) *  n_h0

    ## mean of historical data
    x_h0  <- sum(data_trimed[ , outcome] * w_c) / sum(w_c)
    var_h0 <- sum((data_trimed[ , outcome] - x_h0)^2 * w_c) / (n_h0 - 1)

    ## Construct the informative prior
    if.prior <- RBesT::mixnorm(c(1,  x_h0 , sqrt(var_h0 / n_h0)), sigma = sqrt(var_h0))

  }

  return(if.prior)
}

