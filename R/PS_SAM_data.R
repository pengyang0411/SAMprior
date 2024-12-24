#' @name PS_SAM_data
#'
#' @title Simulated Data for the Construction of Propensity Score-Integrated Informative Priors
#'
#' @description  This dataset demonstrates the construction of a Propensity
#' Score-Integrated (PS) SAM prior. It simulates a two-arm randomized
#' clinical trial (RCT) with a 2:1 randomization ratio between treatment and
#' control arms, considering both binary and continuous endpoints.
#'
#'
#' @docType data
#' @format A data frame with 600 observations.
#' \itemize{
#' \item "A" is the treatment assignment (1 = treated, 0 = control).
#' \item "G" is the study indicator (1 = current, 0 = historical).
#' \item "\eqn{X_1}" is a binary covariate.
#' \item "\eqn{X_2}" is a continuous covariate.
#' \item "\eqn{X_3}" is a continuous covariate.
#' \item "\eqn{Y_{binary}}" is binary outcome.
#' \item "\eqn{Y_{continuous}}" is continuous outcome.
#' }
#'
#' @details
#' The dataset includes:
#' \itemize{
#'   \item Sample size for treatment arm: \eqn{n_t = 200}.
#'   \item Sample size for control arm: \eqn{n_c = 100}.
#'   \item Sample size for historical control study: \eqn{n_h = 300}.
#' }
#'
#' Covariates for the control arm were generated from
#'
#' \deqn{X_1 \sim \text{Ber}(0.5), ~~ X_2 \sim N(0, 1), ~~ X_3 \sim N(0.5, 1),}
#'
#' where \eqn{\text{Ber}(\cdot)} stands for Bernoulli distribution. Covariates for the
#' historical controls were generated from a mixture distribution, with half
#' were generated the same as for the control arm, while the other half were
#' drawn from
#'
#' \deqn{X_1 \sim \text{Ber}(0.8), ~~ X_2 \sim N(-0.4, 1), ~~ X_3 \sim N(-0.2, 1).}
#'
#' For the binary endpoint, \eqn{y_i} were generated from the logit model:
#'
#' \deqn{\text{logit}(\Pr(y_i = 1 | X_{1i}, X_{2i}, X_{3i}, A_i)) = -1.4 - 0.5
#' X_{1i} + X_{2i} + 2 X_{3i} + \lambda A_i,}
#'
#' where \eqn{\lambda} is the treatment effect size, and we let \eqn{\lambda = 0.9}
#' to generate a moderate treatment effect size so that they study has a reasonable
#' power.
#'
#' For the continuous endpoint, \eqn{y_i} were generated from the following
#' normal model:
#'
#' \deqn{y_i = 1.8 X_{1i} + 0.9 X_{2i} - 2 X_{3i} + \lambda A_i + \epsilon_i,}
#'
#' where we let \eqn{\lambda = 1}, and \eqn{\epsilon_i \sim N(0, 3.5^2)}.
#'
#' This dataset enables evaluation of the PS-SAM prior's performance in addressing heterogeneity
#' between the RCT control arm and historical controls.
#'
#' @keywords datasets
#'
#' @examples
#' # Load the dataset
#' data(PS_SAM_data)
#'
#' # View the structure
#' str(PS_SAM_data)
#'
"PS_SAM_data"

