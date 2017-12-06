#' Calculate a Replication Bayes Factor for F-Tests.
#'
#' \code{RBF_ttest} calculates a Replication Bayes Factor for F-Tests from
#' balanced fixed effect, between subject ANOVA designs.
#'
#' The Replication Bayes Factor is a marginal likelihood ratio between two
#' positions (see Verhagen & Wagenmakers, 2014):
#' \itemize{
#'     \item \strong{H0}: The position of a skeptic, who does not place
#'         confidence in the findings of the original study and assumes
#'         \eqn{\theta \approx 0}.
#'     \item \strong{Hr}: The position of a proponent, who expects the original
#'         study to be a faithful estimation of the effect size and for which
#'         the replication offers additional evidence. This is formulated as
#'         \eqn{\theta \approx \theta_{orig}} with \eqn{\theta_{orig}} being
#'         the effect size estimate from the original study.
#' }
#'
#' @references
#' \itemize{
#'     \item \insertRef{Verhagen2014}{ReplicationBF}
#'     \item \insertRef{Harms2016}{ReplicationBF}
#' }
#'
#' @param F.orig F-statistic from the original study.
#' @param df.orig Numeric vector containing the degrees of freedom for the
#'     F-test from the original study.
#' @param N.orig Total number of observations in the original study.
#' @param F.rep F-statistic from the replication study.
#' @param df.rep Numeric vector containing the degrees of freedom for the F-Test
#'     from the replication study.
#' @param N.rep Total number of observations in the replication study.
#' @param M Number of draws from the posterior distribution to approximate the
#'     marginal likelihood.
#' @param store.samples If \emph{TRUE}, the samples of the original's posterior
#'     distribution are stored in the return object.
#'
#' @return An \code{ReplicationBF} object containing the value of the
#'     Replication Bayes Factor in \code{bayesFactor}.
#' @export
#'
#' @examples
#' RBF_Ftest(27.0, c(3, 48), 52, 3.2, c(3, 33), 37)

RBF_Ftest <- function(F.orig, df.orig, N.orig, F.rep, df.rep, N.rep,
                      M = 1e5, store.samples = FALSE)
{
  # Check arguments ------------------------------------------------------------
  if (length(df.orig) != 2 | length(df.rep) != 2)
    stop(paste("Numeric vectors with degrees of freedom has to contain exactly",
         "two values."))

  if (length(F.orig) != 1 | length(F.rep) != 1 | length(N.orig) != 1 |
      length(N.rep) != 1)
    stop("F-statistic and overall N have to be a single numeric value.")

  # Loglikelihood of the original study ----------------------------------------
  loglik <- function(theta, XF, Xdf1, Xdf2, XN) {
    ifelse(theta < 0,
           log(0),
           df(XF, Xdf1, Xdf2, ncp = theta * XN, log = T))
  }

  # MCMC approximation to the original study's posterior -----------------------
  R.utils::captureOutput({
    posterior.sample <- MCMCpack::MCMCmetrop1R(loglik, theta.init = c(0),
                                               logfun = T, burnin = 500, mcmc = M,
                                               thin = 1, verbose = 0,
                                               V = matrix(ncol = 1, nrow = 1,
                                                          c(1)),
                                               XF = F.orig, Xdf1 = df.orig[1],
                                               Xdf2 = df.orig[2], XN = N.orig)
  })

  # Calculate (marginal) likelihoods -------------------------------------------
  likelihood.h0 = df(F.rep, df.rep[1], df.rep[2])
  likelihood.hr = mean(df(F.rep, df.rep[1], df.rep[2],
                          ncp = posterior.sample * N.rep))

  rbf <- likelihood.hr / likelihood.h0

  if (!store.samples)
    posterior.sample = NULL

  # Prepare return object ------------------------------------------------------
  ret.object <- list(
    # Numeric value of the Bayes Factor
    bayesFactor = rbf,

    # Samples from the original study's posterior (empty if store.samples = F)
    posteriorSamples = posterior.sample,

    # String identifying the RBF test
    test = "Replication Bayes Factor for F-tests",

    # User's function call
    functionCall = match.call(),

    # Used approximation method (Normal or MCMC)
    approxMethod = "MCMC",

    # Details of the original study
    originalStudy = list(
      F = F.orig,
      df = df.orig,
      N = N.orig
    ),

    # Details of the replication study
    replicationStudy = list(
      F = F.orig,
      df = df.orig,
      N = N.orig
    )
  )
  class(ret.object) <- "ReplicationBF"

  return(ret.object)
}