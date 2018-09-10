#' Calculate a Replication Bayes Factor for F-Tests.
#'
#' \code{RBF_Ftest} calculates a Replication Bayes Factor for F-Tests from
#' balanced fixed effect, between subject ANOVA designs (Harms, 2018).
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
#' The Bayes factor is estimated through Importance Sampling (Gamerman & Lopes,
#' 2006). The importance density is half-normal with parameters estimated from
#' the posterior distribution. Posterior distribution is sampled using
#' Metropolis-Hastings from \code{MCMCpack::MCMCmetrop1R}.
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

  # Model likelihoods, priors and posteriors -----------------------------------
  loglik <- function(theta, XF, Xdf1, Xdf2, XN) {
    ifelse(theta < 0,
           log(0),
           stats::df(XF, Xdf1, Xdf2, ncp = theta * XN, log = T))
  }

  prior.orig <- function(theta) {
    # We start out with noninformative/uniform prior - even if this is not a
    # proper probability density
    return(log(1))
  }

  posterior.orig <- function(theta, XF, Xdf1, Xdf2, XN) {
    # Non-normalized posterior density by multiplying likelihood and prior
    return(loglik(theta, XF, Xdf1, Xdf2, XN) +
             prior.orig(theta))
  }

  posterior.rep <- function(theta, XF, Xdf1, Xdf2, XN,
                            XForig, Xdf1orig, Xdf2orig, XNorig, XMLorig) {
    # Non-normalized posterior density by multiplying likelihood and
    # prior (i.e. non-normalized posterior of original study) and dividing by
    # marginal likelihood/normalizing constant from original study
    return(loglik(theta, XF, Xdf1, Xdf2, XN) +
             posterior.orig(theta, XForig, Xdf1orig, Xdf2orig, XNorig) -
             log(XMLorig))
  }

  # MCMC approximation settings ------------------------------------------------
  proposal.vcov <- matrix(ncol = 1, nrow = 1, c(1))
  sampling.tune <- 0.5
  sampling.thin <- 1

  # MCMC approximation to the original study's posterior -----------------------
  R.utils::captureOutput({
    posterior.sample.orig <- MCMCpack::MCMCmetrop1R(posterior.orig,
                                                    theta.init = c(0),
                                                    logfun = T, burnin = 500,
                                                    mcmc = M,
                                                    thin = sampling.thin,
                                                    tune = sampling.tune,
                                                    V = proposal.vcov,
                                                    verbose = 0,
                                                    # Observed data (original)
                                                    XF = F.orig,
                                                    Xdf1 = df.orig[1],
                                                    Xdf2 = df.orig[2],
                                                    XN = N.orig)
  })

  # Importance sampling estimate for marginal likelihood of original study -----
  is.mean.est <- mean(posterior.sample.orig)
  is.sd.est <- stats::sd(posterior.sample.orig)
  is.sample.orig <- truncnorm::rtruncnorm(M, mean = is.mean.est, sd = is.sd.est,
                                          a = 0)
  modelevidence.orig <- mean(exp(posterior.orig(theta = is.sample.orig,
                                                XF = F.orig, Xdf1 = df.orig[1],
                                                Xdf2 = df.orig[2],
                                                XN = N.orig)) /
                               truncnorm::dtruncnorm(is.sample.orig,
                                                     mean = is.mean.est,
                                                     sd = is.sd.est,
                                                     a = 0))

  # MCMC approximation to the replication study's posterior --------------------
  R.utils::captureOutput({
    posterior.sample.rep <- MCMCpack::MCMCmetrop1R(posterior.rep,
                                                   theta.init = c(0),
                                                   logfun = T, burnin = 500,
                                                   mcmc = M,
                                                   thin = sampling.thin,
                                                   tune = sampling.tune,
                                                   V = proposal.vcov,
                                                   verbose = 0,
                                                   # Observed data (original)
                                                   XForig = F.orig,
                                                   Xdf1orig = df.orig[1],
                                                   Xdf2orig = df.orig[2],
                                                   XNorig = N.orig,
                                                   XMLorig = modelevidence.orig,
                                                   # Observed data (replication)
                                                   XF = F.rep,
                                                   Xdf1 = df.rep[1],
                                                   Xdf2 = df.rep[2],
                                                   XN = N.rep)
  })

  # Importance sampling estimate for marginal likelihood of original study -----
  is.mean.est <- mean(posterior.sample.rep)
  is.sd.est <- stats::sd(posterior.sample.rep)
  is.sample.rep <- truncnorm::rtruncnorm(M, mean = is.mean.est, sd = is.sd.est,
                                         a = 0)
  modelevidence.rep <- mean(exp(posterior.rep(theta = is.sample.rep,
                                              XForig = F.orig,
                                              Xdf1orig = df.orig[1],
                                              Xdf2orig = df.orig[2],
                                              XNorig = N.orig,
                                              XMLorig = modelevidence.orig,
                                              XF = F.rep, Xdf1 = df.rep[1],
                                              Xdf2 = df.rep[2],
                                              XN = N.rep)) /
                              truncnorm::dtruncnorm(is.sample.rep,
                                                    mean = is.mean.est,
                                                    sd = is.sd.est,
                                                    a = 0))

  # Calculate (marginal) likelihoods -------------------------------------------
  likelihood.h0 <- stats::df(F.rep, df.rep[1], df.rep[2])
  #likelihood.hr_mc <- mean(df(F.rep, df.rep[1], df.rep[2],
  #                           ncp = posterior.sample.orig * N.rep))
  likelihood.hr_is <- modelevidence.rep


  rbf <- likelihood.hr_is / likelihood.h0

  if (!store.samples) {
    posterior.sample.orig <- NULL
    posterior.sample.rep <- NULL
  }

  # Prepare return object ------------------------------------------------------
  ret.object <- list(
    # Numeric value of the Bayes Factor
    bayesFactor = rbf,

    # Samples from the original study's posterior (empty if store.samples = F)
    posteriorSamplesOriginal = posterior.sample.orig,
    posteriorSamplesReplication = posterior.sample.rep,

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