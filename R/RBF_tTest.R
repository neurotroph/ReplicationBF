#' Calculate a Replication Bayes Factor for a two-sample t-Test.
#'
#' \code{RBF_ttest} calculates a Replication Bayes Factor (Verhagen &
#' Wagenmakers, 2014) for t-Tests.
#'
#' The Replication Bayes Factor is a marginal likelihood ratio between two
#' models, characterized by the following positions:
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
#' In contrast to the originally proposed Replication Bayes factor, this method
#' estimates the Bayes factor using importance sampling. This yields different
#' results than the code provided by Verhagen and Wagenmakers (2014), but is
#' more stable and less biased than the Monte Carlo estimate (see Bos, 2002) in
#' cases where original and replication study yield different effect size
#' estimates.
#'
#' @references
#' \itemize{
#'     \item \insertRef{Verhagen2014}{ReplicationBF}
#'     \item \insertRef{Bos2002}{ReplicationBF}
#' }
#'
#' @param t.orig t-statistic of the original study.
#' @param n.orig Numeric vector (2 elements) with cell sizes.
#' @param t.rep t-tstatistic of the replication study.
#' @param n.rep Numeric vector (2 elements) with cell sizes.
#' @param method Either \code{"NormApprox"} for a Normal approximation of the
#'     original study's posterior distribution or \code{"MCMC"} for an
#'     approximation through the Metropolis algorithm.
#' @param M Number of posterior samples (either MCMC or Normal distribution).
#' @param store.samples If \code{TRUE}, the returned object contains \code{M}
#'     samples from the original study's posterior distribution (required for
#'     \code{plot_RBF}).
#' @return An object of \code{ReplicationBF} class.
#'
#' @export
#' @examples
#' \dontrun{
#' # Example 1 from Verhagen & Wagenmakers (2014)
#' # Using a Normal approximation to the original's posterior distribution
#' RBF_ttest(2.18, c(10, 11), 3.06, c(27, 27), method = "NormApprox")
#' RBF_ttest(2.18, c(10, 11), 0.25, c(27, 27), method = "NormApprox")
#' RBF_ttest(2.18, c(10, 11), 2.44, c(16, 17), method = "NormApprox")
#'
#' # Using MCMC to draw samples from the original's posterior distribution
#' RBF_ttest(2.10, c(11, 11), 3.06, c(27, 28), method = "MCMC")
#' RBF_ttest(2.18, c(10, 11), 0.25, c(27, 27), method = "MCMC")
#' RBF_ttest(2.18, c(10, 11), 2.44, c(16, 17), method = "MCMC")
#' }


RBF_ttest <- function(t.orig, n.orig, t.rep, n.rep,
                      method = "NormApprox", M = 1e5, store.samples = FALSE) {
  # Check arguments ------------------------------------------------------------
  if (length(n.orig) != length(n.rep))
    stop("N from original and replication study have different lengths.")
  if (length(n.orig) > 2)
    stop("t-Test is only valid for one or two samples.")
  if (!(method == "NormApprox" | method == "MCMC"))
    stop("Approximation method not supported. Only 'NormApprox' or 'MCMC'.")

  # Calculate degrees of freedom and auxilliary variables ----------------------
  df.orig <- sum(n.orig) - length(n.orig)
  df.rep <- sum(n.rep) - length(n.rep)
  if (length(n.orig) == 1) {
    # One-sample t-Test --------------------------------------------------------
    sqrt.n.orig <- sqrt(n.orig[1])
    sqrt.n.rep <- sqrt(n.rep[1])
  } else if (length(n.orig) == 2) {
    # Two-sample t-Test --------------------------------------------------------
    sqrt.n.orig <- sqrt(1 / (1 / n.orig[1] + 1 / n.orig[2]) )
    sqrt.n.rep <- sqrt(1 / (1 / n.rep[1] + 1 / n.rep[2]) )
  }

  # Model functions for MCMC sampling ------------------------------------------
  loglik <- function(theta, Xt, Xdf, XN) {
    return(stats::dt(Xt, Xdf, ncp = theta*XN, log = T))
  }

  prior.orig <- function(theta) {
    # Using an improper, uniform prior - doesn't matter, because integral stays
    # the same.
    return(log(1))
  }

  posterior.orig <- function(theta, Xt, Xdf, XN) {
    # Non-normalized posterior density by multiplying likelihood and prior
    return(loglik(theta, Xt, Xdf, XN) +
             prior.orig(theta))
  }

  posterior.rep <- function(theta, Xt, Xdf, XN,
                            Xtorig, Xdforig, XNorig, XMLorig) {
    # Non-normalized posterior density by multiplying likelihood and
    # prior (i.e. non-normalized posterior of original study) and dividing by
    # marginal likelihood/normalizing constant from original study
    return(loglik(theta, Xt, Xdf, XN) +
             posterior.orig(theta, Xtorig, Xdforig, XNorig) -
             log(XMLorig))
  }

  # MCMC approximation settings ------------------------------------------------
  proposal.vcov <- matrix(ncol = 1, nrow = 1, c(1))
  sampling.tune <- 0.5
  sampling.thin <- 1


  # Sampling from original study's posterior -----------------------------------
  if (method == "NormApprox") {
    # Calculate Confidence Interval for NCP of noncentral t-Distribution -------
    if (!requireNamespace("MBESS", quietly = TRUE)) {
      stop("To calculate Normal approximation, please install package 'MBESS'.",
           call. = FALSE)
    }
    delta.lower <- MBESS::conf.limits.nct(t.orig, df.orig)$Lower.Limit

    # Calculate paramters for approximated Normal distribution -----------------
    mu.delta <- t.orig / sqrt.n.orig
    sd.delta <- abs( ((t.orig - delta.lower) / stats::qnorm(.025)) / sqrt.n.orig )
    posterior.sample.orig <- stats::rnorm(M, mu.delta, sd.delta)
  } else if (method == "MCMC") {
    # Approximate posterior of original study using MCMC -----------------------
    R.utils::captureOutput({
      posterior.sample.orig <- MCMCpack::MCMCmetrop1R(posterior.orig,
                                                      theta.init = stats::runif(1),
                                                      logfun = TRUE, mcmc = M,
                                                      burnin = 500, verbose = 0,
                                                      thin = sampling.thin,
                                                      tune = sampling.tune,
                                                      V = proposal.vcov,
                                                      Xt = t.orig,
                                                      Xdf = df.orig,
                                                      XN = sqrt.n.orig)
    })

  } else {
    stop("Approximation method not supported.")
  }

  # Importance sampling estimate for marginal likelihood of original study -----
  is.mean.est <- mean(posterior.sample.orig)
  is.sd.est <- stats::sd(posterior.sample.orig)
  is.sample.orig <- stats::rnorm(M, mean = is.mean.est, sd = is.sd.est)
  modelevidence.orig <- mean(exp(posterior.orig(theta = is.sample.orig,
                                                Xt = t.orig, Xdf = df.orig,
                                                XN = sqrt.n.orig)) /
                               stats::dnorm(is.sample.orig, mean = is.mean.est,
                                     sd = is.sd.est))

  # Sample posterior of replication study --------------------------------------
  R.utils::captureOutput({
    posterior.sample.rep <- MCMCpack::MCMCmetrop1R(posterior.rep,
                                                   theta.init = stats::runif(1),
                                                   logfun = TRUE, mcmc = M,
                                                   burnin = 500, verbose = 0,
                                                   thin = sampling.thin,
                                                   tune = sampling.tune,
                                                   V = proposal.vcov,
                                                   # Original data
                                                   Xtorig = t.orig,
                                                   Xdforig = df.orig,
                                                   XNorig = sqrt.n.orig,
                                                   XMLorig = modelevidence.orig,
                                                   # Replication data
                                                   Xt = t.rep,
                                                   Xdf = df.rep,
                                                   XN = sqrt.n.rep)
  })

  # Importance sampling estimate for marginal likelihood of replication --------
  is.mean.est <- mean(posterior.sample.rep)
  is.sd.est <- stats::sd(posterior.sample.rep)
  is.sample.rep <- stats::rnorm(M, mean = is.mean.est, sd = is.sd.est)
  modelevidence.rep <- mean(exp(posterior.rep(theta = is.sample.rep,
                                              # Replication study data
                                              Xt = t.rep, Xdf = df.rep,
                                              XN = sqrt.n.rep,
                                              # Original study data
                                              Xtorig = t.orig, Xdforig = df.orig,
                                              XNorig = sqrt.n.orig,
                                              XMLorig = modelevidence.orig)) /
                              stats::dnorm(is.sample.rep, mean = is.mean.est,
                                     sd = is.sd.est))

  likelihood.h0 <- stats::dt(t.rep, df.rep)
  likelihood.hr_is <- modelevidence.rep

  rbf <- likelihood.hr_is / likelihood.h0
  if (!store.samples) {
    posterior.sample.orig = NULL
    posterior.sample.rep = NULL
  }

  # Prepare return object ------------------------------------------------------
  ret.object <- list(
    # Numeric value of the Bayes Factor
    bayesFactor = rbf,

    # Samples from the original study's posterior (empty if store.samples = F)
    posteriorSamplesOriginal = posterior.sample.orig,
    posteriorSamplesReplication = posterior.sample.rep,

    # String identifying the RBF test
    test = "Replication Bayes Factor for t-tests",

    # User's function call
    functionCall = match.call(),

    # Used approximation method (Normal or MCMC)
    approxMethod = method,

    # Details of the original study
    originalStudy = list(
      t = t.orig,
      n = n.orig
    ),

    # Details of the replication study
    replicationStudy = list(
      t = t.rep,
      n = n.rep
    )
  )
  class(ret.object) <- "ReplicationBF"

  return(ret.object)
}
