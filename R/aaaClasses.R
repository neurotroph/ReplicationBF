#' General S4 class for a single Replication Bayes Factor
#'
#' @slot bayesFactor The numerical value of the Replication Bayes Factor.
#' @slot approxMethod Method used to approximate the marginal likelihoods.
#' @slot posteriorSamplesOriginal A \code{data.frame} containing the samples
#'     from the original study's posterior distribution.
#' @slot posteriorSamplesReplication A \code{data.frame} containing samples
#'     from the replication study's posterior distribution.
#' @slot test A string that contains which test was used.
#' @slot functionCall String containing the function call.
#' @slot originalStudy Data from the original study that went into the analysis.
#' @slot replicationStudy Data from the replication study that was used for the
#'     calculation of the Bayes factor.
#'
#' @export

setClass("ReplicationBF", representation(
  test = "character",
  functionCall = "character",
  approxMethod = "character",
  bayesFactor = "numeric",
  posteriorSamplesOriginal = "numeric",
  posteriorSamplesReplication = "numeric",
  originalStudy = "list",
  replicationStudy = "list"
))