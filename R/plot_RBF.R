#' Plot the results of a Replication Bayes Factor.
#'
#' @param rbf.object Result from \code{\link{RBF_ttest}} or
#'     \code{\link{RBF_Ftest}}.
#'
#' @return None.
#' @export

plot_RBF <- function(rbf.object) {

  # Check arguments ------------------------------------------------------------
  if (class(rbf.object) != "ReplicationBF")
    stop("Please provide a ReplicationBF to plot.")
  if (length(rbf.object$posteriorSample) == 0)
    stop("The ReplicationBF object does not contain posterior/prior samples.")

  stop("Method not yet implemented.")

  prior <- sample(rbf.object$posteriorSample, 1e4, replace = FALSE)
  prior.density <- approxfun(density(prior))


  sqrt.n.rep <- sqrt(1 / (1 / rbf.object$replicationStudy$n[1] +
                            1 / rbf.object$replicationStudy$n[2]) )
  likelihood <- dt(rbf.object$replicationStudy$t,
                   sum(rbf.object$replicationStudy$n) -
                     length(rbf.object$replicationStudy$n),
                   prior*sqrt.n.rep)

  posterior <- (prior.density(prior) * likelihood)

  # Generate plot --------------------------------------------------------------
  plot(density(prior), lwd = 1, lty = 2, col = 1,
       ylab = "Density", xlab = " ", main = rbf.object$test, xlim = c(-1, 3))
  par(new = T)
  plot(density(posterior), lwd = 2, lty = 1, col = 1,
       ylab = "Density", xlab = " ", main = " ", xlim = c(-1, 3))



}