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
  if ((length(rbf.object$posteriorSamplesOriginal) == 0) |
      (length(rbf.object$posteriorSamplesReplication) == 0))
    stop("The ReplicationBF object does not contain posterior/prior samples.")

  #stop("Method not yet implemented.")

  # Approximate density function in order to draw points
  density.orig <- approxfun(density(as.numeric(rbf.object$posteriorSamplesOriginal)))
  dens.orig.null <- density.orig(0)
  density.rep <- approxfun(density(as.numeric(rbf.object$posteriorSamplesReplication)))
  dens.rep.null <- density.rep(0)

  # Generate plot --------------------------------------------------------------
  samples.posterior <- data.frame(orig = as.numeric(rbf.object$posteriorSamplesOriginal),
                             rep = as.numeric(rbf.object$posteriorSamplesReplication))
  plt <- ggplot(data = samples.posterior) +
    stat_density(aes(x = orig, linetype = "prior"), geom = "line") +
    stat_density(aes(x = rep, linetype = "posterior"), geom = "line") +
    geom_point(x = 0, y = dens.orig.null, color = "darkgrey", size = 1.8) +
    geom_point(x = 0, y = dens.rep.null, color = "darkgrey", size = 1.8) +
    scale_x_continuous(name = "Effect Size") +
    scale_y_continuous(name = "Density") +
    scale_linetype_manual(values = c("prior" = "dashed",
                                     "posterior" = "solid"),
                          labels = c("prior" = "Original Study",
                                     "posterior" = "Replication Study"),
                          name = "")

  return(plt)
}