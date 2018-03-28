#' Plot \code{ReplicationBF} object
#'
#' This method can be used to visually represent a Replication Bayes Factor.
#'
#' \code{plot_RBF} generally is a Savage-Dickey ratio representation:
#' Prior and Posterior distributions are plotted and the ratio of the
#' distributions at \eqn{\delta = 0} is calculated (see Wagenmakers et al., 2010).
#'
#' @param rbf.object Object from an \code{RBF_*} call. Has to include Posterior
#'     samples (\code{store.samples = TRUE}).
#' @param use.ggplot If \code{TRUE} uses \code{ggplot2}. \code{FALSE} not yet
#'     supported.
#'
#' @references
#' \itemize{
#'     \item \insertRef{Wagenmakers2010}{ReplicationBF}
#' }
#'
#' @return A \code{ggplot} object containing the plot.
#' @export

plot_RBF <- function(rbf.object, use.ggplot = TRUE) {

  # Check arguments ------------------------------------------------------------
  if (class(rbf.object) != "ReplicationBF")
    stop("Please provide a ReplicationBF to plot.")
  if ((length(rbf.object$posteriorSamplesOriginal) == 0) |
      (length(rbf.object$posteriorSamplesReplication) == 0))
    stop("The ReplicationBF object does not contain posterior/prior samples.")
  if (use.ggplot) {
    # Check for ggplot2 package
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop(paste("To plot the Replication Bayes factor using ggplot2, ",
                 "please install the ggplot2 package.", sep=""))
    }
  } else { stop("Base plots are not yet implemented. Sorry!") }

  # Approximate density function in order to draw points for Savage-Dickey Ratio
  density.orig <- density(as.numeric(rbf.object$posteriorSamplesOriginal))
  density.f.orig <- approxfun(density.orig)
  dens.orig.null <- density.f.orig(0)

  density.rep <- density(as.numeric(rbf.object$posteriorSamplesReplication))
  density.f.rep <- approxfun(density.rep)
  dens.rep.null <- density.f.rep(0)

  # Annotations
  rbf.annotation <- paste0("B[r0] == ", round(rbf.object$bayesFactor, 3))

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
                          name = "") +
    annotate("text",
             x = min(c(samples.posterior$orig, samples.posterior$rep))+1,
             y = max(c(density.orig$y, density.rep$y)),
             hjust = 1, vjust = 1,
             label = rbf.annotation,
             parse = T)

  return(plt)
}
