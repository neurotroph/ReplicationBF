
setMethod('show', 'ReplicationBF', function(object) {
  cat(object$test, "\n")
  cat("========================\n")
  cat("BF_r0 = ", round(object$bayesFactor, 5), "\n")
  cat("BF_0r = ", round(1/object$bayesFactor, 5), "\n\n")
  if (!is.null(object$posteriorSamplesOriginal) &&
      !is.null(object$posteriorSamplesReplication))
    cat("Object contains posterior samples.\n")
  cat("Approximated using ", object$approxMethod, "\n")
  cat("Call:\n")
  cat("  ")
  print(object$functionCall)
})

setMethod('summary', 'ReplicationBF', function(object) {
  show(object)
})

setMethod('print', 'ReplicationBF', function(x) {
  show(x)
})

setMethod('plot', 'ReplicationBF', function(x, ...) {
  plot.ReplicationBF(x, ...)
})

is.ReplicationBF <- function(x)
  inherits(x, "ReplicationBF")

print.ReplicationBF <- function(x) {
  show(x)
}