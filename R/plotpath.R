plotpath <- function(object, labelsize = 0.6,...) {
  if (length(object) > 0) {
    betas <- sapply(object, coefficients, "p")
  }
  # Do not plot the regression coefficients of covariates that are always zero
  remove <- apply(betas, 1, function(bet) all(bet == 0) )

  # Take lambda1, unless all lambda1 are equal. Then take lambda2
  lambda <- sapply(object, function(object) object@lambda1)
  label <- "lambda1"
  if (all(lambda == lambda[1])) {
    lambda <- sapply(object, function(object) object@lambda2)
    label <- "lambda2"
  }
  
  # Adjust the margins to make sure the labels fit
  labwidth <- ifelse(labelsize > 0, max(strwidth(rownames(betas[!remove,]),"inches",labelsize)), 0)
  margins <- par("mai")
  par("mai" = c(margins[1:3], max(margins[3], labwidth*1.4)))
  
  # Plot
  matplot(lambda, t(betas[!remove,,drop=FALSE]), type ="l", ylab = "coefficient", xlab = label, col=rainbow(sum(!remove)), xlim = rev(range(lambda)), ...)
  if (labelsize > 0) {
    take <- which(!remove)
    for (i in 1:sum(!remove)) {
      j <- take[i]
      axis(4, at = betas[j,ncol(betas)], labels = rownames(betas)[j],
      las=1,cex.axis=labelsize, col.axis=rainbow(sum(!remove))[i], lty = (i-1) %% 5 + 1, col = rainbow(sum(!remove))[i])
    }
  }
  
  # Reset the margins
  par("mai"=margins)
  
  return(invisible(NULL))
}
