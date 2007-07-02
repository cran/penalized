plotpath <- function(object, labelsize = 0.6,...) {
  if (length(object) > 0) {
    betas <- sapply(object, coefficients, "p")
  }
  remove <- apply(betas, 1, function(bet) all(bet == 0) )
  lambda1 <- sapply(object, function(object) object@lambda1)
  labwidth <- ifelse(labelsize > 0, max(strwidth(rownames(betas[!remove,]),"inches",labelsize)), 0)
  margins <- par("mai")
  par("mai" = c(margins[1:3], max(margins[3], labwidth*1.4)))
  matplot(lambda1, t(betas[!remove,]), type ="l", ylab = "coefficient", xlab = "L1-penalty",col=rainbow(sum(!remove)), xlim = rev(range(lambda1)), log="x" )
  if (labelsize > 0) {
    take <- which(!remove)
    for (i in 1:sum(!remove)) {
      j <- take[i]
      axis(4, at = betas[j,ncol(betas)], labels = rownames(betas)[j],
      las=1,cex.axis=labelsize, col.axis=rainbow(sum(!remove))[i], lty = (i-1) %% 5 + 1, col = rainbow(sum(!remove))[i])
    }
  }
  par("mai"=margins)
}
