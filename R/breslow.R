# A "breslow" object for storing multiple survival curves on the same time scale
setClass("breslow",
  representation(
    time = "vector",
    curves = "matrix"
  )
)

setMethod("show", "breslow", function(object) {
  cat("A \"breslow\" object with", nrow(object@curves), "survival curve")
  if (nrow(object@curves)>1) cat("s")
  cat(" and", ncol(object@curves), "time points.\n")
})

setMethod("plot", "breslow", function(x, y, ...) {
  plot(0,1,col=0,xlim=range(x@time),ylim=0:1, ylab="",xlab="", ...)
  for (i in 1:nrow(x@curves)) lines(x@time, x@curves[i,], type="s", ...)
  return(invisible(NULL))
})

setMethod("as.matrix", "breslow", function(x, ...) {
  out <- x@curves
  colnames(out) <- x@time
  out
})

setMethod("[", "breslow", function(x, i, j, ... , drop = TRUE) {
  if (missing(i) && missing(j))
    as.matrix(x)[,,drop=drop]
  else if (missing(i))
    as.matrix(x)[,j,drop=drop]
  else if (missing(j))
    as.matrix(x)[i,,drop=drop]
  else
    as.matrix(x)[i,j,drop=drop]
})

setMethod("[[", "breslow", function(x,i,j) {
  x@curves <- x@curves[i,,drop=FALSE]
  x
})

setMethod("time", "breslow", function(x, ...) {
  x@time
})


setMethod("as.list", "breslow", function(x, ...) {
  list(time = x@time, curves = x@curves)
})


