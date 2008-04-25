###################################                                                           
# These functions are not user level functions
# They require a very specific input format and they 
# rely on the functions calling them for input checking
# These functions are not exported in the NAMESPACE file
###################################
 

###################################
# The core lasso and elastic net algorithm
###################################
.lasso <- function(beta, lambda, lambda2=0, positive, X, fit, trace = FALSE, 
  epsilon, maxiter) {
         
  # It is a general function for fitting L1-penalized models
  # possibly with an additional L2-penalty
  # Input:
  #  beta: a vector of length m (say) : starting values
  #  lambda: a vector of length m
  #  fit must be function(beta)
  #   Should return a list with at least:
  #     W:          The weights matrix, or its diagonal
  #     loglik:     The unpenalized loglikelihood, numeric)
  #     residuals:  The residuals
  
  m <- length(beta)
  n <- nrow(X)
  # are we fitting an elastic net?
  enet <- any(lambda2 != 0)   
  # find regression coefficients free of L1-penalty or positivity restraint
  free <- lambda == 0 & !positive      
  
  # initialize
  LL <- -Inf
  penalty <- penalty1 <- penalty2 <- Inf

  active <- !logical(m)
  nvar <- m

  tryNR <- FALSE
  NRfailed <- FALSE   
  whereNR <- NULL

  finished <- FALSE
  newfit <- TRUE
  retain <- 0.05
  cumsteps <- 0
  iter <- 0

  # iterate
  if (trace) cat("# nonzero coefficients:", m)
  while (!finished) {
    nzb <- (beta != 0)
    # calculate the local likelihood fit
    if (newfit) {
      activeX <- X[,nzb, drop=FALSE]
      linpred <- drop(activeX %*% beta[nzb])
      localfit <- fit(linpred)
      # Check for divergence
      if (is.na(localfit$loglik)) {
        if (trace) {
          cat(rep("\b", trunc(log10(nvar))+1), sep ="")
          warning("Model does not converge: please increase lambda.", call.=FALSE)
        }
        converged <- FALSE
        break
      }
      grad <- drop(crossprod(X, localfit$residuals))
      if (enet) {
        grad[active] <- grad[active] - lambda2[active] * beta[active]
      }
      oldLL <- LL
      oldpenalty <- penalty
      LL <- localfit$loglik
      penalty1 <- sum(lambda[active] * abs(beta[active]))
      if (enet) {
        penalty2 <- sum(lambda2[active] * beta[active] * beta[active])
      } else {
        penalty2 <- 0
      }
      penalty <- penalty1 + penalty2
      finishedLL <- (2 * abs(LL - oldLL) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      finishedpen <- (2 * abs(penalty - oldpenalty) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      cumsteps <- 0
    }

    # Calculate the penalized gradient from the likelihood gradient
    direction <- numeric(m)
    direction[nzb] <- grad[nzb] - lambda[nzb] * sign(beta[nzb])
    newb <- (!nzb) & ifelse(positive, grad > lambda, abs(grad) > lambda) 
    direction[newb] <- grad[newb] - lambda[newb] * sign(grad[newb])
    oldactive <- active   
    active <- nzb | newb
    activebeta <- beta[active]
    activedir <- direction[active]

    # check if retaining the old fit of the model does more harm than good
    oldnvar <- nvar
    nvar <- sum(active)
    if ((oldLL - oldpenalty > LL - penalty) || (nvar > 1.1* oldnvar)) {
      retain <- 0.5 * retain
    }

    # check convergence
    finishednvar <- !any(xor(active, oldactive))
    finished <- (finishedLL && finishedpen && finishednvar) || (all(activedir == 0)) || (iter == maxiter)

    if (!finished) {
      iter <- iter+1
      
      # Try Newton-Raphson, using the ridge routine if an L2 penalty is present
      if (tryNR) {
        activeX <- X[,active,drop=FALSE]
        if (enet && nvar > (n+1+sum(free)) ) {
          if (is.null(whereNR) || any(xor(whereNR, active))) {
            whereNR <- active
            P <- .makeP(activeX, lambda2[active], lambda[active], sign(activebeta))
            gams <- .solve(crossprod(t(P)), P %*% activebeta)
            PX <- P %*% t(activeX)
            Pl <- P * matrix(sqrt(lambda2[active]), nrow(P), ncol(P), byrow = TRUE)
            PlP <- crossprod(t(Pl))
          }
          if (is.list(localfit$W)) {
            PXdW <- PX * matrix(sqrt(localfit$W$diagW), nrow(PX), ncol(PX), byrow=TRUE)
            hessian <- - crossprod(t(PXdW)) + crossprod(t(PX %*% localfit$W$P)) - PlP
          } else if (length(localfit$W) > 1) {
            PXW <- PX * matrix(sqrt(localfit$W), nrow(PX), ncol(PX), byrow=TRUE)
            hessian <- -crossprod(t(PXW)) - PlP
          } else {
            hessian <- -crossprod(t(PX)) - PlP
          }
          Pgrad <- P %*% direction[active]
          shg <- drop(.solve(hessian, Pgrad))
          gams <- gams - shg
          NRbeta <- drop(crossprod(P, gams))
        } else { 
          if (is.list(localfit$W)) {
            XdW <- activeX * matrix(sqrt(localfit$W$diagW), nrow(activeX), ncol(activeX))
            hessian <- -crossprod(XdW) + crossprod(crossprod(localfit$W$P, activeX))
          } else if (length(localfit$W) > 1) {
            XW <- activeX * matrix(sqrt(localfit$W), nrow(activeX), ncol(activeX))
            hessian <- -crossprod(XW)
          } else {
            hessian <- -crossprod(activeX)
          } 
          if (enet) diag(hessian) <- diag(hessian) - lambda2[active]
          NRbeta <- activebeta - drop(.solve(hessian, direction[active]))
        }   
        NRfailed <- !all(sign(NRbeta) == sign(activebeta))
        if (!NRfailed) { 
          beta[active] <- NRbeta
          newfit <- TRUE
        } 
      } 

      if (!tryNR || NRfailed) {
        # find the second derivative of the likelihood in the projected direction
        if (newfit) {
          Xdir <- drop(X[,active, drop=F] %*% activedir)
          if (is.list(localfit$W)) {
            curve <- (sum(Xdir * Xdir * localfit$W$diagW) - drop(crossprod(crossprod(localfit$W$P, Xdir)))) / sum(activedir * activedir)
          } else if (length(localfit$W) > 1) {
            curve <- sum(Xdir * Xdir * localfit$W) / sum(activedir * activedir)
          } else {
            curve <- sum(Xdir * Xdir) / sum(activedir * activedir)
          }
          if (enet) {
            curve <- curve + sum(lambda2[active] * activedir * activedir) / sum(activedir * activedir)
          }
          topt <- 1 / curve
        }

        # how far can we go in the calculated direction before finding a new zero?
        tedge <- numeric(m)
        tedge[active] <- -activebeta / activedir
        tedge[tedge <= 0] <- 2 * topt
        tedge[free] <- 2* topt
        mintedge <- min(tedge)

        # recalculate beta
        if (mintedge + cumsteps < topt) {
          beta[active] <- activebeta + mintedge * activedir   
          beta[tedge == mintedge] <- 0  # avoids round-off errors
          cumsteps <- cumsteps + mintedge
          newfit <- (cumsteps > retain * topt) || (nvar == 1)  
          NRfailed <- FALSE
          tryNR <- FALSE
        } else {
          beta[active] <- activebeta + (topt - cumsteps) * activedir
          tryNR <- (cumsteps == 0) && !NRfailed && finishednvar && (enet || nvar < n)
          newfit <- TRUE
        }
      }
    } else {
      converged <- (iter < maxiter)
    }
    if (trace) {
      cat(rep("\b", max(1,trunc(log10(oldnvar))+1)), sep ="")
      cat(nvar)
      flush.console()
    }  
  }
  
  return(list(beta = beta, fit = localfit, penalty = c(L1 = penalty1, L2 = penalty2), iterations = iter, converged = converged))
}



###################################
# Adjusted lasso algorithm
# Tries to prevent large models by first fitting at higher values of lambda         
# Often faster for "pure" lasso
# Not recommended for elastic net
###################################
.steplasso <- function(beta, lambda, lambda2=0, positive, X, fit, trace = FALSE, epsilon, maxiter = Inf) {

  n <- nrow(X)
  finished <- FALSE
  while (!finished) {
    nzb <- (beta != 0)
    lp <- X[,nzb,drop=FALSE] %*% beta[nzb]
    gradient <- drop(crossprod(X[,!nzb,drop=FALSE], fit(lp)$residuals))
    rel <- gradient / lambda[!nzb]
    rel <- rel[rel>0 | !positive[!nzb]]
    if (length(rel) > n) {
      nextlambda <- sort(abs(rel), decreasing = TRUE)[n]
    } else {
      nextlambda <- 1
    }
    if (nextlambda <= 1) {
      finished <- TRUE
    }
    if(!finished) 
      out <- .lasso(beta, nextlambda * lambda, lambda2, positive, X, fit, trace, 1e-4, maxiter)
    else
      out <- .lasso(beta, lambda, lambda2, positive, X, fit, trace, epsilon, maxiter)
    beta <- out$beta
    if (trace && ! finished) cat(rep("\b", 24 + max(1,trunc(log10(sum(beta!=0))+1))), sep = "")
  }
  out
}
  

###################################
# The core ridge algorithm
###################################
.ridge <- function(beta, eta, Lambda, X, fit, trace = FALSE, epsilon = 1e-8, maxiter = 25) {

  if (missing(eta)) eta <- drop(X %*% beta)

  iter <- 0
  oldLL <- -Inf
  finished <- FALSE

  while (!finished)
  {
    localfit <- fit(eta)
    if (is.na(localfit$loglik) || iter == maxiter) {
      if (trace) {
        cat(rep("\b", trunc(log10(iter))+1), sep ="")
        warning("Model does not converge: please increase lambda.", call.=FALSE)
      }
      converged <- FALSE
      break
    }
    if (is.matrix(Lambda)) {
      penalty <- as.numeric(0.5 * sum(beta * (Lambda %*% beta)))
    } else {
      penalty <- as.numeric(0.5 * sum(Lambda * beta * beta))
    }
    LL <- localfit$loglik - penalty
    
    # Check convergence
    finished <- ( 2 * abs(LL - oldLL) / (2 * abs(LL) + 0.1) < epsilon ) | (iter >= maxiter)
    oldLL <- LL
    
    if (!finished) {
      iter <- iter + 1
      if (trace) {
        cat(iter)
        flush.console()
      }
      if (is.matrix(Lambda)) {
        grad <- crossprod(X, localfit$residuals) - Lambda %*% beta
      } else {
        grad <- crossprod(X, localfit$residuals) - Lambda * beta
      }
      if (is.list(localfit$W)) {
        XdW <- X * matrix(sqrt(localfit$W$diagW), nrow(X), ncol(X))
        Hess <-  -crossprod(XdW) + crossprod(crossprod(localfit$W$P, X))
      } else if (length(localfit$W) > 1) {
        XW <- X * matrix(sqrt(localfit$W), nrow(X), ncol(X))
        Hess <- -crossprod(XW) 
      } else {  
        Hess <- -crossprod(X)
      }
      if (is.matrix(Lambda)) {
        Hess <- Hess - Lambda
      } else {
        diag(Hess) <- diag(Hess) - Lambda
      }
      shg <- drop(.solve(Hess, grad))
      beta <- beta - shg       
      eta <- drop(X %*% beta)
      if (is.matrix(Lambda)) {
        penalty <- as.numeric(0.5 * sum(beta * (Lambda %*% beta)))
      } else {
        penalty <- as.numeric(0.5 * sum(Lambda * beta * beta))
      }
      if (trace) cat(rep("\b", trunc(log10(iter))+1), sep ="")
    } else {
      converged <- (iter < maxiter)
    }
  }

  return(list(beta = beta, penalty = c(L1 = 0, L2 = penalty), fit = localfit, iterations = iter, converged = converged))
}


###################################
# Workhorse function for cross-validated likelihood
###################################
.cvl <- function(X, lambda1, lambda2, positive, beta, fit, groups, trace = FALSE, 
  betas = NULL, quit.if.failed = TRUE, ...)  {
                                    
  n <- nrow(X)
  m <- ncol(X)

  # find the right fitting procedure
  useP <- FALSE
  if (all(lambda1 == 0) && !any(positive)) {
    if (m <= n) {
      cvfit <- function(leftout, beta) {
        subfit <- function(lp) fit$fit(lp, leftout)
        .ridge(beta = beta, Lambda = lambda2, X = X[!leftout,,drop = FALSE], fit = subfit, ...)
      }
    } else {
      useP <- TRUE
      P <- .makeP(X, lambda2)
      PX <- P %*% t(X)
      Pl <- P * matrix(sqrt(lambda2), nrow(P), ncol(P), byrow = TRUE)
      cvfit <- function(leftout, beta) {
        subfit <- function(lp) fit$fit(lp, leftout)
        leftoutP <- c(rep(FALSE, nrow(P) - length(leftout)), leftout)
        gams <- .solve(crossprod(t(P[!leftoutP,,drop=FALSE])), P[!leftoutP,,drop=FALSE] %*% beta)
        PlP <- crossprod(t(Pl[!leftoutP,,drop=FALSE]))
        out <- .ridge(beta = gams, Lambda = PlP, X = t(PX[!leftoutP,!leftout,drop = FALSE]), 
          fit = subfit, ...)
        out$beta <- drop(crossprod(P[!leftoutP,,drop=FALSE], out$beta))
        names(out$beta) <- names(beta)
        out
      }
    } 
  } else if (all(lambda2 == 0)) {
    cvfit <- function(leftout, beta) {
      subfit <- function(lp) fit$fit(lp, leftout)
      .steplasso(beta = beta, lambda = lambda1, positive = positive, 
        X = X[!leftout,,drop = FALSE], fit = subfit, ...)
    }
  } else {
    cvfit <- function(leftout, beta) {
      subfit <- function(lp) fit$fit(lp, leftout)
      .lasso(beta = beta, lambda = lambda1, lambda2 = lambda2, positive = positive, 
        X = X[!leftout,,drop = FALSE], fit = subfit, ...)
    }
  }

  # "groups" input lists fold allocation for each subject %in% 1:fold
  fold <- max(groups)
  
  # fit the full model and make an m x fold matrix of beta if necessary
  fullfit <- cvfit(logical(n), beta)
  if (is.null(betas)) {
    if (fullfit$converged) 
      betas <- matrix(fullfit$beta, m, fold)
    else
      betas <- matrix(beta, m, fold)
  } 
  
  # True cross-validation starts here
  if (fold > 1) {
    failed <- FALSE
    predictions <- vector("list", fold)
    cvls <- sapply(1:fold, function(i) {
      if (!failed) {
        if (trace) {
          cat(i)
          flush.console()
        }
        leaveout <- (groups == i)
        foldfit <- cvfit(leaveout, betas[,i])
        lin.pred <- numeric(n)
        lin.pred[leaveout] <- X[leaveout, foldfit$beta != 0, drop=FALSE] %*% 
          foldfit$beta[foldfit$beta != 0]
        lin.pred[!leaveout] <- foldfit$fit$lp0
        predictions[[i]] <<- fit$prediction(lin.pred[leaveout], foldfit$fit$nuisance)
        betas[,i] <<- foldfit$beta
        if (trace) cat(rep("\b", trunc(log10(i))+1), sep ="")
        out <- fit$cvl(lin.pred, leaveout)
        if (quit.if.failed && (is.na(out) || abs(out) == Inf || foldfit$converged == FALSE)) failed <<- TRUE
      } else {
        out <- NA
      }
      
      out
    })
                                      
    if (failed || any(is.na(cvls))) cvls <- -Inf
    
  } else {
    cvls <- NA
    predictions <- NA
    betas <- NA
  }

  list(cvl = sum(cvls), fit = fullfit, betas = betas, predictions = predictions)
}

#######################################
# makes a reduced basis for L2-penalized newton-raphson in the p > n case
#######################################
.makeP <- function(X, lambda2, lambda1 = 0, signbeta) {

  n <- nrow(X)
  p <- ncol(X)
                                                      
  free2 <- (lambda2 == 0)
  free1 <- all(lambda1 == 0)
  m <- sum(free2)
  
  if (free1) {
    P <- matrix(0, n+m, p)
  } else {
    P <- matrix(0, n+m+1, p)
  }

  # First columns: free variables in case of no L2-penalization
  for (i in seq_len(m)) P[i, i] <- 1 
                                                
  # Next n columns: column span of X
  P[m + 1:n, m + seq_len(p-m)] <- X[,!free2,drop=FALSE] * matrix(1/lambda2[!free2], n, p-m, byrow=TRUE)         

  # Additional column due to L1-penalization
  if (!free1) {          
    P[n+m+1,m + seq_len(p-m)] <- (lambda1*signbeta/lambda2)[!free2]
  }
  
  # Numerical stabilization          
  #P <- P / matrix(apply(P, 1, sd), nrow(P), ncol(P), byrow=F)

  return(P)
}

#######################################
# a solve() function that does not complain about near-singularity
# often dangerous, but very useful here
#######################################
.solve <- function(a,b) {
  out <- try(qr.coef(qr(a, LAPACK=TRUE), b))
  if (is(out, "try-error")) stop("Matrix inversion failed. Please increase lambda1 and/or lambda2", call. = FALSE)
  return(out)
}

#######################################
# calculate cross-validation folds
#######################################
.getFolds <- function(fold, n) {
  if(length(fold) == 1) {
    if (fold == n) {
      groups <- 1:n
    } else {
      groups <- sample(n) %% fold + 1
    }
  } else {
    if (length(fold) == n) {
      groups <- fold 
      fold <- max(groups)
    } else {
      stop("incorrect input of \"fold\"", call.=FALSE)
    }
  }
  if (!all(1:fold %in% groups)) stop("incorrect input of \"fold\"", call.=FALSE)
  return(groups)
}