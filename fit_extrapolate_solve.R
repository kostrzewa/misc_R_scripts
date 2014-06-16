# carry out a fit of the form "y=a*x+b" on (bootstrapped) data
# 'y' specified at points 'x'
# parameters:
# 'y': 2 dimensional array with the row index specifying
#      the bootstrap sample and the second index corresponding
#      to the 'x' value that he data point belongs to
# 'x': vector specifying the x-values that the data points
#      'y' correspond to. Length equal to number of columns in 'y'
# 'weights': Weights to be given the different 'x' values in the fit
# 'type': 'nls' or 'lm', specifiying type of fit to be carried out
#         'nls' type required for proper error propagation using
#         'predictNLS'
#
# return value: list of lists containing objects of class 'type

fit_linear_1d <- function(y,x,weights,type="lm",debug=F) {
  if(length(y[1,])!=length(x) || length(y[1,])!=length(weights)){
    stop(sprintf("fit_linear_1d: Mismatch between the lengths of y(%d),x(%d) and weights(%d)!",length(x),length(y[1,]),length(weights)))
  }
  
  n <- length(y[,1])
  
  # return value
  fit <- NULL
  # counter for fits that didn't converge for some reason
  failed_fits <- 0
  
  for( index in 1:n ) {
    dat <- data.frame(y=y[index,],x=x)
    if(debug) cat("Fitting on bootstrap sample", index, "\n")
    if(type == "nls") {
      temp <- try(nls(y~a*x+b,data=dat,start=list(a=0.1,b=1),weights=weights,trace=FALSE))
      if( any(class(temp) == "try-error" )) {
        cat("Fit for bootstrap sample", index, "failed!\n")
        failed_fits <- failed_fits+1
      } else {
        fit[[(index-failed_fits)]] <- temp
      }
    } else if (type == "lm") {
      fit[[index]] <- lm(y~x,data=dat,weights=weights)
    } else {
      stop("fit_linear_1d: unknown fit type",type)
    }
  }
  
  rval <- list(fit=fit,n=(n-failed_fits))
  attr(rval, "class") <- c("fesfit_1d","fesfit_linear", class(rval))
  return( rval )
} 

# linear fit to data 'z' of the form 'z=a*x + b*y + c'
# where 'z' is a data frame with as many rows as there are
# bootstrap samples and columns 'z','x','y' for each

fit_linear_2d <- function(z,weights,type="nls",debug=F) {
  if( type != "nls" ) {
    stop("fit_linear_2d: Only fit type 'nls' currently supported!")
  }

  n <- length(z)
  fit <- NULL
  failed_fits <- 0
  
  for( index in 1:n ) {
    temp <- try(nls(z~a*x+b*y+c,data=z[[index]],start=list(a=0.1,b=0.1,c=1),weights=weights,trace=FALSE))
    if( any(class(temp) == "try-error" )) {
      cat("Fit for bootstrap sample", index, "failed!\n")
      failed_fits <- failed_fits+1
    } else {
      if(debug) print(temp)
      fit[[(index-failed_fits)]] <- temp
    }
  }
  
  rval <- list(fit=fit,n=(n-failed_fits))
  attr(rval, "class") <- c("fesfit_2d","fesfit_linear", class(rval))
  return( rval )
}

fit_linear_1d.new <- function(dat,weights,type="nls",debug=F) {
  if( type != "nls" ) {
    stop("fit_linear_2d: Only fit type 'nls' currently supported!")
  }

  n <- length(dat)
  fit <- NULL
  failed_fits <- 0
  
  for( index in 1:n ) {
    temp <- try(nls(y~a*x+b,data=dat[[index]],start=list(a=0.1,b=0.1),weights=weights,trace=FALSE))
    if( any(class(temp) == "try-error" )) {
      cat("Fit for bootstrap sample", index, "failed!\n")
      failed_fits <- failed_fits+1
    } else {
      fit[[(index-failed_fits)]] <- temp
    }
  }
  
  rval <- list(fit=fit,n=(n-failed_fits))
  attr(rval, "class") <- c("fesfit_1d","fesfit_linear", class(rval))
  return( rval )
}

# extrapolation in 1d on (bootstrapped) fits produced by
# the fit routines in this file

extrapolate_1d <- function(fit,predx,dpredx,debug=FALSE) {
  if( !any( class(fit) == "fesfit_1d" ) ) {
    stop("extrapolate_1d: 'fit' argument needs to be of class 'fesfit_1d'")
  }

  # we want a prediction for some new values of x 
  newy <- NULL
  if(missing(dpredx) || class(fit$fit[[1]]) != "nls" ) {
    newy <- list( y=array(dim=c(fit$n,length(predx))), dy=NA )
    for(index in 1:fit$n) { 
      newy$y[index,] <- predict(fit$fit[[index]],newdata=data.frame(x=predx))
    }
  } else {
    if( class(fit$fit[[1]]) != "nls" ) {
      stop("extrapolate_1d: For doing predict with an error in the predictor variable, the fit type must be \"nls\"!")
    }
    if(debug) cat("Doing predict with errors in predictor variables\n")
    newy <- list( y=array(dim=c(fit$n,length(predx))), 
                  dy=array(dim=c(fit$n,length(predx))) )
    for(index in 1:fit$n) {
      prednls <- predictNLS(fit$fit[[index]],newdata=data.frame(x=predx,dx=dpredx),do.sim=FALSE,interval='prediction')$summary 
      newy$y[index,] <- prednls[,2]
      newy$dy[index,] <- prednls[,4]
    }
  }
  return(newy)
}

extrapolate_2d <- function(fit,predx,predy,dpredx,dpredy,debug=FALSE) {
  if( !any( class(fit) == "fesfit_2d" ) ) {
    stop("extrapolate_2d: 'fit' argument needs to be of class 'fesfit_2d'")
  }
  
  if( class(fit$fit[[1]]) != "nls" ) {
    stop("extrapolate_2d: Only fit type 'nls' currently supported!")
  }
  
  # we want a prediction for some new values of x 
  newz <- list( z=array(dim=c(fit$n,length(predx))), 
                  dz=array(dim=c(fit$n,length(predx))) )
  
  if( missing(dpredx) && missing(dpredy) ) {
    for(index in 1:fit$n) { 
      newz$z[index,] <- predict(fit$fit[[index]],newdata=data.frame(x=predx,y=predy))
      newz$dz[index,] <- rep(0,length(predx))
    }
  } else {
    if( class(fit$fit[[1]]) != "nls" ) {
      stop("extrapolate_2d: For doing predict with an error in the predictor variables, the fit type must be \"nls\"!")
    }
    if(debug) cat("Doing predict with errors in predictor variables\n")
    for(index in 1:fit$n) {
      prednls <- predictNLS(fit$fit[[index]],newdata=data.frame(x=predx,y=predy,dx=dpredx,dy=dpredy),do.sim=FALSE,interval='prediction')$summary 
      newz$z[index,] <- prednls[,2]
      newz$dz[index,] <- prednls[,4]
    }
  }
  return(newz)
}

solve_linear_1d <- function(fit,predy,dpredy) {
  if( !all( c('fesfit_1d','fesfit_linear') %in% class(fit) ) ) {
    stop("extrapolate_1d: 'fit' argument needs to be of class 'fesfit_1d' AND 'fesfit_linear'")
  }

  newx <- NULL
  tempa <- NULL
  if(missing(dpredy)) {
    newx <- list( x=array(dim=c(fit$n,length(predy))), dx=NA )
  } else {
    newx <- list( x=array(dim=c(fit$n,length(predy))), dx=array(dim=c(length(predy))) )
    tempa <- array(dim=c(length(fit$n)))
  } 
  for(index in 1:fit$n) {
    a <- NULL
    b <- NULL
    if( class(x=fit$fit[[index]]) == "lm" ) {
      a <- fit$fit[[index]]$coefficients[2]
      b <- fit$fit[[index]]$coefficients[1]
    } else if ( class(x=fit$fit[[index]]) == "nls" ) {
      a <- fit$fit[[index]]$m$getPars()[1]
      b <- fit$fit[[index]]$m$getPars()[2]
    } else {
      stop("solve_1d: The fit(s) must be of class \"lm\" or \"nls\"!")
    }
    newx$x[index,] <- sapply(X=predy-b,FUN=solve,a=a)
    
    if(!missing(dpredy)) {
      tempa[index] <- a
    }
  }
  if(!missing(dpredy)) {
    # approximately propagate the error from dpredy, this can then later be combined
    # with the bootstrap error in quadrature 
    newx$dx <- dpredy/mean(tempa)
  }
  return(newx)
}

