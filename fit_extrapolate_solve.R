# this is the "fit,extrapolate,solve" function collection for
# doing fits on bootstrap samples of data and extrapolating
# from these fits with essentially arbitrary numbers of paramters

# linear fit of the form y~a1*x1+a2*x2+...+an*xn+c
# to data "dat" which should be a list of data frames
# each component of the list is one bootstrap sample of the data
# The columns of the data frame should be named:
#   y x1 x2 ... xn weight
# where the last column is passed as the weights argument of the nls or lm
# and the intermediate columns correspond to the "x values"
# The fit function is constructed on the fly based on the number of columns
# in dat.

fes_fit_linear <- function(dat,start,type="nls",debug=FALSE,mc=TRUE) {
  lapply.overload <- lapply
  if(mc){
    require("parallel")
    lapply.overload <- mclapply
  }
  if( type != "nls" ) {
    stop("fit_linear: Only fit type 'nls' currently supported!")
  }

  n <- length(dat)
  fit <- NULL
  failed_fits <- 0
  startvals <- list()
  
  # construct a fit function based on the input data
  # and set some starting values
  model <- "y~a1*x1"
  if(missing(start)) {
    startvals[["a1"]] <- 0.1
    startvals[["c"]] <- 1.0
  } else {
    startvals <- start
  }
  
  other_terms <- c(2)
  # if there are more terms, this will take account of them correctly!
  if( length(dat[[1]][1,])-2 >= 2 )
  {
    other_terms <- length(dat[[1]][1,])-2
    for( i in 2:other_terms ) {
      model <- sprintf("%s+a%d*x%d",model,i,i)
      if(missing(start)) {
        startvals[[sprintf("a%s",i)]] <- 0.1 
      }
    }
  }
  model <- sprintf("%s+c",model)
  if(debug) {
    cat("For this fit, the model", model, "is being used.\n")
  }
  
  temp <- lapply.overload(X=dat,FUN=function(x) {
      algorithms <- c("default","port","plinear")
      fit <- NULL
      for(algorithm in algorithms){
        fit <- try(nls(as.formula(model),data=x,start=startvals,weights=x$weight,trace=FALSE,model=TRUE,algorithm=algorithm),silent=TRUE)
        if(!any(class(fit) == "try-error")) break;
      }
      fit
    })
  # remove any fits that failed
  for( index in 1:n ) {
    if( !any(class(temp[[index]])=="try-error") ) {
      fit[[length(fit)+1]] <- temp[[index]]
    }else{
      failed_fits <- failed_fits+1
      cat("fit",index,"failed\n")
    }
  }

  rval <- list(fit=fit,n=(n-failed_fits),model=model)
  attr(rval, "class") <- c("fesfit","fesfit_linear", class(rval))
  return( rval )
}

# use predict (or predictNLS) to extrapolate a model fitted using
# a fes_fit function to some new values of predictor variables

fes_extrapolate <- function(fesfit,pred,debug=FALSE,mc=TRUE) {
  lapply.overload <- lapply
  if(mc){
    require("parallel")
    lapply.overload <- mclapply
  }
  if( !any( class(fesfit) == "fesfit" ) ) {
    stop("fes_extrapolate: 'fesfit' argument needs to be of class 'fesfit'\n")
  }
  
  classes <- NULL 
  for(fit in fesfit$fit){
    classes <- c(classes,class(fit))
  }
  print(unique(classes))

  if( class(fesfit$fit[[1]]) != "nls" ) {
    stop("fes_extrapolate: Only fit type 'nls' currently supported!\n")
  }
  
  # for predictNLS, the model fuction needs to be in the current environment
  # with the same name as it was used in fes_fit
  # the reason is that predicNLS extracts the function from the $m$call member
  # which literally stores the code that was used to pass the model
  # unfortunately, this call then only contains something like
  # "as.formula(model)"
  # and "model" needs to be in the environment if this is to be 
  # evaluated successfully
  # I don't know a solution to this problem yet...
  model<-fesfit$model
  
  # we want a prediction for some new values of the predictor variables
  # in pred 
  predy <- list( y=array(dim=c(fesfit$n,length(pred[,1]))), 
                 dy=array(dim=c(fesfit$n,length(pred[,1]))) )

  require("propagate")
  # suppress some very verbose default output
  sink("/dev/null")
  temp <- lapply.overload(X=fesfit$fit,FUN=function(x) { predictNLS(x,newdata=pred,do.sim=FALSE,interval='prediction')$summary })
  # reset the output!
  sink(NULL)
  sink(NULL)
  # extract content from list
  for( index in 1:fesfit$n ){
    predy$y[index,] <- temp[[index]][,2]
    predy$dy[index,] <- temp[[index]][,4]
  }
  
  return(predy)
}

# taking a fit from a fes_fit function, this function can be supplied
# with some value and error and it will try to find the predictor
# variables at which the function takes these values

# the return value is an array with as many rows as there are bootstrap
# samples which resulted in successful fits
# there are as many columns as there are unknowns, multiplied
# by the number of values which we try to solve for

fes_solve <- function(fesfit,unknown,known,y,dy,interval=c(-10,10),debug=FALSE,mc=TRUE) {
  lapply.overload <- lapply
  if(mc){
    require("parallel")
    lapply.overload <- mclapply
  }
  if( !any( class(fesfit) == "fesfit" ) ) {
    stop("fes_solve: 'fesfit' argument needs to be of class 'fesfit'\n")
  }
  if( missing(unknown) ) {
    stop("fes_solve: 'unknown' argument must be supplied with the name of the predictor variabe which should be solved for")
  }
  if(!missing(known) ) {
    if( any( as.vector( outer(colnames(known),unknown,'==') ) ) ) {
      stop("fes_solve: one or more column names in 'known' match the names of the unknowns, this is surely wrong!\n")
    }
  }
  
  # this needs to be defined in the environment for the extraction below to work
  model <- fesfit$model
  
  find_solutions <- function(fit){
    # extract the rhs of the fit function
    rhs <- as.list(eval(fit$call$formula))[[3]]
    # construct a data frame to act as the environment to evaluate the rhs
    environ <- t(as.data.frame(fit$m$getPars()))
    if(length(known) > 0)
      environ <- cbind(environ,known)
    
    # this is the function which will be minimized with respect to x
    # at the points where it take the value "rootval"
    rootfun <- function(x,rootval) {
      if(length(x) != length(unknown)) {
        stop(sprintf("fes_solve: rootfun supplied with %d unknowns but passed %d values to optimize!\n",length(unknown),length(x)))
      }
      # x is passed as a vector but we need to supply it with names for the
      # unknown(s) that we are looking for
      x.df <- as.data.frame(t(x))
      colnames(x.df) <- unknown
      environ <- cbind(environ,x.df) 
      # evaluate the model with all parameters, knowns anqd numerical values for unknowns
      # the absolute value ensures that we have viable minimum
      rval <- abs(eval(rhs,envir=environ)-rootval)
      return(rval)
    }
    
    findroot <- function(rootval) {
      if(length(unknown)==1) {
        oval <- optimize(f=rootfun, interval=interval, rootval=rootval)
        return( c(oval$minimum,oval$objective) )
      } else {
        oval <- optim(par=rep(0.1,times=length(unknown)), fn=rootfun, rootval=rootval)
        return( c(oval$par,oval$value) )
      }
    }
    
    unlist(lapply(X=y,FUN=findroot))
  }

  temp <- lapply.overload(X=fesfit$fit,FUN=find_solutions)
  # transform list into array
  solutions <- t(array(unlist(temp),dim=c((length(y)*length(unknown)+1),fesfit$n)))

  return(solutions)
}

