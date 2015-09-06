compute_square <- function(x,name=NA,debug=FALSE) {
  rval <- data.frame( val=x$val^2,                                                                                                                                                                             
                 dval=2*x$val*x$dval,
                 name=name )
  if(debug) {
    print(sprintf("compute_square: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_sqrt <- function(x,name=NA,debug=FALSE){
  if(x$val < 0){
    x$val <- abs(x$val)
    cat(sprintf("compute_sqrt: Warning, negative value replaced by absolute value for %s!\n",name))
  }
  rval <- data.frame( val=sqrt(x$val),
                dval=0.5*x$dval/sqrt(x$val),
                name=name )
  if(debug) {
    print(sprintf("compute_sqrt: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}
