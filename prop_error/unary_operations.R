compute_square <- function(x,name=NA,debug=FALSE) {
  rval <- list( val=x$val^2,                                                                                                                                                                             
                 dval=2*x$val*x$dval,
                 name=name )
  if(debug) {
    print(sprintf("compute_square: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}

