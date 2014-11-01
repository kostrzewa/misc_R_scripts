compute_ratio <- function(dividend,divisor,name=NA,debug=FALSE) {
  rval <- list( val=dividend$val / divisor$val, 
                 dval=sqrt( (dividend$dval/divisor$val)^2 + (divisor$dval*dividend$val/divisor$val^2)^2 ), 
                 name=name )
  if(debug) {
    print(sprintf("compute_ratio: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_product <- function(a,b,name=NA,debug=FALSE) {
  rval <- list( val=a$val * b$val, 
                 dval=sqrt( (a$dval*b$val)^2 + (b$dval*a$val)^2 ), 
                 name=name )
  if(debug) {
    print(sprintf("compute_product: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}


