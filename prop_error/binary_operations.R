compute_ratio <- function(dividend,divisor,name=NA,debug=FALSE) {
  rval <- data.frame( val=dividend$val / divisor$val, 
                      dval=sqrt( (dividend$dval/divisor$val)^2 + (divisor$dval*dividend$val/divisor$val^2)^2 ), 
                      name=name )
  if(debug) {
    print(sprintf("compute_ratio: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_product <- function(a,b,name=NA,debug=FALSE) {
  rval <- data.frame( val=a$val * b$val, 
                      dval=sqrt( (a$dval*b$val)^2 + (b$dval*a$val)^2 ), 
                      name=name )
  if(debug) {
    print(sprintf("compute_product: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_sum <- function(a,b,name=NA,debug=FALSE) {
  rval <- data.frame( val=a$val + b$val,
                      dval=sqrt( a$dval^2 + b$dval^2 ),
                      name=name )
  if(debug) {
    cat(sprintf("compute_sum: %s\n",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_difference <- function(pos,neg,name=NA,debug=FALSE) {
  neg$val <- -neg$val
  rval <- compute_sum(a=pos,b=neg,name=name,debug=FALSE)
  if(debug) {
    cat(sprintf("compute_difference: %s\n",as.character(name)))
    print(rval)
  }
  return(rval)
}

