compute_ratio <- function(dividend,divisor,name="",debug=FALSE) {
  ratio <- list( val=dividend$val / divisor$val, 
                 dval=sqrt( (dividend$dval/divisor$val)^2 + (divisor$dval*dividend$val/divisor$val^2)^2 ), 
                 name=name )
  if(debug) {
    print(sprintf("compute_ratio: %s",as.character(name)))
    print(ratio)
  }
  return(ratio)
}

