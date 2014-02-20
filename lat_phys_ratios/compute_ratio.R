compute_ratio <- function(dividend,divisor,name,debug=FALSE) {
  ratio <- list( value=dividend$value / divisor$value, 
                 error=sqrt( (dividend$error/divisor$value)^2 + (divisor$error*dividend$value/divisor$value^2)^2 ), 
                 name=name )
  if(debug) {
    print(sprintf("compute_ratio: %s",as.character(name)))
    print(ratio)
  }
  return(ratio)
}

