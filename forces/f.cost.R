c.mu <- function(mu) {
  8.614/(mu+0.001443)
}

rec.x <- function(x,f,n) { ( f^n*(x-1) + f*(x-1) / (f*(f-1)) ) - 1 }

rec.x2 <- function(x,n) {1+(n-1)*x }

f.cost <- function(x,n1,n,f) {
  # simulated cost for gauge update
  cost <- n1*2^n*c.mu(0.16)
   
  for( i in 1:n ) {
    cost <- cost + n1*2^i*c.mu(rec.x2(n=i,x=x)*0.0008238)
  }
  data.frame(cost=cost,n=n1,x=x)
}
  

