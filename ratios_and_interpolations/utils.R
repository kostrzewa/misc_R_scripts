plot.confband <- function(x,y,dy,col,line=T) {
  xvals <- c(rev(x),x)
  yvals <- c(rev(y-dy),y+dy)
  polygon(x=xvals,y=yvals,col=col,border=NA)
  if(line) lines(x=x,y=y)
}

compute.quant <- function(name,texlabel,dat,res) {
  n <- length(dat$dat)
  
  rval.tsboot <- data.frame( name=rep(name,n), texlabel=rep(texlabel,n), val=dat$dat,
                m11=rep(dat$m1,n), m12=rep(dat$m2,n), m21=rep(NA,n), m22=rep(NA,n) )
  rval <- data.frame( name=name, texlabel=texlabel, val=mean(rval.tsboot$val), dval=sd(rval.tsboot$val),
                      m11=dat$m1, m12=dat$m2, m21=NA, m22=NA )
  
  return( append.result(rval,rval.tsboot,res) )
}

#compute.product <- function(name,texlabel,a,b,res) {
#  n <- length(a$dat)
#  
#  rval.tsboot <- 
#}

compute.expression <- function(expr,envir,m1,m2){
  return( list( dat=eval(expr=expr,envir=envir), m1=m1, m2=m2 ) )
}

compute.ratio <- function(name,texlabel,dividend,divisor,res) {
  n <- length(dividend$dat)
  if( n != length(divisor$dat) ) {
    stop("In compute.ratio, dividend and divisor do not have the same length!")
  }
  rval.tsboot <- data.frame( name=rep(name,n), texlabel=rep(texlabel,n), val=dividend$dat/divisor$dat,
                m11=rep(dividend$m1,n), m12=rep(dividend$m2,n), m21=rep(divisor$m1,n), m22=rep(divisor$m2,n) )
  rval <- data.frame( name=name, texlabel=texlabel, val=mean(rval.tsboot$val), dval=sd(rval.tsboot$val),
                      m11=dividend$m1, m12=dividend$m2, m21=divisor$m1, m22=divisor$m2 )
                      
  return( append.result(val=rval,val.tsboot=rval.tsboot,res) )
}

append.result <- function(val,val.tsboot,res) {
  if( length(res$val) == 0 ) {
    return( list(val=val,val.tsboot=val.tsboot) )
  } else {
    return( list(val=rbind(res$val,val), val.tsboot=rbind(res$val.tsboot,val.tsboot) ) )
  }
}

get.obs.tsboot <- function(obj,type) {
  rval <- list( m1=obj$mu1, m2=obj$mu2, boot=c() )
  if( type == "fps" ) {
    rval$boot <- obj$fps.tsboot
  } else if ( type == "mps" ) {
    rval$boot <- obj$opt.tsboot[1,]
  } else {
    stop("In 'get.obs.tsboot': type '", type,"' is not known. Exiting!")
  }
  return(rval)
}
