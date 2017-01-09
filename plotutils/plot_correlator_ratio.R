plot_correlator_ratio <- function(dividend,divisor,trange,basename="corratio",width=2.6,height=2.6) {
  if(all(dim(dividend) != dim(divisor))){
    stop("Correlation functions have different dimensions, cannot proceed!\n")
  }
  if(missing(trange)){
    trange <- 0:(dim(dividend)[2]-1)
  }
  ind <- trange+1

  # this should be generalised in the future
  val <- apply(X=dividend/divisor,FUN=mean,MARGIN=2)
  dval <- apply(X=dividend/divisor,FUN=sd,MARGIN=2) 

  tikzfiles <- tikz.init(basename=basename,width=width,height=height)
  par(mgp=c(3,0.3,0))
  # set up plot area
  plot(x=trange,y=val[ind],type='n',xlab="",ylab="",main="",tck=0.04,las=1)
  plotwitherror(x=trange,y=val[ind],dy=dval[ind],rep=TRUE)
  mtext(side=2,line=1.2,text="$D(t)/C(t)$")
  mtext(side=1,line=1.2,text="$t/a$")
  tikz.finalize(tikzfiles)
}
