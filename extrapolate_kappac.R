
extrapolate_kappac <- function(datafile,sim=500,n.predict=5000)
{
  kcdat <- read.table(file=datafile,header=T,stringsAsFactors=FALSE,fill=FALSE)

  kcsimdat <- data.frame(matrix(ncol=nrow(kcdat),nrow=sim))
  for(i in 1:nrow(kcdat)){
    kcsimdat[,i] <- rnorm(n=sim,mean=kcdat$kappac[i],sd=kcdat$dkappac[i])
  }

  kcmodel <- lm(kappac~mu,data=kcdat,weights=1/kcdat$dkappac^2)
  kcsimmodel <- apply(X=kcsimdat,MARGIN=1,
                      FUN=function(x){
                            summary(lm(kappac~mu,data=data.frame(kappac=x,mu=kcdat$mu),weights=1/kcdat$dkappac^2))$coefficients[1:2]
                          }
                      )
  #print(kcmodel)
  #print(kcsimmodel)
  
  kc.cov <- cov(t(kcsimmodel))
  print(kc.cov)

  pred.mu <- data.frame(mu=seq(0.0,0.1,length.out=n.predict))
  pred.kc <- predict(kcmodel,newdata=pred.mu)
  dpar <- as.matrix(data.frame(dc=rep(1,times=n.predict),dmu=pred.mu$mu))
  dy <- sqrt(diag( dpar %*% kc.cov %*% t(dpar) ) )
  
  cat(sprintf("kc(mu=0.0) = %.8f (%.8f)\n",pred.kc[1],dy[1]))
  
  tikzfiles <- tikz.init(basename="kappac_extrap",width=5,height=5,lwdUnit=0.8)
  plotwitherror(x=kcdat$mu,y=kcdat$kappac,dy=kcdat$dkappac,type='n',xlab="$a\\mu$",ylab="$\\kappa_c$",xlim=c(0,0.01))
  lines(y=pred.kc,x=pred.mu$mu)
  polygon(y=c(pred.kc+dy,rev(pred.kc-dy)),x=c(pred.mu$mu,rev(pred.mu$mu)),col="#0000FF55",border=NA)
  plotwitherror(rep=TRUE,x=kcdat$mu,y=kcdat$kappac,dy=kcdat$dkappac)
  tikz.finalize(tikzfiles) 
}
