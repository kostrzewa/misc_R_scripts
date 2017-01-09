plot_P_mudelsq <- function(datafile,basename="P_vs_mudelsq",n.sim=500,width=3.2,height=3.2,n.predict=500) {
  Pdat <- read.table(file=datafile,header=TRUE,stringsAsFactors=FALSE)
  Pdat <- data.frame(Pdat,mudelsq=Pdat$mudel^2)

  Pfitdat <- Pdat[which(Pdat$mudel>0.0),]

  simdat <- NULL
  if(n.sim > 0) {
    simdat <- data.frame(matrix(ncol=nrow(Pfitdat),nrow=n.sim))
    for( k in 1:nrow(Pfitdat) ){
      simdat[,k] <- rnorm(n=n.sim,mean=Pfitdat$P[k],sd=Pfitdat$dP[k])
    }
  }

  P.mod <- lm(P~mudelsq,data=Pfitdat,weights=1/Pfitdat$dP^2)
  Psim.mod <- apply(X=simdat,MARGIN=1,
                   FUN=function(x){
                         lm(P~mudelsq,data=data.frame(P=x,mudelsq=Pfitdat$mudelsq),weights=1/Pfitdat$dP^2)$coefficients
                       }
                  )

  P.cov <- cov(t(Psim.mod))
  print(summary(P.mod))
  cat("Errors on fit coefficients:\n")
  print(sqrt(diag(P.cov)))
  
  newdat <- data.frame(mudelsq=seq(0.0,0.2,length.out=n.predict))
  P.pred <- predict(P.mod,newdata=newdat)
  dPdc <- as.matrix(cbind(rep(1,n.predict),newdat$mudelsq))
  
  dP <- sqrt(diag(dPdc %*% P.cov %*% t(dPdc)))
  poly.x <- c(newdat$mudelsq,rev(newdat$mudelsq))
  poly.y <- c(P.pred+dP,rev(P.pred-dP))

  Prange <- range(Pdat$P)+c(-0.00005,0.00005)
  tikzfiles <- tikz.init(basename=basename,width=width,height=height,lwdUnit=0.8)
  par(mgp=c(3,0.3,0))
  # prep plot area
  plot(y=Pdat$P,x=Pdat$mudelsq,type='n',main="",
       xlab="",ylab="$\\langle P \\rangle$",las=1,tck=0.02,
       ylim=Prange)
  mtext(text="$(a\\mu_\\delta)^2$",side=1,line=1.3)
  polygon(x=poly.x,y=poly.y,col="#1B9E7788",border=NA)
  lines(y=P.pred,x=newdat$mudelsq)
  plotwitherror(rep=TRUE,x=Pdat$mudelsq,y=Pdat$P,dy=Pdat$dP,col=Pdat$clr,pch=Pdat$pch)
  legend("topleft",pch=c(15,0,1),x.intersp=0.02,cex=0.95,
         legend=c(sprintf("$\\langle P \\rangle = %.5f + %.3f (a\\mu_\\delta)^2$",P.mod$coefficients[1],P.mod$coefficients[2]),"RHMC","HMC"),
         bty='n',col=c("#1B9E7788","blue","black"),lty=c(1,NA,NA))
  tikz.finalize(tikzfiles)
  

}
