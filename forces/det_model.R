logseq <- function( d1, d2, length.out) exp(log(10)*seq(d1, d2, length.out=length.out))

det_model <- function(datfile,nsteps,basename="det",n.predict=500,boot.R=300,width=6,height=4.2,kappa2mu=0.0){ 
  fdat <- read.table(header=TRUE,file=datfile,stringsAsFactors=FALSE)
  fdat <- fdat[!(fdat$mon %in% c("GAUGE","cloverdet","cloverdetratio1","cloverdetratio2")),]
  MON <- unique(fdat$mon)
  
  pchlist <- c(15,16,17,18,4,7,8,11)
   
  if(missing(nsteps)) {
    nsteps <- length(which(fdat$mon==MON[1]))
  } else {
    fdat <- fdat[1:(length(MON)*nsteps),]
  }
  
  #require("RColorBrewer")
  
  # extract masses from monomial names
  mulist <- strsplit(MON,"det")
  mulist <- lapply(X=mulist,FUN=function(x){ if(length(x)>1) { 
                                                 as.numeric(paste("0",x[2],sep="."))+kappa2mu
                                               } else { NA } } )
  mus <- NULL
  for( mu in mulist ) {
    if( is.na(mu) ) next
    mus <- c(mus,mu)
  }
  
  require("RColorBrewer")
  clrs <- brewer.pal(n=8,name="Dark2")
  mus <- unique(mus)
  mucols <- data.frame(clr=brewer.pal(n=length(mus),name="Dark2"),mu=mus,stringsAsFactors=FALSE)
  mupch <- data.frame(pch=pchlist,mu=mus)

  # assemble data into data frame with some statistical analysis for confidence intervals
  force.df <- NULL

  df.tsboot <- list()
  for(i in 1:boot.R){
    df.tsboot[[i]] <- data.frame()
  }

  for( i in 1:length(MON) ) {
    ind <- which(fdat$mon==MON[i])
    nsteps <- length(fdat$aver[ind])
    aver.tsboot <- tsboot(tseries=fdat$aver[ind],R=boot.R,l=2,sim="geom",statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) })
    max.tsboot <- tsboot(tseries=fdat$max[ind],R=boot.R,l=2,sim="geom",statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) })
    aver <- c(mean(aver.tsboot$t[,1]),mean(aver.tsboot$t[,1]-aver.tsboot$t[,2]),mean(aver.tsboot$t[,3]-aver.tsboot$t[,1]))
    maxi  <- c(mean(max.tsboot$t[,1]),mean(max.tsboot$t[,1]-max.tsboot$t[,2]),mean(max.tsboot$t[,3]-max.tsboot$t[,1]))
    force.df <- rbind(force.df,
                      data.frame(name=MON[i],
                                 aver=aver[1],maver=aver[2],paver=aver[3],
                                 max=maxi[1],mmax=maxi[2],pmax=maxi[3],
                                 mu=mus[i],
                                 mucol=mucols$clr[which(mucols$mu==mus[i])[1]],
                                 mupch=mupch$pch[which(mupch$mu==mus[i])[1]],
                                 stringsAsFactors=FALSE)
                     )
  }
  options(width=250)
  print(force.df)
  
  fmodel <- list()
  fmodel[[1]] <- nls(aver~a/(mu+b),data=force.df,start=list(a=2.0,b=0.1),
                     trace=TRUE,algorithm='port',control=list(maxiter=1000),
                     weights=1/(force.df$paver^2+force.df$maver^2))
  fmodel[[2]] <- nls(max~a/(mu+b),data=force.df,start=list(a=4.0,b=0.1),
                     trace=TRUE,control=list(maxiter=1000),
                     weights=1/(force.df$pmax^2+force.df$mmax^2))
  #fmodel[[3]] <- nls(pmax~a/abs(mu+b),data=force.df,start=list(a=2.2,b=0.1),#,c=0.1),
  #                   trace=TRUE,algorithm='port',control=list(maxiter=1000),
  #                   weights=1/(force.df$pmax^2+force.df$mmax^2))
  
  for(mod in fmodel){
    print(summary(mod))
  }
  #print(summary(fmodel.pmax))
  save(fmodel,file=sprintf("fmodel.%s.Rdata",basename))
 
  predmu <- data.frame(mu=logseq(-4,0,length.out=n.predict))
  
  tikzfiles <- tikz.init(basename=sprintf("%s.fmodel",basename),width=width,height=height,lwdUnit=0.8)
  par(mgp=c(3,0.3,0))

  ylims <- c(4,2*10^3)
  xlims <- c(8e-4,0.9)
  plotwitherror(x=force.df$mu,y=force.df$aver,main="",tck=0.02,xlim=xlims,
                ylab="$\\|F\\|^2$",xlab="",ylim=ylims,las=1,yaxt='n',xaxt='n',log='yx',type='n')
  for(i in 1:length(fmodel)){
    pred.y <- predict(fmodel[[i]],newdata=predmu)
    lines(y=pred.y,x=predmu$mu,col=clrs[i],lty=i)
  }
  plotwitherror(x=force.df$mu,y=force.df$aver,dy=force.df$paver,mdy=force.df$maver,pch=16,rep=TRUE,col=clrs[1])#col=force.df$mu1col,rep=TRUE)
  plotwitherror(x=force.df$mu,y=force.df$max,dy=force.df$pmax,mdy=force.df$mmax,pch=17,rep=TRUE,col=clrs[2])#col=force.df$mu1col,rep=TRUE)
  plotwitherror(x=force.df$mu,y=force.df$pmax,pch=15,rep=TRUE,col=clrs[3])#col=force.df$mu1col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(0:4),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(0:4)),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-4:-1),tck=0.04,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,10^(-4:0))),tck=0.02)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}$")
  legend("topright",pch=c(17,15,16),col=clrs[c(2,3,1)],bty='n',pt.cex=1.3,
         legend=c("$\\|F\\|^2_\\mathrm{max}$","$\\Delta^+\\left(\\|F\\|^2_\\mathrm{max}\\right)$",
                  "$\\|F\\|^2_\\mathrm{av}$"))
  #legend(x=0.23,y=10^(ylims[1]+1),legend=c("max","aver"),pch=c(17,16),bty='n',lty=c(2,1))
  # to get colours and numbers in the same order
  #legend(x=0.23,y=10^ylims[2],legend=sprintf("$\\mu=%.7f$",rev(mucols$mu)),col=rev(mucols$clr),bty='n',pch=rev(mupch)
 
  tikz.finalize(tikzfiles)
  stop()
  
  ylims <- c(-3,3)
  xlims <- c(4e-3,0.21)
  plotwitherror(x=force.df$mu2,y=force.df$max,main="",tck=0.02,xlim=xlims,
                ylab="$\\|F\\|^2_\\mathrm{max}$",xlab="",ylim=10^ylims,las=1,yaxt='n',xaxt='n',log='xy',type='n')
  for(i in 1:length(constmu1)){
    pred.y <- predict(fmodel[[1]],newdata=constmu1[[i]])
    lines(y=pred.y,x=constmu1[[i]]$mu2,col=mu1cols$clr[i],lty=1)
  }
  plotwitherror(rep=TRUE,x=force.df$mu2,y=force.df$max,dy=force.df$pmax,mdy=force.df$mmax,pch=17,col=force.df$mu1col)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-2:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(4:9,1e-3),outer(2:9,1e-2),2e-1),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_2$")
  #legend(x=0.23,y=10^(ylims[1]+1),legend=c("max","aver"),pch=c(17,16),bty='n',lty=c(2,1))
  # to get colours and numbers in the same order
  legend(x=0.23,y=10^ylims[2],legend=sprintf("$\\mu_1=%.7f$",rev(mu1cols$mu1)),fill=rev(mu1cols$clr),bty='n')

  #####
  
  xlims <- c(1e-3,0.15)
  ylims <- c(-5,3)
  plotwitherror(x=force.df$mu1,y=force.df$aver,main="",tck=0.02,
                xlim=xlims,
                ylab="$\\|F\\|_\\mathrm{av}^2$",xlab="",
                ylim=10^ylims,las=1,yaxt='n',xaxt='n',
                log='xy',
                type='n')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel[[1]],newdata=constmu2[[i]])
    lines(y=pred.y,x=constmu2[[i]]$mu1,col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(x=force.df$mu1,y=force.df$aver,dy=force.df$paver,mdy=force.df$maver,pch=force.df$mu2pch,col=force.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-3:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,c(1e-3,1e-2))),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1$")
  # to get colours and numbers in the same order
  legend("topright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),col=mu2cols$clr,pch=mu2pch$pch,bty='n',cex=0.7,pt.cex=1.1)

  xlims <- c(1e-3,0.15)
  ylims <- c(-3,4)
  plotwitherror(x=force.df$mu1,y=force.df$max,main="",tck=0.02,xlim=xlims,
                ylab="$\\|F\\|_\\mathrm{max}^2$",xlab="",ylim=10^ylims,las=1,yaxt='n',xaxt='n',log='xy',type='n')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel[[2]],newdata=constmu2[[i]])
    lines(y=pred.y,x=constmu2[[i]]$mu1,col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(rep=TRUE,x=force.df$mu1,y=force.df$max,dy=force.df$pmax,mdy=force.df$mmax,pch=force.df$mu2pch,col=force.df$mu2col)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-3:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,c(1e-3,1e-2))),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1$")
  # to get colours and numbers in the same order
  legend("topright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),col=mu2cols$clr,pch=mu2pch$pch,bty='n',cex=0.7,pt.cex=1.1)

  tikz.finalize(tikzfiles) 
  
  #########################################

  tikzfiles <- tikz.init(basename=sprintf("%s.fmodel.ratios",basename),width=5,height=3.4,lwdUnit=0.8)
  par(mgp=c(3,0.3,0),mar=c(5,4,4,2))

  # ratios of maximum and average forces 
  constmu1 <- list()
  constmu2 <- list()
  for(mu1 in mu1s){
    constmu1[[length(constmu1)+1]] <- data.frame(mu1=rep(mu1,times=n.predict),mu2=seq(-0.05,10000*mu1,length.out=n.predict))
  }
  for(mu2 in mu2s){
    constmu2[[length(constmu2)+1]] <- data.frame(mu1=seq(mu2/10000,100*mu2,length.out=n.predict),mu2=rep(mu2,times=n.predict))
  }
  
  xlims <- c(-3,3)
  ylims <- c(1,3)
  plotwitherror(x=(force.df$mu2/force.df$mu1)^(-1),y=force.df$aver,main="",tck=0.02,xlim=10^c(xlims[1],(xlims[2]+0.4)),
                ylab="$\\|F\\|_\\mathrm{max}^2/\\|F\\|_\\mathrm{av}^2$",
                xlab="",ylim=10^c(ylims[1],(ylims[2]+0.2)),las=1,yaxt='n',xaxt='n',type='n',log='xy')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel[[2]],newdata=constmu2[[i]])
    pred.y <- pred.y/predict(fmodel[[1]],newdata=constmu2[[i]])
    lines(y=pred.y,x=(constmu2[[i]]$mu2/constmu2[[i]]$mu1)^(-1),col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(x=(force.df$mu2/force.df$mu1)^(-1),y=force.df$max/force.df$aver,
                dy=sqrt( (force.df$pmax/force.df$aver)^2 + (force.df$paver*force.df$max/force.df$aver^2)^2),
                mdy=sqrt( (force.df$mmax/force.df$aver)^2 + (force.df$maver*force.df$max/force.df$aver^2)^2),
                pch=force.df$mu2pch,col=force.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.03,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.015)
  axis(side=1, labels=TRUE, at=10^(xlims[1]:xlims[2]),tck=0.03,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,10^(xlims[1]:xlims[2]))),tck=0.015)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1/\\tilde{\\mu}_2$")
  # to get colours and numbers in the same order
  legend("topright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",rev(mu2cols$mu2)),col=rev(mu2cols$clr),pch=rev(mu2pch$pch),bty='n',cex=0.7,pt.cex=1.1)
  
  xlims <- c(-3,2)
  ylims <- c(-3,2)
  plotwitherror(x=(force.df$mu2/force.df$mu1)^(-1),y=force.df$pmax,main="",tck=0.02,
                #xlim=(10^xlims),
                ylim=10^c(ylims[1],(ylims[2]+0.2)),
                ylab="$\\Delta^+(\\|F\\|_\\mathrm{max}^2)$",
                xlab="",las=1,yaxt='n',xaxt='n',type='n',log='xy') 
  plotwitherror(x=(force.df$mu2/force.df$mu1)^(-1),y=force.df$pmax,
                pch=force.df$pch,col=force.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.03,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.015)
  axis(side=1, labels=TRUE, at=10^(xlims[1]:xlims[2]),tck=0.03,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,10^(xlims[1]:xlims[2]))),tck=0.015)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1/\\tilde{\\mu}_2$")
  # to get colours and numbers in the same order
  legend("bottomleft",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",rev(mu2cols$mu2)),col=rev(mu2cols$clr),pch=rev(mu2pch$pch),bty='n',cex=0.7,pt.cex=1.1)
  
  # return to original parameters 
  constmu1 <- list()
  constmu2 <- list()
  for(mu1 in mu1s){
    constmu1[[length(constmu1)+1]] <- data.frame(mu1=rep(mu1,times=n.predict),mu2=logseq(-4,-0.3,length.out=n.predict))
  }
  for(mu2 in mu2s){
    constmu2[[length(constmu2)+1]] <- data.frame(mu1=logseq(-4,-0.3,length.out=n.predict),mu2=rep(mu2,times=n.predict))
  }
  
  ylims <- c(-5,2)
  plotwitherror(x=force.df$mu1,y=force.df$pmax,main="",tck=0.02,
                #xlim=(10^xlims),
                ylim=10^c(ylims[1],(ylims[2]+0.2)),
                ylab="$\\Delta^+(\\|F\\|_\\mathrm{max}^2)$",
                xlab="",las=1,yaxt='n',xaxt='n',type='n',log='xy')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel.pmax,newdata=constmu2[[i]])
    lines(y=pred.y,x=constmu2[[i]]$mu1,col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(x=force.df$mu1,y=force.df$pmax,
                pch=force.df$mu2pch,col=force.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.03,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.015)
  axis(side=1, labels=TRUE, at=10^(xlims[1]:xlims[2]),tck=0.03,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,10^(xlims[1]:xlims[2]))),tck=0.015)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1$")
  # to get colours and numbers in the same order
  legend("bottomleft",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),col=mu2cols$clr,bty='n',cex=0.7,pt.cex=1.1,pch=mu2pch$pch)
  
  tikz.finalize(tikzfiles)
  stop()
  
  fmodel.dmax <- nls((max-aver)~a*abs(mu2-mu1)^b*abs(mu2/mu1)^c,data=force.df,start=list(a=0.1,b=0.2,c=0.3),trace=TRUE,algorithm='port',control=list(maxiter=1000))

  plotwitherror(x=force.df$mu2,y=force.df$max-force.df$aver,dy=sqrt(force.df$paver^2+force.df$pmax^2),mdy=sqrt(force.df$maver^2+force.df$mmax^2),
                xlab="$\\mu_2$",ylab="$\\max(F^2)-\\bar{F}^2$",tck=0.02,log='y',col=force.df$mu1col,pch=16,las=1)
  for(i in 1:length(constmu1)){
    pred.y <- predict(fmodel.dmax,newdata=constmu1[[i]])
    lines(y=pred.y,x=constmu1[[i]]$mu2,col=mu1cols$clr[i],lty=j)
  }
  # to get colours and numbers in the same order
  legend("bottomright",legend=sprintf("$\\mu_1=%.7f$",rev(mu1cols$mu1)),fill=rev(mu1cols$clr),bty='n')
  
  plotwitherror(x=force.df$mu1,y=force.df$max-force.df$aver,dy=sqrt(force.df$paver^2+force.df$pmax^2),mdy=sqrt(force.df$maver^2+force.df$mmax^2),
                xlab="$\\mu_1$",ylab="$\\max(F^2)-\\bar{F}^2$",tck=0.02,log='y',col=force.df$mu2col,pch=16,las=1,xlim=c(0,0.2))
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel.dmax,newdata=constmu2[[i]])
    lines(y=pred.y,x=constmu2[[i]]$mu1,col=mu2cols$clr[i],lty=j)
  }
  # to get colours and numbers in the same order
  legend("bottomright",legend=sprintf("$\\mu_2=%.7f$",mu2cols$mu2),fill=mu2cols$clr,bty='n')

  tikz.finalize(tikzfiles)
}


