require("propagate")
require("RColorBrewer")

overlay_massprec <- function(datfiles,mus,nsteps,lg,basename,n.predict=500,boot.R=300,width=5,height=3.3) {  
  force.df <- NULL
  ensclrs <- brewer.pal(n=length(datfiles),name="Dark2")
  for(fc in 1:length(datfiles)){
    fdat <- read.table(header=TRUE,file=datfiles[fc],stringsAsFactors=FALSE)
    fdat <- fdat[!fdat$mon%in% c("GAUGE","cloverdet","cloverdetlight","ndcloverrat1a","ndcloverrat2a","ndcloverrat3a"),]
    MON <- unique(fdat$mon)
    
    if(missing(nsteps)) {
      nsteps <- length(which(fdat$mon==MON[1]))
    } else {
      fdat <- fdat[1:(length(MON)*nsteps),]
    }
    
    #require("RColorBrewer")
    #mu1s <- unique(mus$mu1)
    #mu2s <- unique(mus$mu2)
    #mu1cols <- data.frame(clr=brewer.pal(n=length(mu1s),name="Dark2"),mu1=mu1s,stringsAsFactors=FALSE)
    #mu2cols <- data.frame(clr=brewer.pal(n=length(mu2s),name="Dark2"),mu2=mu2s,stringsAsFactors=FALSE)
    
    # assemble data into data frame with some statistical analysis for confidence intervals
    for( i in 1:length(MON) ) {
      cat(sprintf("Analyzing %s\n",MON[i]))
      ind <- which(fdat$mon==MON[i])
      nsteps <- length(fdat$aver[ind])
      aver.tsboot <- tsboot(tseries=fdat$aver[ind],R=boot.R,l=2,sim="geom",statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) })
      max.tsboot <- tsboot(tseries=fdat$max[ind],R=boot.R,l=2,sim="geom",statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) })
      aver <- c(mean(aver.tsboot$t[,1]),mean(aver.tsboot$t[,1]-aver.tsboot$t[,2]),mean(aver.tsboot$t[,3]-aver.tsboot$t[,1]))
      maxi  <- c(mean(max.tsboot$t[,1]),mean(max.tsboot$t[,1]-max.tsboot$t[,2]),mean(max.tsboot$t[,3]-max.tsboot$t[,1]))
      force.df <- rbind(force.df,
                        data.frame(aver=aver[1],maver=aver[2],paver=aver[3],
                                   max=maxi[1],mmax=maxi[2],pmax=maxi[3],
                                   mu1=mus[[fc]]$mu1[i],mu2=mus[[fc]]$mu2[i],kappa2mu=mus[[fc]]$kappa2mu[i],
                                   mu1col=ensclrs[fc],
                                   mu2col=ensclrs[fc],
                                   stringsAsFactors=FALSE)
                       )
    }
  }

  mu1s <- unique(force.df$mu1)
  mu2s <- unique(force.df$mu2)
  #require("RColorBrewer")
  #mu1cols <- data.frame(clr=brewer.pal(n=length(mu1s),name="Spectral"),mu1=mu1s,stringsAsFactors=FALSE)
  #mu2cols <- data.frame(clr=brewer.pal(n=length(mu2s),name="Spectral"),mu2=mu2s,stringsAsFactors=FALSE)

  #force.df <- cbind(force.df,mu1col=c(NA),mu2col=c(NA))
  #for( i in 1:nrow(force.df) ){
  #  force.df$mu1col[i] <- mu1cols$clr[which(mu1cols$mu1==force.df$mu1[i])[1]]
  #  force.df$mu2col[i] <- mu2cols$clr[which(mu2cols$mu2==force.df$mu2[i])[1]]
  #}
  options(width=250)
  print(force.df)

  constmu1 <- list()
  constmu2 <- list()
  for(mu1 in mu1s){
    constmu1[[length(constmu1)+1]] <- data.frame(mu1=rep(mu1,times=n.predict),mu2=seq(0,1.0,length.out=n.predict))
  }
  for(mu2 in mu2s){
    constmu2[[length(constmu2)+1]] <- data.frame(mu1=seq(0,0.6,length.out=n.predict),mu2=rep(mu2,times=n.predict))
  }
  
  load(sprintf("fmodel.%s.Rdata",basename))

  #cA2.30.24 <- data.frame(rho1=c(0.0363,0.0)+0.0008237,rho2=c(0.04,0.008)+0.0008238,clr=rep("blue",2),stringsAsFactors=FALSE)
  #cA2.60.24 <- data.frame(rho1=c(0.011,0.0)+0.0016476,rho2=c(0.06,0.011)+0.0016476,clr=rep("red",2),stringsAsFactors=FALSE)
  #cA2.60.32 <- data.frame(rho1=c(0.08,0.008,0.0)+0.0016476,rho2=c(0.8,0.08,0.008)+0.0016476,clr=rep("magenta",3),stringsAsFactors=FALSE)
  #typ.data <- rbind(cA2.30.24,cA2.60.24,cA2.60.32)
  #typ.y <- predict(rhomodel,newdata=typ.data)

  #rho <- data.frame(rho1=rep(0.0008237+0.014,times=100),rho2=seq(0,0.2,length.out=100))  ,
  #print(cbind(predict(rhomodel,newdata=rho),rho))
  #readline()

  tikzfiles <- tikz.init(basename=sprintf("%s.overlay",basename),width=width,height=height,lwdUnit=0.8)
  par(mgp=c(3,0.3,0))
  ylims <- c(-5,2)
  xlims <- c(-3,-0.7)
  plotwitherror(x=force.df$mu2,y=force.df$aver,main="",tck=0.02,xlim=10^xlims,
                ylab="$\\|F\\|^2$",xlab="",ylim=10^ylims,las=1,yaxt='n',xaxt='n',log='xy',type='n')
  for(j in 1:length(fmodel)){
    pred.y <- predictNLS(fmodel[[j]],newdata=data.frame(mu2=force.df$mu2,mu1=force.df$mu1),alpha=0.3146,interval="confidence",do.sim=TRUE,second.order=FALSE)
    plotwitherror(x=force.df$mu2,y=pred.y$summary[,9],mdy=pred.y$summary[,9]-pred.y$summary[,11],dy=pred.y$summary[,12]-pred.y$summary[,9],pch=j,col=force.df$mu1col,rep=TRUE)
  }
  plotwitherror(x=force.df$mu2,y=force.df$aver,dy=force.df$paver,mdy=force.df$maver,pch=16,col=force.df$mu1col,rep=TRUE)
  plotwitherror(rep=TRUE,x=force.df$mu2,y=force.df$max,dy=force.df$pmax,mdy=force.df$mmax,pch=17,col=force.df$mu1col)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-3:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=outer(2:9,10^(-3:-2)),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_2$")
  legend("topright",legend=c("max","aver"),pch=c(17,16),bty='n')
  legend("topleft",horiz=TRUE,legend=c("prediction","measurement"),pch=c(0,15),bty='n')
  if(!missing(lg))
    legend("bottomright",legend=sprintf("\\textit{%s}",lg),fill=ensclrs,bty='n',cex=0.8)

  #####

  plotwitherror(x=force.df$mu1,y=force.df$aver,main="",tck=0.02,xlim=c(2e-4,0.02),
                ylab="$\\|F\\|^2$",xlab="",ylim=10^ylims,las=1,yaxt='n',log='xy',type='n')
  for(j in 1:length(fmodel)){
    pred.y <- predict(fmodel[[j]],newdata=force.df)
    points(x=force.df$mu1,y=pred.y,pch=j,col=force.df$mu2col)
  }
  plotwitherror(x=force.df$mu1,y=force.df$aver,dy=force.df$paver,mdy=force.df$maver,pch=16,col=force.df$mu2col,rep=TRUE)
  plotwitherror(rep=TRUE,x=force.df$mu1,y=force.df$max,dy=force.df$pmax,mdy=force.df$mmax,pch=17,col=force.df$mu2col)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1$")
  legend("topright",legend=c("max","aver"),pch=c(17,16),pt.cex=1.3,bty='n')
  legend("topleft",horiz=TRUE,legend=c("prediction","measurement"),pch=c(0,15),bty='n')
  if(!missing(lg))
    legend("bottomright",legend=sprintf("\\textit{%s}",lg),fill=ensclrs,bty='n')
  
  ylims <- c(1,5) 
  plotwitherror(x=force.df$mu2,y=force.df$max/force.df$aver,main="",tck=0.02,xlim=10^xlims,
                ylab="$\\|F\\|_\\mathrm{max}^2/\\|F\\|_\\mathrm{av}^2$",xlab="",ylim=10^ylims,las=1,yaxt='n',xaxt='n',log='xy',type='n')
  pred.y <- predict(fmodel[[2]],newdata=force.df)
  pred.y <- pred.y/predict(fmodel[[1]],newdata=force.df)
  points(x=force.df$mu2,y=pred.y,pch=1,col=force.df$mu2col)
  rat <- force.df$max/force.df$aver
  pdrat <- sqrt( (force.df$max*force.df$paver/force.df$aver^2)^2 + (force.df$pmax/force.df$aver)^2 )
  mdrat <- sqrt( (force.df$max*force.df$maver/force.df$aver^2)^2 + (force.df$mmax/force.df$aver)^2 )
  plotwitherror(x=force.df$mu2,y=rat,dy=pdrat,mdy=mdrat,pch=16,col=force.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-3:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=outer(2:9,10^(-3:-2)),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_2$")
  #legend("topright",legend=c("max","aver"),pch=c(17,16),pt.cex=1.3,bty='n')
  legend("topleft",horiz=TRUE,legend=c("prediction","measurement"),pch=c(0,15),bty='n')
  if(!missing(lg))
    legend("topright",legend=sprintf("\\textit{%s}",lg),fill=ensclrs,bty='n',cex=0.8)

  tikz.finalize(tikzfiles)
  stop("stopping here, further code is unfinished!\n")
  fmodel.dmax <- nls((max-aver)~a*(mu2-mu1)^b*(mu2/mu1)^c,data=force.df,start=list(a=0.1,b=0.2,c=0.3),trace=TRUE,algorithm='port',control=list(maxiter=1000))

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


