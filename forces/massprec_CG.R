logseq <- function( d1, d2, length.out) exp(log(10)*seq(d1, d2, length.out=length.out))

massprec_CG <- function(datfile,nsteps,basename="CG",n.predict=500,kappa2mu=0.0002746,boot.R=300,width=6,height=4.2){ 
  fdat <- read.table(header=TRUE,file=datfile,stringsAsFactors=FALSE)
  MON <- unique(fdat$mon)
 
  if(missing(nsteps)) {
    nsteps <- length(which(fdat$mon==MON[1]))
  } else {
    fdat <- fdat[1:(length(MON)*nsteps),]
  }
  
  #require("RColorBrewer")
  
  # extract masses from monomial names
  rholist <- strsplit(MON,"rho")
  rholist <- lapply(X=rholist,FUN=function(x){ if(length(x)>1) { 
                                                 data.frame(mu1=as.numeric(paste("0",x[2],sep=".")),
                                                            mu2=as.numeric(paste("0",x[3],sep="."))) 
                                               } else { NA } } )
  mus <- NULL
  for( rho in rholist ) {
    if( is.na(rho) ) next
    mus <- rbind(mus,rho+kappa2mu)
  }
 
  require("RColorBrewer")
  mu1s <- unique(mus$mu1)
  mu2s <- unique(mus$mu2)
  mu1cols <- data.frame(clr=rainbow(n=length(mu1s)),mu1=mu1s,stringsAsFactors=FALSE)
  mu2cols <- data.frame(clr=brewer.pal(n=length(mu2s),name="Dark2"),mu2=mu2s,stringsAsFactors=FALSE)
  
  # assemble data into data frame with some statistical analysis for confidence intervals
  cg.df <- NULL

  for( i in 1:length(MON) ) {
    ind <- which(fdat$mon==MON[i])
    nsteps <- length(fdat$iter[ind])
    cg.tsboot <- tsboot(tseries=fdat$iter[ind],R=boot.R,l=2,sim="geom",statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) })
    cg <- c(mean(cg.tsboot$t[,1]),mean(cg.tsboot$t[,1]-cg.tsboot$t[,2]),mean(cg.tsboot$t[,3]-cg.tsboot$t[,1]))
    cg.df <- rbind(cg.df,
                    data.frame(name=MON[i],
                               cg=cg[1],mdcg=cg[2],dcg=cg[3],
                               mu1=mus$mu1[i],mu2=mus$mu2[i],
                               mu1col=mu1cols$clr[which(mu1cols$mu1==mus$mu1[i])[1]],
                               mu2col=mu2cols$clr[which(mu2cols$mu2==mus$mu2[i])[1]],
                               stringsAsFactors=FALSE)
                     )
  }
  options(width=250)
  print(cg.df)

  
  fmodel <- nls(cg~a/(mu1+b),data=cg.df,start=list(a=0.1,b=1),trace=TRUE,algorithm='port',control=list(maxiter=1000))
  print(summary(fmodel))
  #save(fmodel,file=sprintf("fmodel.%s.Rdata",basename))
 
  constmu1 <- list()
  constmu2 <- list()
  for(mu1 in mu1s){
    constmu1[[length(constmu1)+1]] <- data.frame(mu1=rep(mu1,times=n.predict),mu2=logseq(-4,-0.3,length.out=n.predict))
  }
  for(mu2 in mu2s){
    constmu2[[length(constmu2)+1]] <- data.frame(mu1=logseq(-4,-0.3,length.out=n.predict),mu2=rep(mu2,times=n.predict))
  }
 
  tikzfiles <- tikz.init(basename=sprintf("%s.massprec",basename),width=width,height=height,lwdUnit=0.8)
  par(mgp=c(3,0.3,0))
  
  ylims <- c(0,3)
  xlims <- c(4e-3,0.21)
  plotwitherror(x=cg.df$mu2,y=cg.df$cg,main="",tck=0.02,
                #xlim=xlims,
                ylab="$\\|F\\|^2_\\mathrm{av}$",xlab="",
                #ylim=10^ylims,
                las=1,yaxt='n',xaxt='n',
                log='xy',
                type='n')

  for(i in 1:length(constmu1)){
    pred.y <- predict(fmodel,newdata=constmu1[[i]])
    lines(y=pred.y,x=constmu1[[i]]$mu2,col=mu1cols$clr[i],lty=1)
  }
  plotwitherror(x=cg.df$mu2,y=cg.df$cg,dy=cg.df$dcg,mdy=cg.df$mdcg,pch=16,col=cg.df$mu1col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-2:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(4:9,1e-3),outer(2:9,1e-2),2e-1),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_2$")
  #legend(x=0.23,y=10^(ylims[1]+1),legend=c("max","aver"),pch=c(17,16),bty='n',lty=c(2,1))
  # to get colours and numbers in the same order
  #legend(x=0.23,y=10^ylims[2],legend=sprintf("$\\mu_1=%.7f$",rev(mu1cols$mu1)),fill=rev(mu1cols$clr),bty='n')
  
  plotwitherror(x=cg.df$mu1,y=cg.df$cg,main="",tck=0.02,
                #xlim=xlims,
                ylab="$\\|F\\|^2_\\mathrm{av}$",xlab="",
                #ylim=10^ylims,
                las=1,yaxt='n',xaxt='n',
                #log='x',
                type='n')

  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel,newdata=constmu2[[i]])
    lines(y=pred.y,x=constmu2[[i]]$mu1,col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(x=cg.df$mu1,y=cg.df$cg,dy=cg.df$dcg,mdy=cg.df$mdcg,pch=16,col=cg.df$mu1col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-2:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(4:9,1e-3),outer(2:9,1e-2),2e-1),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1$")
  #legend(x=0.23,y=10^(ylims[1]+1),legend=c("max","aver"),pch=c(17,16),bty='n',lty=c(2,1))
  # to get colours and numbers in the same order
  #legend(x=0.23,y=10^ylims[2],legend=sprintf("$\\mu_1=%.7f$",rev(mu1cols$mu1)),fill=rev(mu1cols$clr),bty='n')
  
  tikz.finalize(tikzfiles)
  stop()

  ylims <- c(-3,3)
  xlims <- c(4e-3,0.21)
  plotwitherror(x=cg.df$mu2,y=cg.df$max,main="",tck=0.02,xlim=xlims,
                ylab="$\\|F\\|^2_\\mathrm{max}$",xlab="",ylim=10^ylims,las=1,yaxt='n',xaxt='n',log='xy',type='n')
  for(i in 1:length(constmu1)){
    pred.y <- predict(fmodel[[1]],newdata=constmu1[[i]])
    lines(y=pred.y,x=constmu1[[i]]$mu2,col=mu1cols$clr[i],lty=1)
  }
  plotwitherror(rep=TRUE,x=cg.df$mu2,y=cg.df$max,dy=cg.df$pmax,mdy=cg.df$mmax,pch=17,col=cg.df$mu1col)
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
  plotwitherror(x=cg.df$mu1,y=cg.df$cg,main="",tck=0.02,xlim=xlims,
                ylab="$\\|F\\|_\\mathrm{av}^2$",xlab="",ylim=10^ylims,las=1,yaxt='n',xaxt='n',log='xy',type='n')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel[[1]],newdata=constmu2[[i]])
    lines(y=pred.y,x=constmu2[[i]]$mu1,col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(x=cg.df$mu1,y=cg.df$cg,dy=cg.df$dcg,mdy=cg.df$mdcg,pch=16,col=cg.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-3:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,c(1e-3,1e-2))),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1$")
  # to get colours and numbers in the same order
  legend("topright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),fill=mu2cols$clr,bty='n',cex=0.7)

  xlims <- c(1e-3,0.15)
  ylims <- c(-3,4)
  plotwitherror(x=cg.df$mu1,y=cg.df$max,main="",tck=0.02,xlim=xlims,
                ylab="$\\|F\\|_\\mathrm{max}^2$",xlab="",ylim=10^ylims,las=1,yaxt='n',xaxt='n',log='xy',type='n')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel[[2]],newdata=constmu2[[i]])
    lines(y=pred.y,x=constmu2[[i]]$mu1,col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(rep=TRUE,x=cg.df$mu1,y=cg.df$max,dy=cg.df$pmax,mdy=cg.df$mmax,pch=17,col=cg.df$mu2col)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-3:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,c(1e-3,1e-2))),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1$")
  # to get colours and numbers in the same order
  legend("topright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),fill=mu2cols$clr,bty='n',cex=0.7)

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
  plotwitherror(x=(cg.df$mu2/cg.df$mu1)^(-1),y=cg.df$cg,main="",tck=0.02,xlim=10^c(xlims[1],(xlims[2]+0.4)),
                ylab="$\\|F\\|_\\mathrm{max}^2/\\|F\\|_\\mathrm{av}^2$",
                xlab="",ylim=10^c(ylims[1],(ylims[2]+0.2)),las=1,yaxt='n',xaxt='n',type='n',log='xy')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel[[2]],newdata=constmu2[[i]])
    pred.y <- pred.y/predict(fmodel[[1]],newdata=constmu2[[i]])
    lines(y=pred.y,x=(constmu2[[i]]$mu2/constmu2[[i]]$mu1)^(-1),col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(x=(cg.df$mu2/cg.df$mu1)^(-1),y=cg.df$max/cg.df$cg,
                dy=sqrt( (cg.df$pmax/cg.df$cg)^2 + (cg.df$dcg*cg.df$max/cg.df$cg^2)^2),
                mdy=sqrt( (cg.df$mmax/cg.df$cg)^2 + (cg.df$mdcg*cg.df$max/cg.df$cg^2)^2),
                pch=16,col=cg.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.03,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.015)
  axis(side=1, labels=TRUE, at=10^(xlims[1]:xlims[2]),tck=0.03,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,10^(xlims[1]:xlims[2]))),tck=0.015)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1/\\tilde{\\mu}_2$")
  # to get colours and numbers in the same order
  legend("topright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",rev(mu2cols$mu2)),fill=rev(mu2cols$clr),bty='n',cex=0.7)
  
  xlims <- c(-3,2)
  ylims <- c(-3,2)
  plotwitherror(x=(cg.df$mu2/cg.df$mu1)^(-1),y=cg.df$pmax,main="",tck=0.02,
                #xlim=(10^xlims),
                ylim=10^c(ylims[1],(ylims[2]+0.2)),
                ylab="$\\Delta^+(\\|F\\|_\\mathrm{max}^2)$",
                xlab="",las=1,yaxt='n',xaxt='n',type='n',log='xy') 
  plotwitherror(x=(cg.df$mu2/cg.df$mu1)^(-1),y=cg.df$pmax,
                pch=16,col=cg.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.03,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.015)
  axis(side=1, labels=TRUE, at=10^(xlims[1]:xlims[2]),tck=0.03,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,10^(xlims[1]:xlims[2]))),tck=0.015)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1/\\tilde{\\mu}_2$")
  # to get colours and numbers in the same order
  legend("bottomleft",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",rev(mu2cols$mu2)),fill=rev(mu2cols$clr),bty='n',cex=0.7)
  
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
  plotwitherror(x=cg.df$mu1,y=cg.df$pmax,main="",tck=0.02,
                #xlim=(10^xlims),
                ylim=10^c(ylims[1],(ylims[2]+0.2)),
                ylab="$\\Delta^+(\\|F\\|_\\mathrm{max}^2)$",
                xlab="",las=1,yaxt='n',xaxt='n',type='n',log='xy')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel.pmax,newdata=constmu2[[i]])
    lines(y=pred.y,x=constmu2[[i]]$mu1,col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(x=cg.df$mu1,y=cg.df$pmax,
                pch=16,col=cg.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.03,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.015)
  axis(side=1, labels=TRUE, at=10^(xlims[1]:xlims[2]),tck=0.03,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,10^(xlims[1]:xlims[2]))),tck=0.015)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_1$")
  # to get colours and numbers in the same order
  legend("bottomleft",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),fill=mu2cols$clr,bty='n',cex=0.7)
  
  tikz.finalize(tikzfiles)
  stop()
  
  fmodel.dmax <- nls((max-aver)~a*abs(mu2-mu1)^b*abs(mu2/mu1)^c,data=cg.df,start=list(a=0.1,b=0.2,c=0.3),trace=TRUE,algorithm='port',control=list(maxiter=1000))

  plotwitherror(x=cg.df$mu2,y=cg.df$max-cg.df$cg,dy=sqrt(cg.df$dcg^2+cg.df$pmax^2),mdy=sqrt(cg.df$mdcg^2+cg.df$mmax^2),
                xlab="$\\mu_2$",ylab="$\\max(F^2)-\\bar{F}^2$",tck=0.02,log='y',col=cg.df$mu1col,pch=16,las=1)
  for(i in 1:length(constmu1)){
    pred.y <- predict(fmodel.dmax,newdata=constmu1[[i]])
    lines(y=pred.y,x=constmu1[[i]]$mu2,col=mu1cols$clr[i],lty=j)
  }
  # to get colours and numbers in the same order
  legend("bottomright",legend=sprintf("$\\mu_1=%.7f$",rev(mu1cols$mu1)),fill=rev(mu1cols$clr),bty='n')
  
  plotwitherror(x=cg.df$mu1,y=cg.df$max-cg.df$cg,dy=sqrt(cg.df$dcg^2+cg.df$pmax^2),mdy=sqrt(cg.df$mdcg^2+cg.df$mmax^2),
                xlab="$\\mu_1$",ylab="$\\max(F^2)-\\bar{F}^2$",tck=0.02,log='y',col=cg.df$mu2col,pch=16,las=1,xlim=c(0,0.2))
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel.dmax,newdata=constmu2[[i]])
    lines(y=pred.y,x=constmu2[[i]]$mu1,col=mu2cols$clr[i],lty=j)
  }
  # to get colours and numbers in the same order
  legend("bottomright",legend=sprintf("$\\mu_2=%.7f$",mu2cols$mu2),fill=mu2cols$clr,bty='n')

  tikz.finalize(tikzfiles)
}


