logseq <- function( d1, d2, length.out) exp(log(10)*seq(d1, d2, length.out=length.out))

massprec_model_sim <- function(datfile,kappa2mu,basename="massprec",n.predict=500,boot.R=300,width=6,height=4.2,fit="pheno"){ 
  fdat <- NULL
  mus <- NULL
  MON <- NULL
  for( i in 1:length(datfile) ){
    ftemp <- read.table(header=TRUE,file=datfile[i],stringsAsFactors=FALSE)
    ftemp <- ftemp[!(ftemp$mon %in% c("GAUGE","cloverdet","cloverdetratio1","cloverdetratio2")),]
    
    montemp <- unique(ftemp$mon)
    rholist <- strsplit(montemp,"rho")
    rholist <- lapply(X=rholist,FUN=function(x){ if(length(x)>1) { 
                                                   data.frame(mu1=as.numeric(paste("0",x[2],sep=".")),
                                                              mu2=as.numeric(paste("0",x[3],sep=".")),
                                                              kappa2mu=0.0) 
                                                 } else { NA } } )
    
    MON <- c(MON,paste(montemp,kappa2mu[i],sep="."))
    ftemp$mon <- paste(ftemp$mon,kappa2mu[i],sep=".")

    for( rho in rholist ) {
      if( is.na(rho) ) next
      mus <- rbind(mus,rho+kappa2mu[i])
    }
    fdat <- rbind(fdat,ftemp)
  }
  
  require("RColorBrewer")
  mu1s <- sort(unique(mus$mu1),decreasing=TRUE)
  mu2s <- unique(mus$mu2)
  #pchlist <- pchlist[1:length(mu2s)]
  # legend scaling depends on the number of mu2s
  cex <- 0.8 
  if(length(mu2s)>6)
    cex <- 0.7 
  mu1cols <- data.frame(clr=rainbow(n=length(mu1s)),mu1=mu1s,stringsAsFactors=FALSE)
  mu2cols <- data.frame(clr=rainbow(n=length(mu2s)),mu2=mu2s,stringsAsFactors=FALSE)
  #mu2cols <- data.frame(clr=brewer.pal(n=length(mu2s),name="Dark2"),mu2=mu2s,stringsAsFactors=FALSE)
  #mu2pch <- data.frame(pch=pchlist,mu2=mu2s)
  mu2pch <- data.frame(pch=16)

  # assemble data into data frame with some statistical analysis for confidence intervals
  force.df <- NULL

  df.tsboot <- list()
  for(i in 1:boot.R){
    df.tsboot[[i]] <- data.frame()
  }

  for( i in 1:length(MON) ) {
    #if(mus$mu1[i] >= mus$mu2[i]) next
    ind <- which(fdat$mon==MON[i])
    nsteps <- length(fdat$aver[ind])
    aver.tsboot <- tsboot(tseries=fdat$aver[ind],R=boot.R,l=2,sim="geom",statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) })
    max.tsboot <- tsboot(tseries=fdat$max[ind],R=boot.R,l=2,sim="geom",statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) })
    # for possibly more involved error analysis
    #for(R in 1:boot.R){
    #  df.tsboot[[R]] <- rbind(df.tsboot[[R]],data.frame(mu1=mus$mu1[i],mu2=mus$mu2[i],aver=aver.tsboot$t[R,1],max=max.tsboot$t[R,1]))
    #}
    aver <- c(mean(aver.tsboot$t[,1]),mean(aver.tsboot$t[,1]-aver.tsboot$t[,2]),mean(aver.tsboot$t[,3]-aver.tsboot$t[,1]))
    maxi  <- c(mean(max.tsboot$t[,1]),mean(max.tsboot$t[,1]-max.tsboot$t[,2]),mean(max.tsboot$t[,3]-max.tsboot$t[,1]))
    force.df <- rbind(force.df,
                      data.frame(name=MON[i],
                                 aver=aver[1],maver=aver[2],paver=aver[3],
                                 max=maxi[1],mmax=maxi[2],pmax=maxi[3],
                                 mu1=mus$mu1[i],mu2=mus$mu2[i],kappa2mu=mus$kappa2mu[i],
                                 mu1col=mu1cols$clr[which(mu1cols$mu1==mus$mu1[i])[1]],
                                 mu2col=mu2cols$clr[which(mu2cols$mu2==mus$mu2[i])[1]],
                                 mu2pch=16,#mu2pch$pch[which(mu2pch$mu2==mus$mu2[i])[1]],
                                 stringsAsFactors=FALSE)
                     )
  }
  options(width=250)
  print(force.df)
  
  fmodel <- list()
  fmodel.pmax <- NULL
  #fmodel[[1]] <- nls(aver~a*(mu2-mu1)^2*((mu2+1e-11)/(mu1+1e-11))^b,data=force.df,start=list(a=0.1,b=0.2),trace=TRUE,algorithm='port',control=list(maxiter=1000))
  if(fit=="pade") {
    fmodel[[1]] <- nls(aver~(a*abs(mu2-mu1)^2+b*abs(mu2-mu1)^4)/(1+c*abs(mu2-mu1)^2),
                       weights=1/force.df$paver^2,
                       algorithm='port',lower=c(0,0,0),
                       data=force.df,start=list(a=10,b=5,c=2),trace=TRUE,control=list(minFactor=0.00000001,maxiter=10000)) 
                       #trace=TRUE,algorithm='port',control=list(maxiter=10000))  
    fmodel[[2]] <- nls(max~(a*abs(mu2-mu1)^2+b*abs(mu2-mu1)^4)/(1+c*abs(mu2-mu1)^2),
                       #weights=1/force.df$pmax^2,
                       data=force.df,start=list(a=10,b=1,c=0.2),trace=TRUE,control=list(minFactor=0.00000001,maxiter=10000))
    #                   trace=TRUE,algorithm='port',control=list(maxiter=10000))  
    fmodel.pmax <- nls(pmax~(a*abs(mu2-mu1)^2+b*abs(mu2-mu1)^4)/(1+c*abs(mu2-mu1)^2),
                       data=force.df,start=list(a=10,b=1,c=0.2),trace=TRUE,control=list(minFactor=0.00000001,maxiter=10000))
    #                   trace=TRUE,algorithm='port',control=list(maxiter=10000))  
  } else if(fit=="paderat") {
    fmodel[[1]] <- nls(aver~(a*(1-abs(mu2/mu1))^2+b*(1-abs(mu2/mu1))^4)/(1+c*(1-abs(mu2/mu1))^2),
                       weights=1/force.df$paver^2,
                       algorithm='port',lower=c(0,0,0),
                       data=force.df,start=list(a=10,b=5,c=2),trace=TRUE,control=list(minFactor=0.00000001,maxiter=10000)) 
                       #trace=TRUE,algorithm='port',control=list(maxiter=10000))  
    fmodel[[2]] <- nls(max~(a*(1-abs(mu2/mu1))^2+b*(1-abs(mu2/mu1))^4)/(1+c*(1-abs(mu2/mu1))^2),
                       #weights=1/force.df$pmax^2,
                       data=force.df,start=list(a=10,b=1,c=0.2),trace=TRUE,control=list(minFactor=0.00000001,maxiter=10000))
    #                   trace=TRUE,algorithm='port',control=list(maxiter=10000))  
    fmodel.pmax <- nls(aver~(a*(1-abs(mu2/mu1))^2+b*(1-abs(mu2/mu1))^4)/(1+c*(1-abs(mu2/mu1))^2),
                       data=force.df,start=list(a=10,b=1,c=0.2),trace=TRUE,control=list(minFactor=0.00000001,maxiter=10000))
    #                   trace=TRUE,algorithm='port',control=list(maxiter=10000))  

  } else if(fit=="pheno") {
    fmodel[[1]] <- nls(aver~a*abs(mu2-mu1)^2*abs(mu2/mu1)^b,
                       ,data=force.df,start=list(a=0.1,b=0.4),
                       trace=TRUE,control=list(minFactor=0.0000001),
                       weights=1/(force.df$paver^2+force.df$maver^2))
    fmodel[[2]] <- nls(max~a*abs(mu2-mu1)^2/abs((mu2)*(mu1))^b,data=force.df,start=list(a=1.0,b=2.0),
                       trace=TRUE,control=list(maxiter=1000),
                       weights=1/(force.df$pmax^2+force.df$mmax^2))
    fmodel.pmax <- nls(pmax~a*abs(mu2-mu1)^2/abs((mu2+1e-9)*(mu1+1e-9))^b,data=force.df,start=list(a=0.1,b=2.0),#,c=0.1),
                       trace=TRUE,algorithm='port',control=list(maxiter=1000))
  } else if(fit=="phenodim") {
    fmodel[[1]] <- nls(aver~a*abs(mu2-mu1)^2*abs(mu1*mu2/kappa2mu^2)^b,
                       data=force.df,start=list(a=10,b=0.7),
                       trace=TRUE,control=list(minFactor=0.0000001,maxiter=10000))#,#)#,
                       #weights=1/(force.df$paver^2+force.df$maver^2))
    fmodel[[2]] <- nls(max~a*abs(mu2-mu1)^2*abs(mu1*mu2/kappa2mu^2)^b,#abs((mu2)*(mu1)/kappa2mu^2)^b,
                       data=force.df,start=list(a=10,b=0.7),
                       trace=TRUE,control=list(maxiter=1000,minFactor=0.000001),
                       algorithm='port',
                       weights=1/(force.df$pmax^2+force.df$mmax^2)
                       )
    fmodel.pmax <- nls(pmax~a*abs(mu2-mu1)^2*abs(mu1*mu2/kappa2mu^2)^b,
                       #*abs(mu2/mu1)^b,
                       data=force.df,start=list(a=20,b=0.7),#,c=0.1),
                       trace=TRUE,algorithm='port',control=list(maxiter=1000))
  } else if(fit=="phenodim2") {
    fmodel[[1]] <- nls(aver~((mu2+a)/(mu1+b))*(mu2-mu1)^2,
                       data=force.df,start=list(a=1,b=0.2),
                       trace=TRUE,control=list(minFactor=0.0000001,maxiter=10000)#,
                       #algorithm='port'#,
                       )
                       #weights=1/(force.df$paver^2+force.df$maver^2))
    fmodel[[2]] <- nls(max~((mu2+a)/(mu1+b))*(mu2-mu1)^2,
                       data=force.df,start=list(a=20,b=0.2),
                       trace=TRUE,control=list(maxiter=1000,minFactor=0.000001),
                       #algorithm='port',
                       weights=1/(force.df$pmax^2+force.df$mmax^2)
                       )
    fmodel.pmax <- nls(pmax~((mu2+a)/(mu1+b))*(mu2-mu1)^2, 
                       #*abs(mu2/mu1)^b,
                       data=force.df,start=list(a=20,b=0.2),
                       trace=TRUE,algorithm='port',control=list(maxiter=1000))
  } else if(fit=="rat") {
    fmodel[[1]] <- nls(aver~a*(1-abs(mu2/mu1))^b,
                       data=force.df,start=list(a=0.1,b=2),
                       algorithm='port',
                       trace=TRUE,control=list(minFactor=0.0000001,maxiter=10000))#,
                       #weights=1/(force.df$paver^2+force.df$maver^2))
    fmodel[[2]] <- nls(max~a*(1-abs(mu2/mu1))^b,data=force.df,start=list(a=1.0,b=2.0),
                       trace=TRUE,control=list(maxiter=10000,minFactor=0.0000001),
                       weights=1/(force.df$pmax^2+force.df$mmax^2))
    fmodel.pmax <- nls(pmax~a*(1-abs(mu2/mu1)^b),data=force.df,start=list(a=0.1,b=2.0),#,c=0.1),
                       trace=TRUE,algorithm='port',control=list(maxiter=10000,minFactor=0.0000001))
  }
   
  # for more involved error analysis 
  #fmodel.tsboot <- lapply(X=df.tsboot,FUN=function(x) {
  #                                                      nls(aver~a*abs(mu2-mu1)^2*abs((mu2+1e-7)/(mu1+1e-7))^b,data=x,start=list(a=0.1,b=1.0),
  #                                                      trace=FALSE,algorithm='port',control=list(maxiter=1000),
  #                                                      weights=1/(force.df$paver^2+force.df$maver^2)) 
  #                                                    }
  #                       )

  #fmodel.coefs.l <- lapply(X=fmodel.tsboot,FUN=function(x) { summary(x)$coefficients[,1] } )
  #fmodel.coefs.df <- do.call("rbind",fmodel.coefs.l)
  #fmodel.coefs <- rbind(apply(X=fmodel.coefs.df,MARGIN=2,FUN=median),apply(X=fmodel.coefs.df,MARGIN=2,FUN=sd))
  
  for(mod in fmodel){
    print(summary(mod))
    #coefs <- summary(mod)$coefficients[,1]
  }
  print(summary(fmodel.pmax))
  save(fmodel,file=sprintf("fmodel.%s.Rdata",basename))
 
  constmu1 <- list()
  constmu2 <- list()
  for(i in 1:nrow(mus)){
    constmu1[[length(constmu1)+1]] <- data.frame(mu1=rep(mus$mu1[i],times=n.predict),mu2=logseq(-4,-0.3,length.out=n.predict),kappa2mu=rep(mus$kappa2mu[i],times=n.predict))
  }
  for(i in 1:nrow(mus)){
    constmu2[[length(constmu2)+1]] <- data.frame(mu1=logseq(-4,-0.3,length.out=n.predict),mu2=rep(mus$mu2[i],times=n.predict),kappa2mu=rep(mus$kappa2mu[i],times=n.predict))
  }
 
  #cA2.30.24 <- data.frame(rho1=c(0.0363,0.0)+0.0008237,rho2=c(0.04,0.008)+0.0008238,clr=rep("blue",2),stringsAsFactors=FALSE)
  #cA2.60.24 <- data.frame(rho1=c(0.011,0.0)+0.0016476,rho2=c(0.06,0.011)+0.0016476,clr=rep("red",2),stringsAsFactors=FALSE)
  #cA2.60.32 <- data.frame(rho1=c(0.08,0.008,0.0)+0.0016476,rho2=c(0.8,0.08,0.008)+0.0016476,clr=rep("magenta",3),stringsAsFactors=FALSE)
  #typ.data <- rbind(cA2.30.24,cA2.60.24,cA2.60.32)
  #typ.y <- predict(rhomodel,newdata=typ.data)

  #rho <- data.frame(rho1=rep(0.0008237+0.014,times=100),rho2=seq(0,0.2,length.out=100))  ,
  #print(cbind(predict(rhomodel,newdata=rho),rho))
  #readline()
  
  tikzfiles <- tikz.init(basename=sprintf("%s.fmodel",basename),width=width,height=height,lwdUnit=0.8)

  # for the plots in terms of mu2, need the margin
  par(mgp=c(3,0.3,0),mar=c(5,4,4,8)+0.1)

  ylims <- c(-5,1)
  xlims <- c(4e-3,0.21)
  plotwitherror(x=force.df$mu2,y=force.df$aver,main="",tck=0.02,xlim=xlims,
                ylab="$\\|F\\|^2_\\mathrm{av}$",xlab="",ylim=10^ylims,las=1,yaxt='n',xaxt='n',log='xy',type='n')
  for(i in 1:length(constmu1)){
    pred.y <- predict(fmodel[[1]],newdata=constmu1[[i]])
    lines(y=pred.y,x=constmu1[[i]]$mu2,col=mu1cols$clr[i],lty=1)
  }
  plotwitherror(x=force.df$mu2,y=force.df$aver,dy=force.df$paver,mdy=force.df$maver,pch=16,col=force.df$mu1col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-2:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(4:9,1e-3),outer(2:9,1e-2),2e-1),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_2$")
  par(xpd=TRUE)
  #legend(x=0.23,y=10^(ylims[1]+1),legend=c("max","aver"),pch=c(17,16),bty='n',lty=c(2,1))
  # to get colours and numbers in the same order
  legend(x=0.25,y=2*10^ylims[2],legend=sprintf("$\\tilde{\\mu}_1=%.7f$",rev(mu1cols$mu1)),fill=rev(mu1cols$clr),bty='n',cex=0.8)
  par(xpd=FALSE)

  ylims <- c(-3,3)
  xlims <- c(4e-3,0.21)
  plotwitherror(x=force.df$mu2,y=force.df$max,main="",tck=0.02,xlim=xlims,
                ylab="$\\|F\\|^2_\\mathrm{max}$",xlab="",ylim=10^ylims,las=1,yaxt='n',xaxt='n',log='xy',type='n')
  for(i in 1:length(constmu1)){
    pred.y <- predict(fmodel[[2]],newdata=constmu1[[i]])
    lines(y=pred.y,x=constmu1[[i]]$mu2,col=mu1cols$clr[i],lty=1)
  }
  plotwitherror(rep=TRUE,x=force.df$mu2,y=force.df$max,dy=force.df$pmax,mdy=force.df$mmax,pch=17,col=force.df$mu1col)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-2:-1),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(4:9,1e-3),outer(2:9,1e-2),2e-1),tck=0.01)
  mtext(side=1,line=1.5,text="$\\tilde{\\mu}_2$")
  par(xpd=TRUE)
  #legend(x=0.23,y=10^(ylims[1]+1),legend=c("max","aver"),pch=c(17,16),bty='n',lty=c(2,1))
  # to get colours and numbers in the same order
  legend(x=0.25,y=2*10^ylims[2],legend=sprintf("$\\tilde{\\mu}_1=%.7f$",rev(mu1cols$mu1)),fill=rev(mu1cols$clr),bty='n',cex=0.8)
  par(xpd=FALSE)

  #####
  # and default margins for the plots in terms of mu1 
  par(mgp=c(3,0.3,0),mar=c(5,4,4,2)+0.1)
  
  # summary plot in terms of |mu2-mu1|
  xlims <- c(1e-1,1e3)
  ylims <- c(-5,0)
  plotwitherror(x=abs(force.df$mu2/force.df$mu1),y=force.df$aver,main="",tck=0.02,
                xlim=xlims,
                ylab="$\\|F\\|_\\mathrm{av}^2$",xlab="",
                ylim=10^ylims,las=1,yaxt='n',xaxt='n',
                log='xy',
                type='n')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel[[1]],newdata=constmu2[[i]])
    lines(y=pred.y,x=abs(constmu2[[i]]$mu2/constmu2[[i]]$mu1),col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(x=abs(force.df$mu2/force.df$mu1),y=force.df$aver,dy=force.df$paver,mdy=force.df$maver,pch=force.df$mu2pch,col=force.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-5:5),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,10^(-2:3))),tck=0.01)
  mtext(side=1,line=1.5,text="$|\\tilde{\\mu}_2/\\tilde{\\mu}_1|$")
  # to get colours and numbers in the same order
  legend("bottomright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),col=mu2cols$clr,pch=mu2pch$pch,bty='n',cex=cex,pt.cex=1.1)
  
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
  legend("topright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),col=mu2cols$clr,pch=mu2pch$pch,bty='n',cex=cex,pt.cex=1.1)

  # summary plot in terms of |mu2-mu1|
  xlims <- c(1e-1,1e3)
  ylims <- c(-3,3)
  plotwitherror(x=abs(force.df$mu2/force.df$mu1),y=force.df$max,main="",tck=0.02,
                xlim=xlims,
                ylab="$\\|F\\|_\\mathrm{max}^2$",xlab="",
                ylim=10^ylims,las=1,yaxt='n',xaxt='n',
                log='xy',
                type='n')
  for(i in 1:length(constmu2)){
    pred.y <- predict(fmodel[[2]],newdata=constmu2[[i]])
    lines(y=pred.y,x=abs(constmu2[[i]]$mu2/constmu2[[i]]$mu1),col=mu2cols$clr[i],lty=1)
  }
  plotwitherror(x=abs(force.df$mu2/force.df$mu1),y=force.df$max,dy=force.df$pmax,mdy=force.df$mmax,pch=force.df$mu2pch,col=force.df$mu2col,rep=TRUE)
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  axis(side=1, labels=TRUE, at=10^(-5:5),tck=0.02,las=1)
  axis(side=1, labels=FALSE, at=c(outer(2:9,10^(-2:3))),tck=0.01)
  mtext(side=1,line=1.5,text="$|\\tilde{\\mu}_2/\\tilde{\\mu}_1|$")
  # to get colours and numbers in the same order
  legend("bottomright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),col=mu2cols$clr,pch=mu2pch$pch,bty='n',cex=cex,pt.cex=1.1)
  
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
  legend("topright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),col=mu2cols$clr,pch=mu2pch$pch,bty='n',cex=cex,pt.cex=1.1)
  
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
  legend("topright",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",rev(mu2cols$mu2)),col=rev(mu2cols$clr),pch=rev(mu2pch$pch),bty='n',cex=cex,pt.cex=1.1)
  
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
  legend("bottomleft",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",rev(mu2cols$mu2)),col=rev(mu2cols$clr),pch=rev(mu2pch$pch),bty='n',cex=cex,pt.cex=1.1)
  
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
  legend("bottomleft",legend=sprintf("$\\tilde{\\mu}_2=%.7f$",mu2cols$mu2),col=mu2cols$clr,bty='n',cex=cex,pt.cex=1.1,pch=mu2pch$pch)
  
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


