source("~/code/R/misc_R_scripts/plotutils/plot_multihist.R")
require("parallel")

plot_forces <- function(datfile,nsteps,lg,basename="forces",width=10,height=8,ylims=c(-5,3),summary=TRUE,analysis=FALSE,n.predict=500) {
  fdat <- read.table(header=TRUE,file=datfile,stringsAsFactors=FALSE)
  mon <- unique(fdat$mon)
 
  if(missing(lg)){
    lg <- mon
  }
  
  if(missing(nsteps)) {
    nsteps <- length(which(fdat$mon==mon[1]))
  } else {
    fdat <- fdat[1:(length(mon)*nsteps),]
  }
  
  require("RColorBrewer")
  clr <- rep(brewer.pal(name="Set1",n=9),2)
  
  if(summary){
    tikzfiles <- tikz.init(basename=basename,width=width,height=height)
    # we will draw legends outside the plot area on the right, make some space there
    par(mar=c(5,4,4,12),xpd=TRUE,mgp=c(3,0.3,0))
    # prepare plot area
    plot(NULL,
         xlab="",ylab="$\\|F\\|^2$",main="",las=1,tck=0.02,
         ylim=10^ylims,xlim=c(1,nsteps),log='y',yaxt='n',
        )
    # draw major and minor tick marks
    axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
    axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
    mtext(side=1, line=1.1, text="$t_\\mathrm{MD}$") 
     
    for( i in 1:length(mon) ) {
      # in case there is some slight disagreement among the number of steps, we determine them for each monomial separately
      nsteps <- length(fdat$aver[fdat$mon==mon[i]])
      lines(lty=1,col=clr[i],x=1:nsteps,y=fdat$aver[fdat$mon==mon[i]],lwd=3)
      lines(lty=3,col=clr[i],x=1:nsteps,y=fdat$max[fdat$mon==mon[i]],lwd=3)
    }
    legend(lty=1,lwd=3,col=clr[1:length(mon)],legend=lg,x=1.05*nsteps,y=10^ylims[2],bty='n')
    legend(lty=c(3,1),legend=c("max","average"),col="black",lwd=3,x=1.05*nsteps,y=1,bty='n')
    
    # and plots for each monomial seperately
    for( i in 1:length(mon) ) {
      nsteps <- length(fdat$aver[fdat$mon==mon[i]])
      ind <- which(fdat$mon==mon[i])
      yrange <- range(c(fdat$aver[ind],fdat$max[ind]))
      plot(NULL,
           xlab="$t_\\mathrm{MD}$",ylab="$\\bar{F}^2$",main="",las=1,tck=0.02,
           ylim=yrange,xlim=c(1,nsteps)
          )
      lines(lty=1,col=clr[i],x=1:nsteps,y=fdat$aver[ind],lwd=3)
      lines(lty=3,col=clr[i],x=1:nsteps,y=fdat$max[ind],lwd=3)
      
      legend(lty=1,lwd=3,col=clr[1:length(mon)],legend=mon,x=1.05*nsteps,y=yrange[2],bty='n')
      legend(lty=c(3,1),legend=c("max","average"),col="black",lwd=3,x=1.05*nsteps,y=(yrange[2]-yrange[1])/2,bty='n')
    }
    tikz.finalize(tikzfiles)
      
    tikzfiles <- tikz.init(basename=sprintf("%s.hist",basename),width=3.3,height=3.3)
    hst <- list()
    for( i in 1:length(mon) ) {
      ind <- which(fdat$mon==mon[i])
      hst[[length(hst)+1]] <- fdat$max[ind]
      # remove outliers from individual histogram plots so that the plots are clear and there is a sufficient number of breaks
      aver.qt <- quantile(fdat$aver[ind],probs=c(0.01,0.99))
      max.qt <- quantile(fdat$max[ind],probs=c(0.01,0.99))
      
      dat.aver <- c(fdat$aver[ind])[ which( fdat$aver[ind] >= aver.qt[1] & fdat$aver[ind] <= aver.qt[2] ) ]
      dat.max <- c(fdat$max[ind])[ which( fdat$max[ind] >= max.qt[1] & fdat$max[ind] <= max.qt[2] ) ]
      hst.aver <- hist( dat.aver, #fdat$aver[ind],
                       plot=FALSE, breaks=50)
      hst.max <- hist( dat.max, #fdat$max[ind],
                      plot=FALSE, breaks=50)
      plot(hst.aver,col=clr[i],main="",xlab="$\\mathrm{avg}(F^2)$")
      plot(hst.max,col=clr[i],main="",xlab="$\\max(F^2)$")
    }
    tikz.finalize(tikzfiles)
    plot_multihist(dat=hst,basename=sprintf("%s.multihist",basename),factor=50)
  } 
  if(analysis) {
    av <- NULL
    for( i in 1:length(mon) ) {
      if(any(mon[i]==c("GAUGE","rho10rho10","cloverdet","cloverdetratio1","cloverdetratio2","cloverdetratio3","cloverdetratio4","cloverdetratio5"))) next
      ind <- which(fdat$mon==mon[i])
      nsteps <- length(fdat$aver[ind])
      aver.tsboot <- tsboot(tseries=fdat$aver[ind],R=3*nsteps,l=2,sim="geom",statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) })
      max.tsboot <- tsboot(tseries=fdat$max[ind],R=3*nsteps,l=2,sim="geom",statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) })
      aver <- c(mean(aver.tsboot$t[,1]),mean(aver.tsboot$t[,1]-aver.tsboot$t[,2]),mean(aver.tsboot$t[,3]-aver.tsboot$t[,1]))
      maxi  <- c(mean(max.tsboot$t[,1]),mean(max.tsboot$t[,1]-max.tsboot$t[,2]),mean(max.tsboot$t[,3]-max.tsboot$t[,1]))
      av <- rbind(av,data.frame(name=mon[i],aver=aver[1],maver=aver[2],paver=aver[3],max=maxi[1],mmax=maxi[2],pmax=maxi[3],stringsAsFactors=FALSE))
      #cat(sprintf("%20s %.4g -%.4g +%.4g \t %.4g -%.4g +%.4g\n",mon[i],aver[1],aver[2],aver[3],maxi[1],maxi[2],maxi[3]))
    }
    options(width=160)
 
    av2 <- NULL
    for( i in 1:length(mon) ) {
      if(!any(mon[i]==c("GAUGE","cloverdet",sprintf("cloverdetratio%d",1:5)))) next
      cat(sprintf("Bootstrapping %s\n",mon[i]))
      ind <- which(fdat$mon==mon[i])
      nsteps <- length(fdat$aver[ind])
      R <- 3*nsteps
      if(nsteps > 1000) R <- nsteps
      aver.tsboot <- tsboot(tseries=fdat$aver[ind],R=R,l=2,sim="geom",
                            statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) },
                            parallel="multicore",ncpus=8)
      max.tsboot <- tsboot(tseries=fdat$max[ind],R=R,l=2,sim="geom",
                           statistic=function(x){ quantile(x,probs=c(0.5,0.1573,0.8427)) },
                           parallel="multicore",ncpus=8)
      aver <- c(mean(aver.tsboot$t[,1]),mean(aver.tsboot$t[,1]-aver.tsboot$t[,2]),mean(aver.tsboot$t[,3]-aver.tsboot$t[,1]))
      maxi  <- c(mean(max.tsboot$t[,1]),mean(max.tsboot$t[,1]-max.tsboot$t[,2]),mean(max.tsboot$t[,3]-max.tsboot$t[,1]))
      print(aver)
      print(maxi)
      av2 <- rbind(av2,data.frame(name=mon[i],aver=aver[1],maver=aver[2],paver=aver[3],max=maxi[1],mmax=maxi[2],pmax=maxi[3],stringsAsFactors=FALSE))
      #cat(sprintf("%20s %.4g -%.4g +%.4g \t %.4g -%.4g +%.4g\n",mon[i],aver[1],aver[2],aver[3],maxi[1],maxi[2],maxi[3]))
    }

    print(av2)
    
    xy <- NULL
    RHO1 <- c(0.001,0.002,0.004,0.008,0.016)+0.0002746
    RHO2 <- c(0.4,0.2,0.1,0.08,0.04,0.02)+0.0002746
    rho1cols <- rev(brewer.pal(name="Dark2",n=length(RHO1)))
    rho2cols <- brewer.pal(name="Dark2",n=length(RHO2))
    
    for(i in 1:length(RHO1)){
      rho1 <- RHO1[i]
      for(j in 1:length(RHO2)){
        rho2 <- RHO2[j]
        xy <- rbind(xy,data.frame(rho1=rho1,rho2=rho2,rho1col=rho1cols[i],rho2col=rho2cols[j],stringsAsFactors=FALSE))
      }
    }
    
    rho1nd <- list()
    rho2nd <- list()
    for(rho1 in RHO1){
      rho1nd[[length(rho1nd)+1]] <- data.frame(rho1=rep(rho1,times=n.predict),rho2=seq(0,1,length.out=n.predict))
    }
    for(rho2 in RHO2){
      rho2nd[[length(rho2nd)+1]] <- data.frame(rho1=seq(0,0.2,length.out=n.predict),rho2=rep(rho2,times=n.predict))
    }

    av <- cbind(xy,av)
    options(width=250)
    print(av)

    #rhomodel <- nls(aver~a*rho2^b*exp(c*rho1/sqrt(rho2)),data=av,start=list(a=0.1,b=0.1,c=0.1),trace=TRUE,algorithm='port',control=list(maxiter=1000))
    rhomodel <- nls(aver~a*(rho2-rho1)^2*(rho1/rho2)^b,data=av,start=list(a=0.1,b=0.2),trace=TRUE,algorithm='port',control=list(maxiter=1000))
    print(summary(rhomodel))
    coefs <- summary(rhomodel)$coefficients[,1]

    #typ.data <- data.frame(rho1=c(0.01,0.001,0.0),rho2=c(0.1,0.01,0.001))
    u <- 3.1
    #u.rho <- data.frame(rho1=u^
    #print(cbind(u.rho,u.y))
#    cA2.09.48 <- data.frame(rho1=c(0.011,0.0)+0.0016476,rho2=c(0.06,
    cA2.30.24 <- data.frame(rho1=c(0.0363,0.0)+0.0008237,rho2=c(0.04,0.008)+0.0008238,clr=rep("blue",2),stringsAsFactors=FALSE)
    cA2.60.24 <- data.frame(rho1=c(0.011,0.0)+0.0016476,rho2=c(0.06,0.011)+0.0016476,clr=rep("red",2),stringsAsFactors=FALSE)
    cA2.60.32 <- data.frame(rho1=c(0.08,0.008,0.0)+0.0016476,rho2=c(0.8,0.08,0.008)+0.0016476,clr=rep("magenta",3),stringsAsFactors=FALSE)
    typ.data <- rbind(cA2.30.24,cA2.60.24,cA2.60.32)
    typ.y <- predict(rhomodel,newdata=typ.data)

    #rho <- data.frame(rho1=rep(0.0008237+0.014,times=100),rho2=seq(0,0.2,length.out=100))  
    #print(cbind(predict(rhomodel,newdata=rho),rho))
    #readline()

    tikzfiles <- tikz.init(basename=sprintf("%s.rho",basename),width=4,height=8,lwdUnit=0.8)
    par(mgp=c(3,0.3,0))
    ylims <- c(-5,3)
    plotwitherror(x=xy$rho2,y=av$av,main="",tck=0.02,xlim=c(0,max(xy$rho2)),
                  ylab="$F^2$",xlab="$\\rho_2$",ylim=10^ylims,las=1,yaxt='n',log='y',type='n')
    for(i in 1:length(rho1nd)){
      pred.y <- predict(rhomodel,newdata=rho1nd[[i]])
      lines(y=pred.y,x=rho1nd[[i]]$rho2,col=rho1cols[i])
    }
    plotwitherror(x=xy$rho2,y=av$aver,dy=av$paver,mdy=av$maver,pch=16,col=xy$rho1col,rep=TRUE)
    plotwitherror(rep=TRUE,x=xy$rho2,y=av$max,dy=av$pmax,mdy=av$mmax,pch=17,col=xy$rho1col)
    plotwitherror(x=c(0.3,0.03)+0.002746,y=av2$aver,dy=av2$paver,mdy=av2$maver,pch=16,cex=1.2,rep=TRUE)
    plotwitherror(x=c(0.3,0.03)+0.002746,y=av2$max,dy=av2$pmax,mdy=av2$mmax,pch=17,cex=1.2,rep=TRUE)
    points(x=typ.data$rho2,y=typ.y,pch=5,cex=2,col=typ.data$clr)
    axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
    axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
    legend("topright",legend=c("max","aver"),pch=c(17,16),pt.cex=1.3,bty='n')
    # to get colours and numbers in the same order
    legend("bottomright",legend=sprintf("$\\rho_1=%.3f$",rev(RHO1)),fill=rev(rho1cols),bty='n')

    #####

    plotwitherror(x=xy$rho1,y=av$aver,main="",tck=0.02,xlim=c(0,0.15),
                  ylab="$F^2$",xlab="$\\rho_1$",ylim=10^ylims,las=1,yaxt='n',log='y',type='n')
    for(i in 1:length(rho2nd)){
      pred.y <- predict(rhomodel,newdata=rho2nd[[i]])
      lines(y=pred.y,x=rho2nd[[i]]$rho1,col=rho2cols[i])
    }
    
    plotwitherror(x=xy$rho1,y=av$aver,dy=av$paver,mdy=av$maver,pch=16,col=xy$rho2col,rep=TRUE)
    plotwitherror(rep=TRUE,x=xy$rho1,y=av$max,dy=av$pmax,mdy=av$mmax,pch=17,col=xy$rho2col)
    plotwitherror(x=c(0.03,0.0)+0.002746,y=av2$aver,dy=av2$paver,mdy=av2$maver,pch=16,cex=1.2,rep=TRUE)
    plotwitherror(x=c(0.03,0.0)+0.002746,y=av2$max,dy=av2$pmax,mdy=av2$mmax,pch=17,cex=1.2,rep=TRUE)
    points(x=typ.data$rho1,y=typ.y,pch=5,cex=2,col=typ.data$clr)
    axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
    axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
    legend("topright",legend=c("max","aver"),pch=c(17,16),pt.cex=1.3,bty='n')
    legend("bottomleft",legend=sprintf("$\\rho_2=%.3f$",RHO2),fill=rho2cols,bty='n')

    plotwitherror(x=xy$rho2,y=av$max-av$aver,dy=sqrt(av$paver^2+av$pmax^2),mdy=sqrt(av$maver^2+av$mmax^2),
                  xlab="$\\rho_2$",ylab="$\\max(F^2)-\\bar{F}^2$",tck=0.02,log='y',col=xy$rho1col,pch=16,las=1)
    # to get colours and numbers in the same order
    legend("bottomright",legend=sprintf("$\\rho_1=%.3f$",rev(RHO1)),fill=rev(rho1cols),bty='n')
    
    plotwitherror(x=xy$rho1,y=av$max-av$aver,dy=sqrt(av$paver^2+av$pmax^2),mdy=sqrt(av$maver^2+av$mmax^2),
                  xlab="$\\rho_1$",ylab="$\\max(F^2)-\\bar{F}^2$",tck=0.02,log='y',col=xy$rho2col,pch=16,las=1)
    legend("bottomleft",legend=sprintf("$\\rho_2=%.3f$",RHO2),fill=rho2cols,bty='n')


    cA2.30.24 <- data.frame(rho1=c(0.008,0.0)+0.0008237,rho2=c(0.04,0.008)+0.0008238,clr=rep("blue",2),stringsAsFactors=FALSE)
    cA2.60.24 <- data.frame(rho1=c(0.011,0.0)+0.0016476,rho2=c(0.06,0.011)+0.0016476,clr=rep("red",2),stringsAsFactors=FALSE)
    cA2.60.32 <- data.frame(rho1=c(0.08,0.008,0.0)+0.0016476,rho2=c(0.8,0.08,0.008)+0.0016476,clr=rep("magenta",3),stringsAsFactors=FALSE)
    typ.data <- rbind(cA2.30.24,cA2.60.24,cA2.60.32)
    typ.y <- predict(rhomodel,newdata=typ.data)

    rho1 <- sort(unique(typ.data$rho1))
    rho2 <- sort(unique(typ.data$rho2))

    rho1cols <- brewer.pal(name="Dark2",n=length(rho1))
    rho2cols <- brewer.pal(name="Dark2",n=length(rho2))

    plotwitherror(x=xy$rho2,y=av$aver,main="",tck=0.02,xlim=c(0,0.4),
                  ylab="$F^2$",xlab="$\\rho_2$",ylim=10^ylims,las=1,yaxt='n',log='y',type='n')
    for(i in 1:length(rho1)){
      newdata <- data.frame(rho1=rep(rho1[i],times=n.predict),rho2=seq(0,0.5,length.out=n.predict))
      pred.y <- predict(rhomodel,newdata=newdata)
      lines(y=pred.y,x=newdata$rho2,col=rho1cols[i])
      #print(data.frame(pred.y[1:100],newdata$rho2[1:100],newdata$rho2[1:100]-0.0008238))
      #readline()
    }
    points(x=typ.data$rho2,y=typ.y,pch=5,cex=2,col=typ.data$clr)
    axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
    axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
    legend("topright",legend=c("max","aver"),pch=c(17,16),pt.cex=1.3,bty='n')
    legend("bottomright",legend=sprintf("$\\rho_1=%.7f$",rho1),fill=rho1cols,bty='n')
    
    plotwitherror(x=xy$rho1,y=av$aver,main="",tck=0.02,xlim=c(0,0.1),
                  ylab="$F^2$",xlab="$\\rho_1$",ylim=10^ylims,las=1,yaxt='n',log='y',type='n')
    for(i in 1:length(rho2)){
      newdata <- data.frame(rho2=rep(rho2[i],times=n.predict),rho1=seq(0,0.5,length.out=n.predict))
      pred.y <- predict(rhomodel,newdata=newdata)
      lines(y=pred.y,x=newdata$rho1,col=rho2cols[i])
    }
    points(x=typ.data$rho1,y=typ.y,pch=5,cex=2,col=typ.data$clr)
    axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
    axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
    legend("bottomright",legend=c("max","aver"),pch=c(17,16),pt.cex=1.3,bty='n')
    legend("topright",legend=rev(sprintf("$\\rho_2=%.7f$",rho2)),fill=rev(rho2cols),bty='n')

    tikz.finalize(tikzfiles)
  } # if(analysis)
}


