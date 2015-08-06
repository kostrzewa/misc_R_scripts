plot_forces <- function(datfile,basename="forces",width=10,height=8,ylims=c(-5,3)) {
  fdat <- read.table(header=TRUE,file=datfile,stringsAsFactors=FALSE)
  mon <- unique(fdat$mon)
  
  nsteps <- length(which(fdat$mon==mon[1]))

  require("RColorBrewer")
  clr <- brewer.pal(name="Paired",n=12)

  tikzfiles <- tikz.init(basename=basename,width=width,height=height)
  # we will draw legends outside the plot area on the right, make some space there
  par(mar=c(5,4,4,12),xpd=TRUE)
  # prepare plot area
  plot(NULL,
       xlab="$t_\\mathrm{MD}$",ylab="$\\bar{F}^2$",main="",las=1,tck=0.02,
       ylim=10^ylims,xlim=c(1,nsteps),log='y',yaxt='n',
      )
  # draw major and minor tick marks
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)

  for( i in 1:length(mon) ) {
    lines(lty=1,col=clr[i],x=1:nsteps,y=fdat$aver[fdat$mon==mon[i]],lwd=3)
    lines(lty=3,col=clr[i],x=1:nsteps,y=fdat$max[fdat$mon==mon[i]],lwd=3)
  }
  legend(lty=1,lwd=3,col=clr[1:length(mon)],legend=mon,x=1.05*nsteps,y=10^ylims[2],bty='n')
  legend(lty=c(3,1),legend=c("max","average"),col="black",lwd=3,x=1.05*nsteps,y=1,bty='n')
  
  # and plots for each monomial seperately
  for( i in 1:length(mon) ) {
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
  for( i in 1:length(mon) ) {
    ind <- which(fdat$mon==mon[i])
    hst.max <- hist(fdat$max[ind],plot=FALSE,breaks=50) 
    plot(hst.max,col=clr[i],main="",xlab="$\\max(F^2)$")
  }
  tikz.finalize(tikzfiles)
}


