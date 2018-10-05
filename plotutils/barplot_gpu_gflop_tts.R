barplot_gpu_gflops_tts <- function(datfile,basename="gpuperf",dbg=FALSE) {
  require("RColorBrewer")
  
  dat <- read.table(datfile,stringsAsFactors=FALSE,header=TRUE)

  precs <- unique(dat$prec)
  devs <- unique(dat$devs)
  gflops <- matrix(nrow=length(precs),ncol=length(devs))

  for(i in 1:length(devs)){
    gflops[,i] <- dat[dat$devs==devs[i],]$gflops
  }
  
  pal <- brewer.pal(n=length(precs),name="Blues")
  tikzFiles <- tikz.init(basename=basename,width=4.0,height=4.0)
  par(mgp=c(2.5,0.3,0.5))
  par(mar=c(5,4,4,5)+.1)
  mids <- barplot(gflops,beside=TRUE,las=1,ylab="",xlab="",
                  col=pal,
                  ylim=c(0,3000),xaxs='i',yaxs='i',
                  tck=0.015)
  mtext(side=2,line=2.5,"Gflop/s")
  mtext(side=1,line=0.5,(devs/4),at=mids[c(2,5,8)])
  mtext(side=1,line=1.5,text="Jureca K80 GPU Nodes")
  legend('topleft',cex=0.77,legend=c("CG double (Gflop/s)",sprintf("MixCG %s (Gflop/s)",precs[2:3]),"speed-up"),
         pch=c(21,22,23,NA),pt.cex=1.5,col=c(pal,'black'),lty=c(rep(NA,3),1),lw=3,bty='n',pt.bg=c(pal,NA))
  xlims <- par("usr")[1:2]
  
  # add points and lines for speedup factors
  par(new=TRUE)
  speedup <- NULL
  for( prec in precs ){
    idx <- which(dat$prec==prec)
    speedup <- c(speedup,max(dat[,'tts'])/dat[idx,'tts'])
  }
  dat <- cbind(dat,speedup)
  if(dbg) print(dat)

  plot(x=NA,y=NA,t='n',
       lwd=3,xlab='',ylab='',xaxt='n',yaxt='n',
       #xlim=0.925*(xlims+0.9*xlims[1]),ylim=c(0,range(dat$speedup)[2]*1.3),
       xlim=xlims,ylim=c(0,1.3*max(dat$speedup)),xaxs='i',yaxs='i',
       #1.3*range(dat$tts)-50,
       bty='n')
  abline(h=1,lty=3)
  for(i in 1:3){
    #lines(x=mids[c(i,i+3,i+6)],y=dat[dat$prec==precs[i],]$speedup,lwd=3)
    #points(x=mids[c(i,i+3,i+6)],y=dat[dat$prec==precs[i],]$speedup,lwd=3,pch=20+i,bg=pal[i],cex=1.2)
    lines(x=mids[c(2,5,8)],y=dat[dat$prec==precs[i],]$speedup,lwd=3)
    points(x=mids[c(2,5,8)],y=dat[dat$prec==precs[i],]$speedup,lwd=3,pch=20+i,bg=pal[i],cex=1.2)
  }
  axis(4,las=1,at=0:9,tck=0.015)
  mtext(side=4,"speed-up",line=1.3)
  tikz.finalize(tikzFiles)
}
  
