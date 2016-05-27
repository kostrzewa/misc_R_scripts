plot_multires <- function(dat,basename,lg,hl,stride=1,width=4.5,height=4.5,shift=TRUE,...){
  start.idx <- which(dat$iters[,1]==1)
  end.idx <- length(dat$iters[,1])
  if(length(start.idx)>1){
    end.idx <- NULL
    for( i in 2:length(start.idx) ){
      end.idx <- c(end.idx,start.idx[i]-1)
    }
    end.idx <- c(end.idx,length(dat$iters[,1]))
  }

  clr <- c("red","blue")
  if(length(start.idx)>2){
    require("RColorBrewer")
    clr <- brewer.pal(n=length(start.idx),name="Dark2")
  }

  tikzfiles <- tikz.init(basename=basename,width=width,height=height)
  par(mgp=c(3,0.3,0))
  xlims <- c(1,max(end.idx))
  if(shift==FALSE) xlims <- c(1,max(dat$iters[,1]))
  plot(x=NA,y=NA,ylim=c(min(dat$iters[,2]),max(dat$iters[,2])),xlim=xlims,
       tck=0.02,yaxt='n',
       log='y',ylab="$\\|r\\|^2$",xlab="",...)
  
  label <- "$N_\\mathrm{iter}$"
  if(shift==TRUE) label <- sprintf("%s %s",label," (shifted)")
  text(labels=label,x=0.8*xlims[2],y=min(dat$iters[,2]),adj=c(0,0))
  
  
  for( i in 1:length(start.idx) ){
    x <- dat$iters[ seq(from=start.idx[i],to=end.idx[i],by=stride), 1 ]
    if(shift==TRUE) x <- x + (start.idx[i]-1)
    y <- dat$iters[ seq( from=start.idx[i],to=end.idx[i],by=stride), 2 ]
    lines(x=x,y=y,col=clr[i],lty=i,lwd=3)
    points(x=max(x),y=dat$real.res[i],col=clr[i],pch=4,cex=2)
  }
  axis(side=2, labels=TRUE, at=10^(-30:30),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(-30:(30-1))),tck=0.01)
  if(!missing(hl)) { abline(h=hl,lty=2); text(x=0,y=hl,labels=sprintf("$\\|r\\|^2=$%s",hl),adj=c(0,1.7)) }

  if(!missing(lg)){
    legend(x="topright",legend=lg,lty=1:length(start.idx),col=clr[1:length(start.idx)],lwd=3,bty='n')
  }

  tikz.finalize(tikzfiles)
}
