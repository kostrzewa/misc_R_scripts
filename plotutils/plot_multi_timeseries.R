source("~/code/R/misc_R_scripts/plotutils/plot_multihist.R")

plot_multi_timeseries <- function(dat,basename,ylims,hl,width=4.5,height=4.5,name="observable",dbg=FALSE) {

  if(missing(ylims)) ylims <- c(0,0)
  xlims <- c(0,0)
  
  # determine plot boundaries from data frame and prepare raw data for multihist
  rdat <- list()
  lg <- NULL
  for( d in dat ){
    r <- range(d$ts)
    if(missing(ylims)){
      if( r[1] < ylims[1] ) ylims[1] <- r[1]
      if( r[2] > ylims[2] ) ylims[2] <- r[2]
    }
    if( d$t[1] < xlims[1] ) xlims[1] <- d$t[1]
    if( d$t[2] > xlims[2] ) xlims[2] <- d$t[2]
    rdat[[length(rdat)+1]] <- d$ts[d$idx[1]:d$idx[2]]
    if(dbg){
      print(d$name)
      print(d$idx[2]-d$idx[1])
      print(d$t[2]-d$t[1])
    }
    lg <- c(lg,d$name)
  }

  # set up colours
  clr <- c("red","blue")
  if(length(dat)>2){
    require("RColorBrewer")
    clr <- brewer.pal(n=length(dat),name="Set2")
  }
  

  # set up the plot region
  tikzfiles <- tikz.init(basename=basename,width=width,height=height)
  #op <- par(family="Palatino",cex.main=0.8,font.main=1)
  par(mgp=c(1.0,0.2,0))
  plot(x=NA,y=NA,t='n',xlim=xlims,ylim=ylims,ylab="",xlab="",las=1,tck=0.02,yaxp=c(ylims,10))
  for( i in 1:length(dat) ){
    lines(y=rdat[[i]], x=seq(from=dat[[i]]$t[1],to=dat[[i]]$t[2],by=dat[[i]]$stepsize),
          col=clr[i], lty=i, lwd=2 )
  }
  text(labels=name,x=xlims[1],y=ylims[2])
  text(labels="$t_\\mathrm{MD}$",x=xlims[2],y=ylims[1])
  # add a horizontal line to the plot, if desired (also add to legend)
  if(!missing(hl)) { 
    abline(h=hl,lty=2); 
    legend(x="topright", legend=c(lg,sprintf("%s=%s",name,hl)), pch=NA, col=c(clr[1:length(dat)],"black"), lty=c(1:length(dat),2), lwd=2, bty='n')
  } else {  
    legend(x="topright", legend=lg, pch=NA, col=clr[1:length(dat)], lty=1:length(dat), lwd=2, bty='n')
  }

  plot_multihist(dat=rdat,lg=lg,cols=clr,init=FALSE,label.x=name,main="",lg.pos="topleft",xlim.probs=c(0.01,1.0))

  tikz.finalize(tikzfiles)

}
