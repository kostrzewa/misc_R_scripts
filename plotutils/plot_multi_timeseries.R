source("~/code/R/misc_R_scripts/plotutils/plot_multihist.R")

plot_multi_timeseries <- function(dat,basename,ylims,hl,yticks,width=4.5,height=4.5,name="observable",dbg=FALSE) {

  if(dbg==TRUE) for( d in dat ) { print(summary(d)) }

  # determine plot boundaries from data frame and prepare raw data for multihist
  rdat <- list()
  xlims <- NULL
  lg <- NULL
  for( d in dat ){
    r <- range(d$ts)
    xlims <- c(xlims,d$t)
    
    rdat[[length(rdat)+1]] <- d$ts[d$idx[1]:d$idx[2]]
    if(dbg){
      print(d$name)
      print(d$idx[2]-d$idx[1])
      print(d$t[2]-d$t[1])
      print(c(mean(rdat[[length(rdat)]]),var(rdat[[length(rdat)]])))
    }
    lg <- c(lg,d$name)
  }
  
  if(missing(ylims)) ylims <- range(unlist(rdat))
  xlims <- range(xlims)
  if(dbg==TRUE) print(ylims); print(xlims)

  # set up colours
  clr <- c("red","blue")
  if(length(dat)>2){
    require("RColorBrewer")
    n <- length(dat)
    # if using colourblind palette, better use even number of colours
    # because this avoids really pastel hues 
    #if(n %% 2 != 0) n <- n+1
    clr <- brewer.pal(n=n,name="Dark2")
  }
  

  # set up the plot region
  tikzfiles <- tikz.init(basename=basename,width=width,height=height)
  #op <- par(family="Palatino",cex.main=0.8,font.main=1)
  par(mgp=c(1.0,0.2,0))
  if(!missing(yticks)){
    plot(x=NA,y=NA,t='n',xlim=xlims,ylim=ylims,ylab="",xlab="",las=1,tck=0.02,yaxp=c(ylims,yticks))
  } else {
    plot(x=NA,y=NA,t='n',xlim=xlims,ylim=ylims,ylab="",xlab="",las=1,tck=0.02)
  }
  for( i in 1:length(dat) ){
    lines(y=rdat[[i]], x=seq(from=dat[[i]]$t[1],to=dat[[i]]$t[2],by=dat[[i]]$stepsize),
          col=clr[i], lty=1, lwd=2 )
  }
  text(labels=name,x=xlims[1],y=ylims[2],adj=c(0,1))
  text(labels="$t_\\mathrm{MD}$",x=xlims[2],y=ylims[1],adj=c(1,0))
  # add a horizontal line to the plot, if desired (also add to legend)
  if(!missing(hl)) { 
    abline(h=hl,lty=2); 
    legend(x="topright", legend=c(lg,sprintf("%s=%s",name,hl)), pch=NA, col=c(clr[1:length(dat)],"black"), lty=c(rep(1,times=length(dat)),2), lwd=3, bty='n')
  } else {  
    legend(x="topright", legend=lg, pch=NA, col=clr[1:length(dat)], lty=1, lwd=3, bty='n')
  }  

  plot_multihist(dat=rdat,lg=lg,cols=clr,init=FALSE,label.x=name,main="",lg.pos="topleft",xlim.probs=c(0.001,0.999))

  tikz.finalize(tikzfiles)

}
