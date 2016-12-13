# function to plot timeseries data, a corresponding histogram
# and an error shading for an error analysis via uwerr

plot_timeseries <- function(dat,trange,pdf.filename,
                            ylab,plotsize,titletext,hist.by,
                            name="",xlab="$t_\\mathrm{MD}$",hist.probs=c(0.0,1.0),errorband_color=rgb(0.6,0.0,0.0,0.6),stepsize=1,
                            uwerr.S=2,periodogram=FALSE,debug=FALSE,uw.summary=TRUE,...) {
  xdat <- seq(1,length(dat),stepsize)
  if(!missing(trange)) xdat <- seq(trange[1],trange[2],stepsize)

  yrange <- range(dat)
  
  uw.data <- uwerrprimary(dat,S=uwerr.S)
  if(debug) {
    print(paste("uw.",name,sep=""))
    print(summary(uw.data))
  }
  
  tikzfiles <- tikz.init(basename=pdf.filename,width=plotsize,height=plotsize)
  op <- par(family="Palatino",cex.main=0.8,font.main=1)
  par(mgp=c(2,1.0,0))

  # plot the timeseries
  plot(x=xdat,xlim=range(xdat),y=dat,ylab=ylab,t='l',xlab=xlab,main=titletext,...)

  rect(xleft=range(xdat)[1],
       xright=range(xdat)[2],
       ytop=uw.data$value+uw.data$dvalue,
       ybottom=uw.data$value-uw.data$dvalue,border=FALSE,col=errorband_color)
  abline(h=uw.data$value,col="black")                                                                                                   
  # plot the corresponding histogram
  hist.data <- NULL
  
  hist.breaks <- floor( ( max(dat)-min(dat) ) / uw.data$dvalue )
  if(!missing(hist.by)){
    hist.breaks <- floor( ( max(dat)-min(dat) ) / hist.by )
  } else {
    if(hist.breaks < 10 || hist.breaks > 150){
      hist.breaks <- 70
    }
  }
  
  hist.data <- hist(dat,xlim=quantile(dat,probs=hist.probs),main=titletext,xlab=ylab, breaks=hist.breaks)
  rect(ytop=max(hist.data$counts),
       ybottom=0,
       xright=uw.data$value+uw.data$dvalue,
       xleft=uw.data$value-uw.data$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.data$value,col="black")                                                                                             

  # and a periodogram
  if(periodogram)
  {
    spec.pgram(x=dat,main=paste(ylab,paste("raw periodogram",titletext)))
  }

  # and the uwerr plots
  if(uw.summary){
    plot(uw.data,main=paste(ylab,paste("UWErr analysis",titletext)),x11=FALSE,plot.hist=FALSE)
  }
   
  tikz.finalize(tikzfiles)
  #dev.off()

  return(t(data.frame(val=uw.data$value, dval=uw.data$dvalue, tauint=uw.data$tauint, dtauint=uw.data$dtauint, Wopt=uw.data$Wopt)))
}
