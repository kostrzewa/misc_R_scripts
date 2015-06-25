# function to plot timeseries data, a corresponding histogram
# and an error shading for an error analysis via uwerr

plot_timeseries <- function(dat,trange,pdf.filename,
                            ylab,name,plotsize,filelabel,titletext,errorband_color=rgb(0.6,0.0,0.0,0.6),stepsize=1,
                            hist.breaks=30,uwerr.S=5,periodogram=FALSE,debug=FALSE,...) {
  xdat <- seq(trange[1],trange[2],stepsize)
  yrange <- range(dat)
      
  uw.data <- uwerrprimary(dat,S=uwerr.S)
  if(debug) {
    print(paste("uw.",name,sep=""))
    summary(uw.data)
  }
  
  tikzfiles <- tikz.init(basename=pdf.filename,width=plotsize,height=plotsize)
  #pdf(pdf.filename,width=plotsize,height=plotsize,title=filelabel)
  op <- par(family="Palatino",cex.main=0.6,font.main=1)
  par(mgp=c(2,1,0))

  # plot the timeseries
  plot(x=xdat,xlim=trange,y=dat,ylim=yrange,ylab=ylab,t='l',xlab=expression(t[MD]),main=titletext,...)

  rect(xleft=trange[1],
       xright=trange[2],
       ytop=uw.data$value+uw.data$dvalue,
       ybottom=uw.data$value-uw.data$dvalue,border=FALSE,col=errorband_color)
  abline(h=uw.data$value,col="black")                                                                                                   
  # plot the corresponding histogram
  hist.data <- hist(dat,xlim=yrange,main=paste("histogram",titletext),xlab=ylab, breaks=hist.breaks)
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
  plot(uw.data,main=paste(ylab,paste("UWErr analysis",titletext)),x11=FALSE,plot.hist=FALSE)
   
  tikz.finalize(tikzfiles)
  #dev.off()

  return(c(val=uw.data$value, dval=uw.data$dvalue, tauint=uw.data$tauint, dtauint=uw.data$dtauint, Wopt=uw.data$Wopt))
}
