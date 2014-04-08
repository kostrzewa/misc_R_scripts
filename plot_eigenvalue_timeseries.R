# function to plot timeseries of eigenvlues, including minimum and maximum eigenvalue bands
# as found in the monomial_0x.data files produced by tmLQCD

plot_eigenvalue_timeseries <- function(dat,trange,stepsize=1,pdf.filename,
                            ylab,plotsize,filelabel,titletext,errorband_color=rgb(0.6,0.0,0.0,0.6)) {
  xdat <- seq(trange[1],trange[2],stepsize)
  yrange <- range(dat[,2:5])
      
  uw.min_ev <- uwerrprimary(dat[,2])
  uw.max_ev <- uwerrprimary(dat[,3])
  print("uw.eval.min_ev")
  summary(uw.min_ev)
  summary(uw.max_ev)

  pdf(pdf.filename,width=plotsize,height=plotsize,title=filelabel)
  op <- par(family="Palatino",cex.main=0.6,font.main=1)
  par(mgp=c(2,1,0))

  # plot the timeseries
  plot(x=xdat,xlim=trange,y=dat[,2],ylim=yrange,t='l',ylab=ylab,xlab=expression(t[HMC]),main=titletext,log='y')
  lines(x=xdat,y=dat[,3])
 
  ## add the approximation interval
  lines(x=xdat,y=dat[,4],lty=2,col="darkgreen")
  lines(x=xdat,y=dat[,5],lty=2,col="darkgreen")
  
  # plot the corresponding histograms with error bands and mean values
  hist.min_ev <- hist(dat[,2],main=paste("min. eval",titletext),xlab="min. eval")
  rect(ytop=max(hist.min_ev$counts),
       ybottom=0,
       xright=uw.min_ev$value+uw.min_ev$dvalue,
       xleft=uw.min_ev$value-uw.min_ev$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.min_ev$value,col="black")                                                                                                   

  hist.max_ev <- hist(dat[,3],main=paste("max. eval",titletext),xlab="max. eval")
  rect(ytop=max(hist.max_ev$counts),
       ybottom=0,
       xright=uw.max_ev$value+uw.max_ev$dvalue,
       xleft=uw.max_ev$value-uw.max_ev$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.max_ev$value,col="black")                                                                                                   

  dev.off()
}
