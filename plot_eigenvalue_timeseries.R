# function to plot timeseries of eigenvlues, including minimum and maximum eigenvalue bands
# as found in the monomial_0x.data files produced by tmLQCD

source("~/code/R/misc_R_scripts/plotutils/plot_multihist.R")

plot_eigenvalue_timeseries <- function(dat,trange,stepsize=1,pdf.filename,
                            ylab,plotsize,filelabel,titletext,errorband_color=rgb(0.6,0.0,0.0,0.6),
                            debug=FALSE) {
  xdat <- seq(trange[1],trange[2],stepsize)
  yrange <- range(dat[,2:5])
      
  uw.min_ev <- uwerrprimary(dat[,2])
  uw.max_ev <- uwerrprimary(dat[,3])
  
  if(debug){
    print("uw.eval.min_ev")
    summary(uw.min_ev)
    print("uw.eval.max_ev")
    summary(uw.max_ev)
  }

  #pdf(pdf.filename,width=plotsize,height=plotsize,title=filelabel)
  #op <- par(family="Palatino",cex.main=0.6,font.main=1)
  tikzfiles <- tikz.init(basename=pdf.filename,width=plotsize,height=plotsize)
  par(mgp=c(2,1,0))

  # plot the timeseries
  plot(x=xdat,xlim=trange,y=dat[,2],ylim=yrange,t='l',ylab=ylab,xlab=expression(t[MD]),main=titletext,log='y',tcl=0.02)
  lines(x=xdat,y=dat[,3])
 
  ## add the approximation interval
  lines(x=xdat,y=dat[,4],lty=2,col="darkgreen")
  lines(x=xdat,y=dat[,5],lty=2,col="darkgreen")
  
  # plot the corresponding histograms with error bands and mean values
  hist.min_ev <- hist(dat[,2],main=paste("min. eval",titletext),xlab="min. eval",tcl=0.02)
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

  tikz.finalize(tikzfiles)

  return(list(mineval=c(val=uw.min_ev$value, dval=uw.min_ev$dvalue, tauint=uw.min_ev$tauint, dtauint=uw.min_ev$dtauint, Wopt=uw.min_ev$Wopt),
              maxeval=c(val=uw.max_ev$value, dval=uw.max_ev$dvalue, tauint=uw.max_ev$tauint, dtauint=uw.max_ev$dtauint, Wopt=uw.max_ev$Wopt) ) )
              
}
