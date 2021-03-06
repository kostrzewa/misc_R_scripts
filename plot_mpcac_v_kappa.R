# function to plot familar mpcac versus 1/2kappa plot in a pretty fashion with legend
# the functon takes a data file as an argument and can be passed variable parameters 
# from the 'par' family which are in turn passed to the plotting function
# the legend position is determined from 'xlim' and 'ylim' automatically

# the data file must have a header containing (*without* quotes or #):
# "kappa mu mpcac dmpcac colour pch"
# with the corresponding values listed below, one set per row
# colour is a string for the colour name such as "red" or "blue"

plot_mpcac_v_kappa <- function(datafile,
                               width=3.4, height=3.4,
                               basename = "mpcac_v_kappa",
                               debug=FALSE,neg=TRUE,fit=TRUE,
                               sim=500,n.predict=2000,
                               mpcac.ylim=c(-0.001,0.001),
                               mu.ylim=c(-0.001,0.001),
                               kappalegendpos="topleft",
                               ...)
{
  pcacdat <- read.table(file=datafile,header=T,stringsAsFactors=FALSE,fill=FALSE)
  pcacdat <- cbind(oneov2k=1/(2*pcacdat$kappa),pcacdat)
  if(debug) {
    print(pcacdat)
  }
  
  simdata <- NULL
  if(sim > 0) {
    simdata <- data.frame(matrix(ncol=nrow(pcacdat),nrow=sim))
    for( k in 1:nrow(pcacdat) ){
      simdata[,k] <- rnorm(n=sim,mean=pcacdat$mpcac[k],sd=pcacdat$dmpcac[k])
    }
  }
  
  lmodel <- NULL
  lmodel.sim <- NULL
  param.cov <- NULL
  pred.y <- NULL
  pred.x <- data.frame(oneov2k=seq(2.0,4.0,length.out=n.predict))
  if(fit) {
    lmodel <- lm(mpcac ~ oneov2k, data=pcacdat, weights=1/pcacdat$dmpcac^2)
    pred.y <- predict(lmodel,newdata=pred.x)
    print(summary(lmodel))
    if(sim>0){
      lmodel.sim <- t(apply(X=simdata,
                            MARGIN=1,
                            FUN=function(x){ 
                                    summary(lm(mpcac~oneov2k,
                                               data=data.frame(mpcac=x,
                                                               oneov2k=pcacdat$oneov2k),
                                               weights=1/pcacdat$dmpcac^2))$coefficients[1:2]
                                }
                         ) )
      param.cov <- cov(lmodel.sim)
      print(param.cov)
      print(sqrt(diag(param.cov)))
    }
  }

  tikzfiles <- tikz.init(basename=basename,width=width,height=height,lwdUnit=0.8)
  
  par(mgp=c(3,0.3,0))
  # prepare plot area
  plotwitherror(x=pcacdat$oneov2k, y=pcacdat$mpcac, dy=pcacdat$dmpcac, 
    xlab="", ylab="$ a m_\\mathrm{PCAC} $", col=pcacdat$colour, 
    pch=pcacdat$pch,las=1,tck=0.02, ylim=mpcac.ylim, t='n', ... )
  mtext(side=1,text="$ 1/2\\kappa $",line=1.3)
  abline(h=0,col="gray",lty=2)
  if(fit) {
    if(sim>0) {
      dpar <- as.matrix(data.frame(dint=rep(1,times=n.predict),dslope=pred.x$oneov2k))
      dy <- sqrt(diag( dpar %*% param.cov %*% t(dpar) ))
      poly.x <- c(pred.x$oneov2k,rev(pred.x$oneov2k))
      poly.y <- c(pred.y+dy,rev(pred.y-dy))
      polygon(x=poly.x,y=poly.y,col="#0000FF55",border=NA)
    }
    lines(x=pred.x$oneov2k,y=pred.y)
  }
  plotwitherror(x=pcacdat$oneov2k, y=pcacdat$mpcac, dy=pcacdat$dmpcac, 
                col=pcacdat$colour, pch=pcacdat$pch, rep=TRUE )
  lns <- rep(NA,length(unique(pcacdat$colour)))
  pch <- rep(15,length(unique(pcacdat$colour)))
  clrs <- unique(pcacdat$colour)

  lg <- sprintf("$a\\mu = %.6f$", unique( pcacdat$mu ))
  legend( x="topleft", col=clrs, legend=lg, 
          pch=pch, bty='n', pt.cex=1.3, ncol=2 )
  if(fit){
    legend(x = "topright",
           bty = 'n',
           pch = NA,
           lty = 1,
           col = "#0000FF55",
           legend = sprintf("$am_\\mathrm{PCAC}=%.3f+%.3f \\frac{1}{2\\kappa}$",lmodel$coefficients[1],lmodel$coefficients[2]))
    legend(x = "bottomleft",
           pch = NA,
           bty = 'n',
           legend = sprintf("$\\kappa_\\mathrm{cr} = %s$",
                   tex.catwitherror(x = -0.5*lmodel$coefficients[2]/lmodel$coefficients[1],
                                    dx = sd(-0.5*lmodel.sim[,2]/lmodel.sim[,1]),
                                    digits = 2,
                                    with.dollar = FALSE)
                   )
           )
  }
  legend( x="bottomright", col="black", legend=sprintf("$L/a = %d$", unique( pcacdat$L ) ), 
          pch=unique(pcacdat$pch), bty='n' )
  
  plot(x=pcacdat$mu, y=pcacdat$mpcac, dy=pcacdat$dmpcac, t='n', ylim=mu.ylim,
    xlab="", ylab="$ am_\\mathrm{PCAC} $", pch=pcacdat$pch,las=1,tck=0.02)
  abline(h=0,col="gray",lty=2)
  plotwitherror(rep=TRUE, x=pcacdat$mu, y=pcacdat$mpcac, dy=pcacdat$dmpcac,
    xlab="", ylab="$ am_\\mathrm{PCAC} $", col=pcacdat$kappacolour, 
    pch=pcacdat$pch )
  mtext(side=1,text="$ a\\mu_\\ell $",line=1.3)
  legend( x=kappalegendpos, col=unique(pcacdat$kappacolour), legend=sprintf("$\\kappa = %.7f$", unique( pcacdat$kappa )), 
          pch=15, bty='n', pt.cex=1.3 )
  legend( x="bottomright", col="black", legend=sprintf("$L/a = %d$", unique( pcacdat$L ) ), 
          pch=unique(pcacdat$pch), bty='n' )

  plotwitherror(x=pcacdat$oneov2k, y=pcacdat$mpi^2, dy=pcacdat$dmpi*pcacdat$mpi*2,
    xlab="", ylab="$ (a M_\\mathrm{PS} )^2 $", col=pcacdat$colour, 
    pch=pcacdat$pch,las=1,tck=0.02, ... )
  mtext(side=1,text="$ 1/2\\kappa $",line=1.3)
  legend( x="topright", col=unique(pcacdat$colour), legend=sprintf("$a\\mu = %.6f$", unique( pcacdat$mu )), 
          pch=15, bty='n', pt.cex=1.3 )
  legend( x="bottomright", col="black", legend=sprintf("$L/a = %d$", unique( pcacdat$L ) ), 
          pch=unique(pcacdat$pch), bty='n' )
  
  plotwitherror(x=pcacdat$oneov2k, y=pcacdat$fpi, dy=pcacdat$dfpi, 
    xlab="", ylab="$ a f_\\mathrm{PS} $", col=pcacdat$colour, 
    pch=pcacdat$pch,las=1,tck=0.02, ... )
  mtext(side=1,text="$ 1/2\\kappa $",line=1.3)
  legend( x="topright", col=unique(pcacdat$colour), legend=sprintf("$a\\mu = %.6f$", unique( pcacdat$mu )), 
          pch=15, bty='n', pt.cex=1.3 )
  legend( x="bottomright", col="black", legend=sprintf("$L/a = %d$", unique( pcacdat$L ) ), 
          pch=unique(pcacdat$pch), bty='n' )
  
  plotwitherror(x=pcacdat$mpcac, dx=pcacdat$dmpcac, y=pcacdat$mpi^2, dy=pcacdat$dmpi*pcacdat$mpi*2,
    xlab="", ylab="$ (a M_\\mathrm{PS} )^2 $", col=pcacdat$colour, 
    pch=pcacdat$pch,las=1,tck=0.02 )
  mtext(side=1,text="$ am_\\mathrm{PCAC} $",line=1.3)
  legend( x="topleft", col=unique(pcacdat$colour), legend=sprintf("$a\\mu = %.6f$", unique( pcacdat$mu )), 
          pch=15, bty='n', pt.cex=1.3 )
  legend( x="bottomright", col="black", legend=sprintf("$L/a = %d$", unique( pcacdat$L ) ), 
          pch=unique(pcacdat$pch), bty='n' )
  
  plotwitherror(x=pcacdat$mpcac, dx=pcacdat$dmpcac, y=pcacdat$fpi, dy=pcacdat$dfpi, 
    xlab="", ylab="$ a f_\\mathrm{PS} $", col=pcacdat$colour, 
    pch=pcacdat$pch,las=1,tck=0.02 )
  mtext(side=1,text="$ am_\\mathrm{PCAC} $",line=1.3)
  legend( x="topleft", col=unique(pcacdat$colour), legend=sprintf("$a\\mu = %.6f$", unique( pcacdat$mu )), 
          pch=15, bty='n', pt.cex=1.3 )
  legend( x="bottomright", col="black", legend=sprintf("$L/a = %d$", unique( pcacdat$L ) ), 
          pch=unique(pcacdat$pch), bty='n' )
  
  plotwitherror(x=pcacdat$mu, y=pcacdat$mpi^2, dy=pcacdat$dmpi*pcacdat$mpi*2,
    xlab="", ylab="$ (a M_\\mathrm{PS} )^2 $", col=pcacdat$kappacolour, 
    pch=pcacdat$pch,las=1,tck=0.02 )
  mtext(side=1,text="$ a \\mu_{\\ell} $",line=1.3)
  legend( x="topleft", col=unique(pcacdat$kappacolour), legend=sprintf("$\\kappa = %.7f$", unique( pcacdat$kappa )), 
          pch=15, bty='n', pt.cex=1.3 )
  legend( x="bottomright", col="black", legend=sprintf("$L/a = %d$", unique( pcacdat$L ) ), 
          pch=unique(pcacdat$pch), bty='n' )
  
  plotwitherror(x=pcacdat$mu, y=pcacdat$fpi, dy=pcacdat$dfpi,
    xlab="", ylab="$ a f_\\mathrm{PS} $", col=pcacdat$kappacolour, 
    pch=pcacdat$pch,las=1,tck=0.02 )
  mtext(side=1,text="$ a \\mu_{\\ell} $",line=1.3)
  legend( x="topleft", col=unique(pcacdat$kappacolour), legend=sprintf("$\\kappa = %.7f$", unique( pcacdat$kappa )), 
          pch=15, bty='n', pt.cex=1.3 )
  legend( x="bottomright", col="black", legend=sprintf("$L/a = %d$", unique( pcacdat$L ) ), 
          pch=unique(pcacdat$pch), bty='n' )
  
  plotwitherror(x=pcacdat$mu, y=pcacdat$P, dy=pcacdat$dP,
    xlab="", ylab="$ \\langle P \\rangle $", col=pcacdat$kappacolour, 
    pch=pcacdat$pch,las=1,tck=0.02 )
  mtext(side=1,text="$ a \\mu_{\\ell} $",line=1.3)
  legend( x="topright", col=unique(pcacdat$kappacolour), legend=sprintf("$\\kappa = %.7f$", unique( pcacdat$kappa )), 
          pch=15, bty='n', pt.cex=1.3 )
  legend( x="bottomright", col="black", legend=sprintf("$L/a = %d$", unique( pcacdat$L ) ), 
          pch=unique(pcacdat$pch), bty='n' )
  
  tikz.finalize(tikzfiles)
  
}
