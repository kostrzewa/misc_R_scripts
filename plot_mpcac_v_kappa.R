# function to plot familar mpcac versus 1/2kappa plot in a pretty fashion with legend
# the functon takes a data file as an argument and can be passed variable parameters 
# from the 'par' family which are in turn passed to the plotting function
# the legend position is determined from 'xlim' and 'ylim' automatically

# the data file must have a header containing (*without* quotes or #):
# "kappa mu mpcac dmpcac colour pch"
# with the corresponding values listed below, one set per row
# colour is a string for the colour name such as "red" or "blue"

plot_mpcac_v_kappa <- function(datafile,interval,debug=F,...)
{
  pcacdat <- read.table(file=datafile,header=T,stringsAsFactors=FALSE)
  pcacdat <- cbind(oneov2k=1/(2*pcacdat$kappa),pcacdat)
  if(debug) {
    print(pcacdat)
  }
  
  if(!missing(interval)){
    mpcacmod <- lm(mpcac~oneov2k, data=pcacdat, weights=(1/pcacdat$dmpcac)^2)
    rootfun <- function(x,coefs) { coefs[1] + coefs[2]*x }
    kappa_c <- uniroot(f=rootfun,interval=interval,coefs=mpcacmod$coefficients)
    cat("estimate of kappa_c: ", 0.5*(1/kappa_c$root),"\n")
  }
  
  par(family="Times")

  plotwitherror(x=( 1/(2*pcacdat$kappa) ), y=pcacdat$mpcac, dy=pcacdat$dmpcac, 
    xlab=expression(paste("1/2",kappa)), ylab=expression(~am[PCAC]), col=pcacdat$colour, 
    pch=pcacdat$pch+pcacdat$offsetpch, ... )

  abline(h=0,lty=2)

  # get plot boundaries
  lims <- par("usr")

  # attempt to extract some coordinates for the legend 
  # from variable parameter list
  var_params <- list(...)
  if(debug){
    print(var_params)
  }
  
  if( 'xlim' %in% names(var_params) ){
    legend.xpos <- var_params$xlim[1]
  } else {
    print("legend x position NOT set")
    legend.xpos <- 0
  }

  if( 'ylim' %in% names(var_params) ){
    legend.ypos <- var_params$ylim[2]
  } else {
    print("legend y position NOT set")
    legend.ypos <- lims[4]
  }

  legend( x=legend.xpos, y=legend.ypos, col=unique(pcacdat$colour), legend=paste( "mu =", unique( pcacdat$mu ) ), 
          pch=unique(pcacdat$pch+pcacdat$offsetpch) )
  
}
