# function to plot familar mpcac versus 1/2kappa plot in a pretty fashion with legend
# the functon takes a data file as an argument and can be passed variable parameters 
# from the 'par' family which are in turn passed to the plotting function
# the legend position is determined from 'xlim' and 'ylim' automatically

# the data file must have a header containing (without quotes or #):
# "kappa mu mpcac dmpcac colour"
# with the corresponding values listed below, one set per row

plot_mpcac_v_kappa <- function(datafile,debug=F,...)
{
  pcacdat <- read.table(file=datafile,header=T)
  if(debug) {
    print(pcacdat)
  }

  par(family="Times")

  plotwitherror(x=( 1/(2*pcacdat$kappa) ), y=pcacdat$mpcac, dy=pcacdat$dmpcac, 
    xlab=expression(paste("1/2",kappa)), ylab=expression(~am[PCAC]), col=pcacdat$colour, ... )

  abline(h=0,lty=2)

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
    legend.ypos <- 0
  }

  legend( x=legend.xpos, y=legend.ypos, col=unique(pcacdat$colour), legend=paste( "mu =", unique( pcacdat$mu ) ), pch=1 )

}
