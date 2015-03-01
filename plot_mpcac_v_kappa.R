# function to plot familar mpcac versus 1/2kappa plot in a pretty fashion with legend
# the functon takes a data file as an argument and can be passed variable parameters 
# from the 'par' family which are in turn passed to the plotting function
# the legend position is determined from 'xlim' and 'ylim' automatically

# the data file must have a header containing (*without* quotes or #):
# "kappa mu mpcac dmpcac colour pch"
# with the corresponding values listed below, one set per row
# colour is a string for the colour name such as "red" or "blue"


plot_mpcac_v_kappa <- function(datafile,debug=F,...)
{
  pcacdat <- read.table(file=datafile,header=T,stringsAsFactors=FALSE,fill=TRUE)
  pcacdat <- cbind(oneov2k=1/(2*pcacdat$kappa),pcacdat)
  if(debug) {
    print(pcacdat)
  }

  predict.newdata <- data.frame(oneov2k=1/(2*seq(0.1,0.2,length.out=200)))

  
  models <- list()
  predictions <- list()
  kappa_c <- data.frame(min=c(),med=c(),max=c())

  # if we have negative and positive mpcac masses (at least two of either),
  # we should do two fits because the negative and positive slopes are different
  no_mpcac <- c(length(which(pcacdat$mpcac>0)),length(which(pcacdat$mpcac<0)) )
  sign <- c(1,-1)
  if(no_mpcac[1]>=2 | no_mpcac[2]>=2){
    for( i in 1:2 ){
      if( no_mpcac[i] < 2 ) next;
      rws <- which(sign[i]*pcacdat$mpcac>0)
      mpcacmod <- lm(mpcac~oneov2k, data=pcacdat[rws,], weights=(1/pcacdat$dmpcac[rws])^2)
      prediction <- predict(mpcacmod,newdata=predict.newdata,interval="confidence",level=0.68)
      models[[length(models)+1]] <- mpcacmod
      predictions[[length(predictions)+1]] <- prediction
      # interpolate the model and the confidence bands (if available) to find 1/2k corresponding to mpcac=0.0
      lwr <- NA
      upr <- NA
      if(!any(is.na(prediction[,2]))) {
        lwr <- approx(x=prediction[,2],y=predict.newdata$oneov2k,xout=0.0)$y
        upr <- approx(x=prediction[,3],y=predict.newdata$oneov2k,xout=0.0)$y
      }
      kappa_c <- rbind(kappa_c, 
                       data.frame( lwr=lwr,
                                   fit=approx(x=prediction[,1],y=predict.newdata$oneov2k,xout=0.0)$y,
                                   upr=upr ) )
    }
  } else if(sum(no_mpcac) >= 2) {
    mpcacmod <- lm(mpcac~oneov2k, data=pcacdat, weights=(1/pcacdat$dmpcac)^2)
    models[[length(models)+1]] <- mpcacmod
    prediction <- predict(mpcacmod,newdata=predict.newdata,interval="confidence",level=0.68)
    predictions[[length(predictions)+1]] <- prediction
    # interpolate the model and the confidence bands (if available) to find 1/2k corresponding to mpcac=0.0
    lwr <- NA
    upr <- NA
    if(!any(is.na(prediction[,2]))) {
      lwr <- approx(x=prediction[,2],y=predict.newdata$oneov2k,xout=0.0)$y
      upr <- approx(x=prediction[,3],y=predict.newdata$oneov2k,xout=0.0)$y
    }
    kappa_c <- rbind(kappa_c, 
                     data.frame( lwr=lwr,
                                 fit=approx(x=prediction[,1],y=predict.newdata$oneov2k,xout=0.0)$y,
                                 upr=upr ) )
  }
  if(nrow(kappa_c)>0) {
    cat("Estimates of kappa_c\n")
    print(0.5*(1/kappa_c))
    cat("\n")
  }
   
  par(family="Times")

  plotwitherror(x=( 1/(2*pcacdat$kappa) ), y=pcacdat$mpcac, dy=pcacdat$dmpcac, 
    xlab=expression(paste("1/2",kappa)), ylab=expression(~am[PCAC]), col=pcacdat$colour, 
    pch=pcacdat$pch+pcacdat$offsetpch, ... )

  if(nrow(kappa_c)>0){
    cols <- c("magenta","cyan")
    alpha.cols <- col2rgb(cols,alpha=TRUE)
    alpha.cols[4,] <- 40
    for( i in 1:nrow(kappa_c) ) {
      if(!any(is.na(predictions[[i]][,2]))){
        poly.x <- c(predict.newdata$oneov2k,rev(predict.newdata$oneov2k))
        poly.y <- c(predictions[[i]][,2],rev(predictions[[i]][,3]))
        polygon(x=poly.x,y=poly.y,border=NA,
                col=rgb(red=alpha.cols[1,i],green=alpha.cols[2,i],blue=alpha.cols[3,i],alpha=alpha.cols[4,i],maxColorValue=255))
      }
      points(y=0,x=kappa_c$fit[i],col=cols[i],pch=16)
      abline(models[[i]],col=cols[i])
    }
  }

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
