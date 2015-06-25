# function to plot familar mpcac versus 1/2kappa plot in a pretty fashion with legend
# the functon takes a data file as an argument and can be passed variable parameters 
# from the 'par' family which are in turn passed to the plotting function
# the legend position is determined from 'xlim' and 'ylim' automatically

# the data file must have a header containing (*without* quotes or #):
# "kappa mu mpcac dmpcac colour pch"
# with the corresponding values listed below, one set per row
# colour is a string for the colour name such as "red" or "blue"

plot_mpcac_v_kappa <- function(datafile,debug=F,sim=500,neg=TRUE,n.predict=5000,...)
{
  pcacdat <- read.table(file=datafile,header=T,stringsAsFactors=FALSE,fill=FALSE)
  pcacdat <- cbind(oneov2k=1/(2*pcacdat$kappa),pcacdat)
  if(debug) {
    print(pcacdat)
  }

  predict.newdata <- data.frame(oneov2k=seq(1/0.32,1/0.25,length.out=n.predict))

  
  models <- list()
  predictions <- list()
  coefs.cov <- list()
  coefs <- list()
  kappa_c <- data.frame(lwr=c(),fit=c(),upr=c())

  testpoint <- data.frame(x=0,y=0,dy=0)

  # if we have negative and positive mpcac masses (at least two of either),
  # we should do two fits because the negative and positive slopes are different
  no_mpcac <- c(length(which(pcacdat$mpcac>0)),length(which(pcacdat$mpcac<0)) )
  if(sum(no_mpcac)>=2){
    l.rws <- list()
    l.rws[[1]] <- 1:sum(no_mpcac)
    if(no_mpcac[1]>=2 | no_mpcac[2]>=2){
      l.rws[[1]] <- which(pcacdat$mpcac>0)
      l.rws[[2]] <- which(pcacdat$mpcac<0)
    }
    
    for( rws in l.rws ){
      if( length(rws) < 2 ) next;
      simdata <- data.frame(matrix(ncol=length(rws),nrow=sim))
      for( k in 1:length(rws) ) {
        simdata[,k] <- rnorm(n=sim,mean=pcacdat$mpcac[rws[k]],sd=pcacdat$dmpcac[rws[k]])
      }
      mpcacmod <- apply(X=simdata,MARGIN=1,
                        FUN=function(x) {
                              lm(mpcac~oneov2k,data=data.frame(mpcac=x,oneov2k=pcacdat$oneov2k[rws]),weights=1/(pcacdat$dmpcac[rws])^2)
                            } )

      coefs.tmp <- lapply(X=mpcacmod,FUN=function(x) { x$coefficients })
      coefs.tmp <- array(data=unlist(coefs.tmp),dim=c(2,sim))

      prediction <- lapply(X=mpcacmod,FUN=function(x) { predict(x,newdata=predict.newdata,interval="none") })

      prediction <- array(data=unlist(prediction),dim=c(n.predict,sim))

      models[[length(models)+1]] <- mpcacmod
      predictions[[length(predictions)+1]] <- data.frame(val=apply(X=prediction,MARGIN=1,FUN=mean),err=apply(X=prediction,MARGIN=1,FUN=sd))
      coefs[[length(coefs)+1]] <- apply(X=coefs.tmp,MARGIN=1,FUN=mean)
      coefs.cov[[length(coefs.cov)+1]] <- cov(t(coefs.tmp))

      kappa_estimate <- apply(X=prediction,MARGIN=2,FUN=function(x) { approx(x=x,y=predict.newdata$oneov2k,xout=0.0)$y } )
      kappa_c <- rbind(kappa_c, quantile(kappa_estimate,probs=c(0.1573,0.5,0.8427) ) )
    }
  } 
  if(nrow(kappa_c)>0) {
    cat("Estimates of kappa_c\n")
    kc <- 0.5/rev(kappa_c)
    print(sprintf("k_c = %f (+%f -%f)", kc[,2], kc[,3]-kc[,2], kc[,2]-kc[,1]))
    cat("\n")
  }
  
  cat("Fit coefficients and the sqrt of their variance-covariance matrices\n")
  for(idx in 1:length(coefs) ) {
    print(coefs[[idx]])
    print(sqrt(coefs.cov[[idx]]))
  }
  
  par(family="Times")
  
  # set up plot area
  plot(x=( 1/(2*pcacdat$kappa) ), y=pcacdat$mpcac,  
    xlab=expression(paste("1/2",kappa)), ylab=expression(~am[PCAC]), col=pcacdat$colour, 
    pch=pcacdat$pch+pcacdat$offsetpch, type='n', ... )
  
  if(nrow(kappa_c)>0){
    cols <- c("magenta","cyan")
    alpha.cols <- col2rgb(cols,alpha=TRUE)
    alpha.cols[4,] <- 40
    for( i in 1:nrow(kappa_c) ) {
      poly.x <- c(predict.newdata$oneov2k,rev(predict.newdata$oneov2k))

      coef.deriv <- matrix(nrow=2,ncol=length(predict.newdata$oneov2k))
      coef.deriv[1,] <- 1
      coef.deriv[2,] <- predict.newdata$oneov2k
      yvar <- t(coef.deriv) %*% coefs.cov[[i]] %*% coef.deriv
      dy <- sqrt(diag(yvar))

      poly.y <- c(predictions[[i]]$val-dy,rev(predictions[[i]]$val+dy))
      polygon(x=poly.x,y=poly.y,border=NA,
              col=rgb(red=alpha.cols[1,i],green=alpha.cols[2,i],blue=alpha.cols[3,i],alpha=alpha.cols[4,i],maxColorValue=255))
      points(y=0,x=kappa_c[i,2],col=cols[i],pch=16)
      lines(x=predict.newdata$oneov2k,y=predictions[[i]]$val,col=cols[i])
    }
  }

  abline(h=0,lty=2)
  
  plotwitherror(x=( 1/(2*pcacdat$kappa) ), y=pcacdat$mpcac, dy=pcacdat$dmpcac, 
    xlab=expression(paste("1/2",kappa)), ylab=expression(~am[PCAC]), col=pcacdat$colour, 
    pch=pcacdat$pch+pcacdat$offsetpch, rep=TRUE, ... )

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

  legend( x="topleft", col=unique(pcacdat$colour), legend=paste( "mu =", unique( pcacdat$mu ) ), 
          pch=unique(pcacdat$pch+pcacdat$offsetpch), bty='n' )
  
}
