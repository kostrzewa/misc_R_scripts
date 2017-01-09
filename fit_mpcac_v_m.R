# function to plot familar mpcac versus 1/2kappa plot in a pretty fashion with legend
# the functon takes a data file as an argument and can be passed variable parameters 
# from the 'par' family which are in turn passed to the plotting function
# the legend position is determined from 'xlim' and 'ylim' automatically

# the data file must have a header containing (*without* quotes or #):
# "kappa mu mpcac dmpcac colour pch"
# with the corresponding values listed below, one set per row
# colour is a string for the colour name such as "red" or "blue"

fit_mpcac_v_m <- function(datafile,kappac.only=TRUE,n.param=2,debug=FALSE,sim=500,neg=TRUE,n.predict=5000,draw.pred.1=TRUE,...)
{
  pcacdat <- read.table(file=datafile,header=T,stringsAsFactors=FALSE,fill=FALSE)
  pcacdat <- cbind(oneov2k=1/(2*pcacdat$kappa),pcacdat)
  if(debug) {
    print(pcacdat)
  }
  
  # attempt to determine c2*a^2/B
  if(!is.null(pcacdat$kappac) && kappac.only==FALSE ) {
    m <- 1/(2*pcacdat$kappa) - 1/(2*pcacdat$kappac)
    dm <- pcacdat$dkappac/(2*pcacdat$kappac^2)
  
    mpcacdat <- data.frame(mpcac=pcacdat$mpcac,dmpcac=pcacdat$dmpcac,m=m,dm=dm,mu=pcacdat$mu)
    mpcacmodel <- NULL
    mpcacmodel.sim <- NULL
    if(n.param==5){
      mpcacmodel <- nls(mpcac ~ b*(m+d*mu^2) - c*( b*(m+d*mu^2) / sqrt( (b*(m+d*mu^2))^2 + (e*(mu + f*m))^2) ), start=list(b=1,c=-1,d=0.01,e=0.01,f=0.01),
                        data=mpcacdat,weights=1/pcacdat$dmpcac^2, model=TRUE,algorithm='port',trace=TRUE,control=list(maxiter=10000))
    } else if(n.param==3){
      mpcacmodel <- nls(mpcac ~ b*(m+d*mu) - c*( (m+d*mu) / sqrt( (m+d*mu)^2+mu^2) ) ,start=list(b=1,c=-1,d=0.01),
                        data=mpcacdat,weights=1/pcacdat$dmpcac^2, model=TRUE,algorithm='port',trace=TRUE,control=list(maxiter=10000))
      #mpcacmodel <- nls(mpcac ~ a*m - c*( m / sqrt( m^2+mu^2)) + d*mu ,start=list(a=1,c=-1,d=0.01),
      #                  data=mpcacdat,weights=1/pcacdat$dmpcac^2, model=TRUE,algorithm='port',trace=TRUE,control=list(maxiter=10000))
    } else if(n.param==2){
      mpcacmodel <- nls(mpcac ~ a*m - c*( m / sqrt(m^2+mu^2) ) ,start=list(a=1,c=-1),
                        data=mpcacdat,weights=1/pcacdat$dmpcac^2, model=TRUE,algorithm='port',trace=TRUE,control=list(maxiter=10000))
    } else {
      mpcacmodel <- nls(mpcac ~ m - c*( m / sqrt(m^2+mu^2) ) ,start=list(c=-1),data=mpcacdat,weights=1/pcacdat$dmpcac^2, model=TRUE,algorithm='port',trace=TRUE,control=list(maxiter=10000))
    }

    #mpcaclinmodel <- lm(mpcac ~ m - 1 ,data=mpcacdat,weights=1/pcacdat$dmpcac^2,model=TRUE)
    
    if(sim>0){
      cat("Computing errors via simulations\n")
       
      # repeat fits on simulated data set with representative errors 
      simdata <- list()
      for( i in 1:sim ){
        simdata[[length(simdata)+1]] <- data.frame(m=rnorm(mean=m,sd=dm,n=length(m)),
                                                   mpcac=rnorm(mean=pcacdat$mpcac,sd=pcacdat$dmpcac,n=length(pcacdat$mpcac)),
                                                   mu=pcacdat$mu
                                                  )
      }

      mpcacmodel.sim <- lapply(X=simdata,FUN=function(x) {
                                                          if(n.param==5){
                                                             summary(nls(mpcac ~ b*(m+d*mu^2) - c*( b*(m+d*mu^2) / sqrt( (b*(m+d*mu^2))^2 + (e*(mu + f*m))^2) ), 
                                                                         start=list(b=1,c=-1,d=0.01,e=0.01,f=0.01),
                                                                         data=x,weights=1/pcacdat$dmpcac^2, model=TRUE,
                                                                         algorithm='port',
                                                                         trace=TRUE,control=list(maxiter=10000))
                                                                    )$coefficients[,1]
                                                          } else if(n.param==3){ 
                                                            #summary(nls(mpcac ~ a*(m+d*mu) - c*( a*(m+d*mu) / sqrt( (a*(m+d*mu))^2+mu^2) ) ,
                                                            #            start=list(a=1,c=-1,d=0.01),data=x,weights=1/pcacdat$dmpcac^2,algorithm='port',
                                                            #            control=list(maxiter=10000))
                                                            #       )$coefficients[,1] 
                                                            summary(nls(mpcac ~ b*(m+d*mu) - c*( (m+d*mu) / sqrt( (m+d*mu)^2+mu^2 ) ) ,
                                                                        start=list(b=1,c=-1,d=0.01),data=x,weights=1/pcacdat$dmpcac^2,algorithm='port',
                                                                        control=list(maxiter=10000))
                                                                   )$coefficients[,1] 
                                                          } else if (n.param==2) {
                                                            summary(nls(mpcac ~ a*(m+d*mu) - c*( (m+d*mu) / sqrt((m+d*mu)^2+mu^2) ) ,
                                                                        start=list(a=1,c=-1,d=0.01),data=x,weights=1/pcacdat$dmpcac^2,algorithm='port')
                                                                   )$coefficients[,1]
                                                          } else {
                                                            summary(nls(mpcac ~ m - c*( m / sqrt(m^2+mu^2) ) ,
                                                                        start=list(c=-1),data=x,weights=1/pcacdat$dmpcac^2)
                                                                   )$coefficients[,1] 
                                                          }
                                                         } 
                              )
      mpcacmodel.sim <- t(array(unlist(mpcacmodel.sim),dim=c(length(mpcacmodel.sim[[1]]),sim)))
      
      #mpcaclinmodel.sim <- lapply(X=simdata,FUN=function(x) { summary(lm(mpcac ~ m - 1,data=x,weights=1/pcacdat$dmpcac^2,model=TRUE))$coefficients[1] } )
      #mpcaclinmodel.sim <- t(array(unlist(mpcaclinmodel.sim),dim=c(length(mpcaclinmodel.sim[[1]]),sim)))

      #print(mpcacmodel.sim)
      #print(mpcaclinmodel.sim)
      if(n.param==2){
        # two param
        cat(sprintf("b=%.8f(%.8f) c=%.8f(%.8f)\n",mean(mpcacmodel.sim[,1]),sd(mpcacmodel.sim[,1]),mean(mpcacmodel.sim[,2]),sd(mpcacmodel.sim[,2])))
      } else {
        # one param
        cat(sprintf("c=%.8f(%.8f)\n",mean(mpcacmodel.sim[,1]),sd(mpcacmodel.sim[,1])))#,mean(mpcacmodel.sim[,2]),sd(mpcacmodel.sim[,2])))
      }
      mpcacmodel.cov <- cov(mpcacmodel.sim)
      print(sqrt(mpcacmodel.cov))
    }
    #cat(sprintf("b=%.8f(%.8f)\n",mean(mpcaclinmodel.sim),sd(mpcaclinmodel.sim)))
    
    #print(mpcacmodel)
    #print(summary(mpcacmodel))
    cat("\n##############################################\n")
    cat(" ChiPT weighted chisq/df:", sum(summary(mpcacmodel)$residuals^2)/summary(mpcacmodel)$df[2], "\n")
    #print(summary(mpcacmodel)$residuals)
    print(summary(mpcacmodel)$coefficients)
    #cat("\n Linear model weighted chisq/df:", sum(summary(mpcaclinmodel)$residuals^2)/summary(mpcaclinmodel)$df[2], "\n")
    #print(summary(mpcaclinmodel)$residuals)
    #print(summary(mpcaclinmodel)$coefficients)
    #require("propagate")
    cat("##############################################\n\n")
  
    b <- summary(mpcacmodel)$coefficients[1,1]
    c <- summary(mpcacmodel)$coefficients[2,1]
    d <- summary(mpcacmodel)$coefficients[3,1]
    dpar <- function(b,c,d,mu,m){
      as.matrix(data.frame(db=m+d*mu,
                           dc=-(m+d*mu)/sqrt((m+d*mu)^2+mu^2),
                           dd=b*mu + c*mu*(d*mu+m)^2/((d*mu+m)^2+mu^2)^(3/2)-
                              c*mu/sqrt((d*mu+m)^2+mu^2)))
    }
    pred.m <- seq(-0.1,0.1,0.00005)
    newdata.1 <- data.frame(m=pred.m,mu=rep(0.0009,times=length(pred.m)))
    newdata.2 <- data.frame(m=pred.m,mu=rep(0.003,times=length(pred.m)))
    newdata.3 <- data.frame(m=pred.m,mu=rep(0.006,times=length(pred.m)))
    newdata.4 <- data.frame(m=pred.m,mu=rep(0.0,times=length(pred.m)))

    #pred <- predictNLS(mpcacmodel,newdata=newdata,interval="confidence",alpha=0.3146,do.sim=FALSE)
    pred.1 <- predict(mpcacmodel,newdata=newdata.1,interval="none")
    pred.2 <- predict(mpcacmodel,newdata=newdata.2,interval="none")
    pred.3 <- predict(mpcacmodel,newdata=newdata.3,interval="none")
    #pred.lin <- predict(mpcaclinmodel,newdata=data.frame(m=pred.m),interval='none')
    
    if(sim>0){
      # compute errors
      dpar.3 <- NULL
      dpar.2 <- NULL
      dpar.1 <- NULL
      if(n.param==3){
        dpar.1 <- dpar(b=b,c=c,d=d,m=pred.m,mu=newdata.1$mu)
        dpar.2 <- dpar(b=b,c=c,d=d,m=pred.m,mu=newdata.2$mu)
        dpar.3 <- dpar(b=b,c=c,d=d,m=pred.m,mu=newdata.3$mu)
      } else {
        dpar.2 <- as.matrix(data.frame(dc=-pred.m/sqrt(pred.m^2+newdata.2$mu^2)))
      }
      pred.dy.2 <- abs(sqrt(diag(dpar.2 %*% mpcacmodel.cov %*% t(dpar.2))))
      poly.x.2 <- c(pred.m,rev(pred.m))
      poly.y.2 <- c(pred.2+pred.dy.2,rev(pred.2-pred.dy.2))
      
      pred.dy.1 <- abs(sqrt(diag(dpar.1 %*% mpcacmodel.cov %*% t(dpar.1))))
      poly.x.1 <- c(pred.m,rev(pred.m))
      poly.y.1 <- c(pred.1+pred.dy.1,rev(pred.1-pred.dy.1))
      
      pred.dy.3 <- abs(sqrt(diag(dpar.3 %*% mpcacmodel.cov %*% t(dpar.3))))
      poly.x.3 <- c(pred.m,rev(pred.m))
      poly.y.3 <- c(pred.3+pred.dy.3,rev(pred.3-pred.dy.3))
    }

    #pdf("ampcac_v_am.pdf")
    tikzfiles <- tikz.init(basename="ampcac_v_am",width=4.8,height=4.5,lwdUnit=0.8)
    
    lims <- data.frame(ymin=c(-0.05,-0.015,-0.0015),ymax=c(0.05,0.015,0.0015),xmin=c(-0.05,-0.007,-0.0015),xmax=c(0.05,0.007,0.0015))
    
    for( r in 1:nrow(lims) ){
    # prepare plot region
      plotwitherror(x=m,dx=dm,y=pcacdat$mpcac,dy=pcacdat$dmpcac,xlab="$am_0-am_\\mathrm{cr}$",ylab="$am_\\mathrm{PCAC}$",t='n',
                    ylim=c(lims$ymin[r],lims$ymax[r]),
                    xlim=c(lims$xmin[r],lims$xmax[r])
                   )
      #for(i in 1:sim){
      #  points(x=simdata[[i]]$m,y=simdata[[i]]$mpcac,pch='.')
      #}
      #lines(x=newdata$m,y=pred$summary[,1])
      if(sim>0) polygon(x=poly.x.2,y=poly.y.2,col="grey",border=NA)
      if(draw.pred.1){
        if(sim>0) polygon(x=poly.x.1,y=poly.y.1,col="#0099FF55",border=NA)
        lines(x=newdata.1$m,y=pred.1,lty=2,col="blue",lwd=1.2)
        
        # for mu==0, the discontinuity means that the polygon can't be drawn properly...
        if(sim>0) polygon(x=poly.x.3,y=poly.y.3,col="#00FF0055",border=NA)
        lines(x=newdata.3$m,y=pred.3,lty=4,col="green",lwd=1.2)
      }
      lines(x=newdata.2$m,y=pred.2,lty=1,lwd=1.2)
      #lines(x=pred.m,y=pred.lin,col="red",lty=3,lwd=1.2)
  
      #lines(x=newdata.3$m,y=pred.3,lty=3)
      #lines(x=newdata$m,y=pred$summary[,5])
      #lines(x=newdata$m,y=pred$summary[,6])
      abline(h=0,lty=3)
      abline(v=0,lty=3)
      plotwitherror(x=m,dx=dm,y=pcacdat$mpcac,dy=pcacdat$dmpcac,pch=pcacdat$pch,col=pcacdat$colour,rep=TRUE,cex=1.3)
      #legend with line labels, looks a bit bad because the labels are intersected by the 0,0 guides
      #legend(x="topleft",
      #       legend=c(sprintf("$a\\mu_\\ell = %.4f$",unique(pcacdat$mu)),sprintf("$b\\cdot am-c\\cdot ( am / \\sqrt{ (am)^2 + (%.4f)^2 } )$",c(0.003,0.0009)),"$d\\cdot am$"),
      #       col=c(unique(pcacdat$colour),"black","blue","red"),
      #       pch=c(rep(15,length(unique(pcacdat$mu))),NA,NA,NA),
      #       lty=c(rep(NA,length(unique(pcacdat$mu))),1,2,3),
      #       bty='n',pt.cex=1.5,lwd=2)
      legend(x="topleft",legend=sprintf("$a\\mu_\\ell = %.4f$",unique(pcacdat$mu)),col=unique(pcacdat$colour),pch=15,bty='n',pt.cex=1.5)
      legend(x="bottomright",legend=sprintf("$L/a=%d$",unique(pcacdat$L)),col="black",pch=unique(pcacdat$pch),bty='n')
    }
    
    plotwitherror(y=pcacdat$mpi^2,dy=2*pcacdat$mpi*pcacdat$dmpi,x=pcacdat$mpcac,dx=pcacdat$dmpcac,ylab="$(aM_{\\pi^\\pm})^2$",xlab="$m_\\mathrm{PCAC}$")
    plotwitherror(y=pcacdat$fpi,dy=pcacdat$dfpi,x=pcacdat$mpcac,dx=pcacdat$dmpcac,ylab="$f_{\\pi^\\pm}$",xlab="$m_\\mathrm{PCAC}$")
    #plotwitherror(y=pcacdat$P,dy=pcacdat$dP,x=pcacdat$mpcac,dx=pcacdat$dmpcac,ylab="$\\langle P \\rangle$",xlab="$m_\\mathrm{PCAC}$")

    tikz.finalize(tikzfiles)
    
    readline()
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
    if(neg==FALSE && (no_mpcac[1]>=2 | no_mpcac[2]>=2)){
      l.rws[[1]] <- which(pcacdat$mpcac>0)
      l.rws[[2]] <- which(pcacdat$mpcac<0)
    }
    
    for( rws in l.rws ){
      print(rws)
      if( length(rws) < 2 ) next;
      simdata <- data.frame(matrix(ncol=length(rws),nrow=sim))
      for( k in 1:length(rws) ) {
        print("simulating")
        simdata[,k] <- rnorm(n=sim,mean=pcacdat$mpcac[rws[k]],sd=pcacdat$dmpcac[rws[k]])
      }
      # TODO: add estimate of systematic error due to point removal if there are enough degrees of freedom to do so
      mpcacmod <- lm(mpcac~oneov2k,data=pcacdat[rws,],weights=1/(pcacdat$dmpcac[rws])^2)
      prediction <- predict(mpcacmod,newdata=predict.newdata,interval="confidence",level=0.6854)
      pred <- data.frame(val=prediction[,1],mdval=prediction[,1]-prediction[,2],dval=prediction[,3]-prediction[,2])

      sim.mpcacmod <- apply(X=simdata,MARGIN=1,
                        FUN=function(x) {
                              lm(mpcac~oneov2k,data=data.frame(mpcac=x,oneov2k=pcacdat$oneov2k[rws]),weights=1/(pcacdat$dmpcac[rws])^2)
                            } )

      coefs.tmp <- lapply(X=sim.mpcacmod,FUN=function(x) { x$coefficients })
      coefs.tmp <- array(data=unlist(coefs.tmp),dim=c(2,sim))

      sim.prediction <- lapply(X=sim.mpcacmod,FUN=function(x) { predict(x,newdata=predict.newdata,interval="none") })

      sim.prediction <- array(data=unlist(sim.prediction),dim=c(n.predict,sim))

      models[[length(models)+1]] <- sim.mpcacmod
      predictions[[length(predictions)+1]] <- data.frame(val=apply(X=sim.prediction,MARGIN=1,FUN=median),
                                                         mdval = sqrt( apply(X=sim.prediction,MARGIN=1,FUN=sd)^2  +
                                                                       pred$mdval^2 ),
                                                         dval = sqrt( apply(X=sim.prediction,MARGIN=1,FUN=sd)^2   +
                                                                      pred$dval^2 ) 
                                                        )
                                                         
      coefs[[length(coefs)+1]] <- apply(X=coefs.tmp,MARGIN=1,FUN=mean)
      coefs.cov[[length(coefs.cov)+1]] <- cov(t(coefs.tmp))

      print("sim calculating kappa_estimate")
      sim.kappa_estimate <- apply(X=sim.prediction,MARGIN=2,FUN=function(x) { approx(x=x,y=predict.newdata$oneov2k,xout=0.0)$y } )
      print("calculating kappa_estimate")
      kappa_estimate <- approx(x=prediction[,1],y=predict.newdata$oneov2k,xout=0.0)
      print("test")
      print(kappa_estimate)
      #stop()
      sim.kappa_c <- quantile(sim.kappa_estimate,probs=c(0.1573,0.5,0.8427))
      kappa_c <- rbind(kappa_c, quantile(sim.kappa_estimate,probs=c(0.1573,0.5,0.8427) ) )
    }
  } 
  if(nrow(kappa_c)>0) {
    cat("Estimates of kappa_c\n")
    kc <- 0.5/rev(kappa_c)
    print(sprintf("k_c = %.8f (+%.8f -%.8f)", kc[,2], kc[,3]-kc[,2], kc[,2]-kc[,1]))
    cat("\n")
  }
  
  cat("Fit coefficients and the sqrt of their variance-covariance matrices\n")
  for(idx in 1:length(coefs) ) {
    print(coefs[[idx]])
    print(sqrt(coefs.cov[[idx]]))
  }
  
  pdf("mpcac_v_kappa.pdf")
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

  dev.off()
  
}
