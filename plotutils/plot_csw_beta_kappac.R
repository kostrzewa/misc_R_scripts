plot_csw_beta_kappac <- function(datafile, input,
                                 cswrange=c(1.0,2.5),betarange=c(1.6,1.8),kapparange=c(0.125,0.150),
                                 n.predict=1000, n.sim=200) {

  require("RColorBrewer")

  kcdat <- read.table(file=datafile,header=TRUE,colClasses="numeric")
  kcdat <- kcdat[order(kcdat$csw,kcdat$beta),]
  kcdat <- cbind(kcdat,logkappac=log(kcdat$kappac))
  kcdat <- cbind(kcdat,dlogkappac=kcdat$dkappac/kcdat$kappac)
  save(kcdat,file="csw_beta_kappac.Rdata")

  # prepare some simulated data with gaussian distributions
  simdata <- data.frame(matrix(ncol=length(kcdat$kappac),nrow=n.sim))
  for( k in 1:length(kcdat$kappac) ) {
    simdata[,k] <- log(rnorm(n=n.sim,mean=kcdat$kappac[k],sd=kcdat$dkappac[k]))
  }
  
  # csw,beta pairs for the prediction 
  newdata <- data.frame(csw=seq(cswrange[1]-0.5,cswrange[2]+0.5,length.out=n.predict),
                        beta=seq(betarange[1]-0.5,betarange[2]+0.5,length.out=n.predict))

  # model for log(kappa_c)
  kappacmodel <- lm(logkappac ~ csw + beta, data=kcdat, weights=1/kcdat$dlogkappac^2)
  print(anova(kappacmodel))
  print(summary(kappacmodel))
  # and confidence intervals which are calculated by removing points
  kc.confidence <- predict(kappacmodel,newdata=newdata,interval="confidence",level=0.6854)
  # bring into the same format as sim.confidence below
  kc.confidence <- kc.confidence[,c(2,1,3)]

  # and repeated on all the simulated data
  sim.kappacmodel <- apply(X=simdata,MARGIN=1,
                        FUN=function(x) {
                              lm(logkappac ~ csw + beta, 
                                 data=data.frame(logkappac=x,csw=kcdat$csw,beta=kcdat$beta), 
                                 weights=1/kcdat$dlogkappac^2)
                            } )
  
  # for the fits on the simulated data sets, do the prediction
  sim.prediction <- lapply(X=sim.kappacmodel, FUN=function(x) { predict(x,newdata=newdata,interval="none") } )
  sim.prediction <- array(data=unlist(sim.prediction),dim=c(n.predict,n.sim))
  
  # 68.57 confidence band for log(kappac) 
  sim.confidence <- t(apply(X=sim.prediction,MARGIN=1,FUN=quantile,probs=c(0.1573,0.5,0.8427)))
  
  # obtain some confidence bands for the fit coefficients 
  sim.coefs <- lapply(X=sim.kappacmodel, FUN=function(x) { x$coefficients })
  sim.coefs <- t(array(unlist(sim.coefs),dim=c(3,n.sim)))
  coefs.cov <- cov(sim.coefs)

  coefs <- apply(X=sim.coefs,MARGIN=2,quantile,probs=c(0.1573,0.5,0.8427))
  cat(sprintf("Intercept:    %f (+ %f - %f)\n", coefs[2,1], coefs[3,1]-coefs[2,1], coefs[2,1]-coefs[1,1]))
  cat(sprintf("alpha(csw):   %f (+ %f - %f)\n", coefs[2,2], coefs[3,2]-coefs[2,2], coefs[2,2]-coefs[1,2]))
  cat(sprintf("b(beta):      %f (+ %f - %f)\n", coefs[2,3], coefs[3,3]-coefs[2,3], coefs[2,3]-coefs[1,3]))
  print(sqrt(cov(sim.coefs)))
  
  # extract the median coefficients only 
  coefs <- coefs[2,]

  save(kappacmodel,file="kappac_csw_beta.model.Rdata")
  k.cfs <- kappacmodel$coefficients[2:3]

  # one symbol per csw value
  csw <- unique(kcdat$csw)
  syms <- c()
  sym <- 0
  for(i in csw){
    syms <- c(syms,rep(sym, times=length(which(kcdat$csw==i))) )
    sym <- sym + 1
  }
  
  # one colour per mass value
  mu <- unique(kcdat$mu)
  clr <- c()
  clr.idx <- 1
  clrs <- brewer.pal(n=9,name="Set1")
  for(i in mu){
    clr <- c(clr,rep(clrs[clr.idx], times=length(which(kcdat$mu==i))) )
    clr.idx <- clr.idx+1
    # we need to wrap around if there are too many
    if(clr.idx==9) clr.idx <- 1
  }

  tikzfiles <- tikz.init("csw_beta_kappac",width=3.5,height=3.5,lwdUnit=0.7)
  
  # the uncertainty in x comes from the slopes only
  #coef.deriv <- matrix(nrow=2,ncol=n.predict)
  #coef.deriv[1,] <- newdata$csw
  #coef.deriv[2,] <- newdata$beta
  #xvar <- t(coef.deriv) %*% coefs.cov[2:3,2:3] %*% coef.deriv
  #dx <- sqrt(diag(xvar))

  poly.x <- -(coefs[2]*newdata$csw+coefs[3]*newdata$beta)
  poly.x <- c(poly.x,rev(poly.x))

  # the uncertainty in y we can take from the simulated confidence interval
  # plus the confidence interval of the model (which includes factors such as the removal of points)
  
  dyp <- sqrt( (kc.confidence[,3]-kc.confidence[,2])^2 + (sim.confidence[,3]-sim.confidence[,2])^2 )
  dym <- sqrt( (kc.confidence[,2]-kc.confidence[,1])^2 + (sim.confidence[,2]-sim.confidence[,1])^2 ) 
   
  poly.y <- c(kc.confidence[,2]+dyp,rev(kc.confidence[,2]-dym))

    
  # and finally account for the uncertainty of the x position of the
  # data points (since it depends on the model)
  coef.deriv <- matrix(nrow=2,ncol=length(kcdat$csw))
  coef.deriv[1,] <- kcdat$csw
  coef.deriv[2,] <- kcdat$beta
  xvar <- t(coef.deriv) %*% coefs.cov[2:3,2:3] %*% coef.deriv
  dx <- sqrt(diag(xvar))
  xpts <- -(coefs[2]*kcdat$csw+coefs[3]*kcdat$beta)

  # compute plot limits based on the supplied beta and csw ranges
  xlims <- -c(coefs[2]*cswrange[1]+coefs[3]*betarange[1], coefs[2]*cswrange[2]+coefs[3]*betarange[2] )

  # prepare plot area
  plot(x=xpts,y=kcdat$kappac,#log='y',
       xlim=xlims, ylim=kapparange,
       xlab=sprintf("$ %.3f~c_\\mathrm{sw} + %.3f~\\beta $",-coefs[2],-coefs[3]),
       ylab="$ \\kappa_c $", 
       pch=syms, col=clr, type='n', las=1 )
  
  # draw error band and fit line
  polygon(x=poly.x,y=exp(poly.y),col="grey65", border=NA)
  lines(x=-c(coefs[2]*newdata$csw+coefs[3]*newdata$beta),y=exp(kc.confidence[,2]),lty=3)
 
  # to compare, let's try a simple linear fit
  #lin.kappacmodel <- lm(kappac ~ csw + beta, data=kcdat, weights=1/kcdat$dkappac^2)
  #lin.prediction <- predict(lin.kappacmodel,interval="none",newdata=newdata)
  #lines(y=lin.prediction,x=-(coefs[2]*newdata$csw+coefs[3]*newdata$beta),col="darkgreen",lty=3,lwd=2)
  
  # add points on top
  plotwitherror(x=xpts,y=kcdat$kappac,dy=kcdat$dkappac,
       xlim=xlims,ylim=kapparange,
       pch=syms, col=clr, rep=TRUE, las=1 )
  
  csw.legend <- sprintf("$ c_\\mathrm{sw} = %s $", csw)
  legend(x="bottomleft",legend=csw.legend,pch=unique(syms), bty='n')
  mu.legend <- sprintf("$ a\\mu = %s $", mu )
  legend(x="topright",legend=mu.legend,pch=15,col=unique(clr), bty='n',pt.cex=1.5)

  tikz.finalize(tikzfiles)

  if(!missing(input)){
    kc.confidence <- predict(kappacmodel,newdata=input,interval="confidence",level=0.6854)
    # bring into the same format as sim.confidence below
    kc.confidence <- exp(kc.confidence[,c(2,1,3)])
    # for the fits on the simulated data sets, do the prediction
    sim.prediction <- lapply(X=sim.kappacmodel, FUN=function(x) { predict(x,newdata=input,interval="none") } )
    sim.prediction <- array(data=unlist(sim.prediction),dim=c(length(input$csw),n.sim))
    # 68.54 confidence band for log(kappac) 
    sim.confidence <- exp(t(apply(X=sim.prediction,MARGIN=1,FUN=quantile,probs=c(0.1573,0.5,0.8427))))
    dyp <- sqrt( (kc.confidence[,3]-kc.confidence[,2])^2 + 
                 (sim.confidence[,3]-sim.confidence[,2])^2 
               )
    dym <- sqrt( (kc.confidence[,2]-kc.confidence[,1])^2 + 
                 (sim.confidence[,2]-sim.confidence[,1])^2 
               )

    for( i in 1:length(input$beta) ){
      cat(sprintf("kappa_c(beta=%f,csw=%f): %f (+ %f - %f)\n", input$beta[i], input$csw[i],
                                            kc.confidence[i,2], dyp[i], dym[i] ) 
         )
    }
  }
}
  
