plot_csw_beta_a <- function(datafile, input,
                                 cswrange=c(1.0,2.5),betarange=c(1.6,1.8),arange=c(0.05,0.15),
                                 n.predict=1000, n.sim=200) {

  require("RColorBrewer")

  adat <- read.table(file=datafile,header=TRUE,colClasses="numeric")
  adat <- adat[order(adat$csw,adat$beta),]
  adat <- cbind(adat,loga=log(adat$a))
  adat <- cbind(adat,dloga=adat$da/adat$a)
  save(adat,file="csw_beta_a.Rdata")

  # for the fit we take into account the number of trajectories used for the estimate of w0
  # the inverse error on a and the inverse absolute value of ampcac, so that 
  # critical ensembles are assigned a larger weight
  #weights <- adat$Nt/abs(adat$mpcac)/adat$dloga^2
  weights <- (adat$Nt/max(adat$Nt))*(max(abs(adat$mpcac))/abs(adat$mpcac))/adat$dloga^2

  # prepare some simulated data with gaussian distributions
  simdata <- data.frame(matrix(ncol=length(adat$a),nrow=n.sim))
  for( k in 1:length(adat$a) ) {
    simdata[,k] <- log(rnorm(n=n.sim,mean=adat$a[k],sd=adat$da[k]))
  }
  
  # csw,beta pairs for the prediction 
  newdata <- data.frame(csw=seq(cswrange[1]-0.5,cswrange[2]+0.5,length.out=n.predict),
                        beta=seq(betarange[1]-0.5,betarange[2]+0.5,length.out=n.predict))

  # model for log(a)
  amodel <- lm(loga ~ csw + beta, data=adat, weights=weights)
  print(anova(amodel))
  
  # and confidence intervals which are calculated by removing points
  a.confidence <- predict(amodel,newdata=newdata,interval="confidence",level=0.6854)
  # bring into the same format as sim.confidence below
  a.confidence <- a.confidence[,c(2,1,3)]

  # and repeated on all the simulated data
  sim.amodel <- apply(X=simdata,MARGIN=1,
                        FUN=function(x) {
                              lm(loga ~ csw + beta, 
                                 data=data.frame(loga=x,csw=adat$csw,beta=adat$beta), 
                                 weights=weights)
                            } )
  
  # for the fits on the simulated data sets, do the prediction
  sim.prediction <- lapply(X=sim.amodel, FUN=function(x) { predict(x,newdata=newdata,interval="none") } )
  sim.prediction <- array(data=unlist(sim.prediction),dim=c(n.predict,n.sim))
  
  # 68.57 confidence band for log(a) 
  sim.confidence <- t(apply(X=sim.prediction,MARGIN=1,FUN=quantile,probs=c(0.1573,0.5,0.8427)))
  
  # obtain some confidence levels for the fit coefficients 
  sim.coefs <- lapply(X=sim.amodel, FUN=function(x) { x$coefficients })
  sim.coefs <- t(array(unlist(sim.coefs),dim=c(3,n.sim)))
  # we will use the variance-covariance matrix below
  coefs.cov <- cov(sim.coefs)

  coefs <- apply(X=sim.coefs,MARGIN=2,quantile,probs=c(0.1573,0.5,0.8427))
  cat(sprintf("Intercept:    %f (+ %f - %f)\n", coefs[2,1], coefs[3,1]-coefs[2,1], coefs[2,1]-coefs[1,1]))
  cat(sprintf("alpha(csw):   %f (+ %f - %f)\n", coefs[2,2], coefs[3,2]-coefs[2,2], coefs[2,2]-coefs[1,2]))
  cat(sprintf("b(beta):      %f (+ %f - %f)\n", coefs[2,3], coefs[3,3]-coefs[2,3], coefs[2,3]-coefs[1,3]))
  print(sqrt(cov(sim.coefs)))
  
  # extract the median coefficients only 
  coefs <- coefs[2,]

  save(amodel,file="a_csw_beta.model.Rdata")
  save(sim.amodel,file="a_csw_beta.sim.model.Rdata")

  # one symbol per csw value
  csw <- unique(adat$csw)
  syms <- c()
  sym <- 0
  for(i in csw){
    syms <- c(syms,rep(sym, times=length(which(adat$csw==i))) )
    sym <- sym + 1
  }
  
  # one colour per mass value
  mu <- unique(adat$mu)
  clr <- c()
  clr.idx <- 1
  clrs <- brewer.pal(n=9,name="Set1")
  for(i in mu){
    clr <- c(clr,rep(clrs[clr.idx], times=length(which(adat$mu==i))) )
    clr.idx <- clr.idx+1
    # we need to wrap around if there are too many
    if(clr.idx==(length(clrs)+1) ) clr.idx <- 1
  }

  tikzfiles <- tikz.init("csw_beta_a",width=3.5,height=3.5,lwdUnit=0.7)
  
   
  # the uncertainty in x comes from the slopes only
  # but we don't take it into account here
  # although this is not strictly true for the band shown here
  # it does mean however that the points have an uncertain horizontal position
  #coef.deriv <- matrix(nrow=2,ncol=n.predict)
  #coef.deriv[1,] <- newdata$csw
  #coef.deriv[2,] <- newdata$beta
  #xvar <- t(coef.deriv) %*% coefs.cov[2:3,2:3] %*% coef.deriv
  #dx <- sqrt(diag(xvar))

  poly.x <- -(coefs[2]*newdata$csw+coefs[3]*newdata$beta)
  
  poly.x <- c(poly.x,rev(poly.x))

  # the uncertainty in y we can take from the simulated confidence interval
  # plus the confidence interval of the model (which includes factors such as the removal of points)
  
  dyp <- sqrt( (a.confidence[,3]-a.confidence[,2])^2 + 
               (sim.confidence[,3]-sim.confidence[,2])^2 )
  dym <- sqrt( (a.confidence[,2]-a.confidence[,1])^2 + 
               (sim.confidence[,2]-sim.confidence[,1])^2 ) 
   
  poly.y <- c(a.confidence[,2]+dyp,rev(a.confidence[,2]-dym))

    
  # and finally account for the uncertainty of the x position of the
  # data points (since it depends on the model)
  coef.deriv <- matrix(nrow=2,ncol=length(adat$csw))
  coef.deriv[1,] <- adat$csw
  coef.deriv[2,] <- adat$beta
  xvar <- t(coef.deriv) %*% coefs.cov[2:3,2:3] %*% coef.deriv
  dx <- sqrt(diag(xvar))
  xpts <- -(coefs[2]*adat$csw+coefs[3]*adat$beta)

  # compute plot limits based on the supplied beta and csw ranges
  xlims <- -c(coefs[2]*cswrange[1]+coefs[3]*betarange[1], coefs[2]*cswrange[2]+coefs[3]*betarange[2] )

  # prepare plot area
  plot(x=xpts,y=adat$a,#log='y',
       xlim=xlims, ylim=arange,
       xlab=sprintf("$ %.3f~c_\\mathrm{sw} + %.3f~\\beta $",-coefs[2],-coefs[3]),
       ylab="$ a~[\\mathrm{fm}] $", 
       pch=syms, col=clr, type='n', yaxt='n' )

  axis(side=2, at=seq(from=arange[1],to=arange[2],by=0.01), las=2) 
   
  # draw error band and fit line
  polygon(x=poly.x,y=exp(poly.y),col="grey65", border=NA)
  lines(x=-(coefs[2]*newdata$csw+coefs[3]*newdata$beta),y=exp(a.confidence[,2]),lty=3)

  # to compare, let's try a simple linear fit
  #lin.amodel <- lm(a ~ csw + beta, data=adat, weights=1/adat$da^2)
  #print(anova(lin.amodel))
  #lin.prediction <- predict(lin.amodel,interval="none",newdata=newdata)
  #lines(y=lin.prediction,x=-(coefs[2]*newdata$csw+coefs[3]*newdata$beta),col="darkgreen",lty=3,lwd=2)
  
  # add points on top
  plotwitherror(x=xpts,
                y=adat$a,dy=adat$da,
                pch=syms, col=clr, rep=TRUE )
  
  csw.legend <- sprintf("$ c_\\mathrm{sw} = %s $", csw)
  legend(x="bottomleft",legend=csw.legend,pch=unique(syms), bty='n')
  mu.legend <- sprintf("$ a\\mu = %s $", mu )
  legend(x="topright",legend=mu.legend,pch=15,col=unique(clr), bty='n', pt.cex=1.5)

  tikz.finalize(tikzfiles)

  if(!missing(input)){
    a.confidence <- predict(amodel,newdata=input,interval="confidence",level=0.6854)
    # bring into the same format as sim.confidence below
    a.confidence <- exp(a.confidence[,c(2,1,3)])
    # for the fits on the simulated data sets, do the prediction
    sim.prediction <- lapply(X=sim.amodel, FUN=function(x) { predict(x,newdata=input,interval="none") } )
    sim.prediction <- array(data=unlist(sim.prediction),dim=c(length(input$csw),n.sim))
    # 68.54 confidence band for log(kappac) 
    sim.confidence <- exp(t(apply(X=sim.prediction,MARGIN=1,FUN=quantile,probs=c(0.1573,0.5,0.8427))))
    dyp <- sqrt( (a.confidence[,3]-a.confidence[,2])^2 + 
                 (sim.confidence[,3]-sim.confidence[,2])^2 
               )
    dym <- sqrt( (a.confidence[,2]-a.confidence[,1])^2 + 
                 (sim.confidence[,2]-sim.confidence[,1])^2 
               )

    for( i in 1:length(input$beta) ){
      cat(sprintf("a(beta=%f,csw=%f): %f (+ %f - %f) fm\n", input$beta[i], input$csw[i],
                                            a.confidence[i,2], dyp[i], dym[i] ) 
         )
    }
  }
}
  
