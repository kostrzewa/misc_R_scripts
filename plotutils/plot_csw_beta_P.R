plot_csw_beta_P <- function(datafile,
                                 cswrange=c(1.0,2.5),betarange=c(1.6,1.8),Prange=c(0.45,0.65),
                                 n.predict=1000, n.sim=200, a=0.080, a.col=rgb(0.0,1.0,0.0,0.6) ) {

  require("RColorBrewer")

  Pdat <- read.table(file=datafile,header=TRUE,stringsAsFactors=FALSE, fill=TRUE)
  Pdat <- Pdat[order(Pdat$csw,Pdat$beta),]
  Pdat <- cbind(Pdat,logP=log(Pdat$P))
  Pdat <- cbind(Pdat,dlogP=Pdat$dP/Pdat$P)
  save(Pdat,file="csw_beta_P.Rdata")

  # for the fit we take into account the number of trajectories used for the estimate of w0
  # the inverse error on a and the inverse absolute value of ampcac, so that 
  # critical ensembles are assigned a larger weight
  #weights <- Pdat$Nt/abs(Pdat$mpcac)/Pdat$dlogP^2
  weights <- (max(abs(Pdat$mpcac))/abs(Pdat$mpcac))/Pdat$dP^2

  # prepare some simulated data with gaussian distributions
  simdata <- data.frame(matrix(ncol=length(Pdat$P),nrow=n.sim))
  for( k in 1:length(Pdat$P) ) {
    simdata[,k] <- rnorm(n=n.sim,mean=Pdat$P[k],sd=Pdat$dP[k])
  }
  
  # csw,beta pairs for the prediction 
  newdata <- data.frame(csw=seq(cswrange[1]-0.5,cswrange[2]+0.5,length.out=n.predict),
                        beta=seq(betarange[1]-0.5,betarange[2]+0.5,length.out=n.predict))

  # model for log(a)
  Pmodel <- lm(P ~ csw + beta, data=Pdat, weights=weights, trace=TRUE)
  #print(anova(Pmodel))
  
  # and confidence intervals which are calculated by removing points
  P.confidence <- predict(Pmodel,newdata=newdata,interval="confidence",level=0.6854)
  # bring into the same format as sim.confidence below
  P.confidence <- P.confidence[,c(2,1,3)]

  # and repeated on all the simulated data
  sim.Pmodel <- apply(X=simdata,MARGIN=1,
                        FUN=function(x) {
                              lm(P ~ csw + beta, 
                                 data=data.frame(P=x,csw=Pdat$csw,beta=Pdat$beta), 
                                 weights=weights)
                            } )
  
  # for the fits on the simulated data sets, do the prediction
  sim.prediction <- lapply(X=sim.Pmodel, FUN=function(x) { predict(x,newdata=newdata,interval="none") } )
  sim.prediction <- array(data=unlist(sim.prediction),dim=c(n.predict,n.sim))
  
  # 68.57 confidence band for log(a) 
  sim.confidence <- t(apply(X=sim.prediction,MARGIN=1,FUN=quantile,probs=c(0.1573,0.5,0.8427)))
  
  # obtain some confidence levels for the fit coefficients 
  sim.coefs <- lapply(X=sim.Pmodel, FUN=function(x) { x$coefficients })
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

  save(Pmodel,file="P_csw_beta.model.Rdata")
  P.cfs <- Pmodel$coefficients[2:3]

  # one symbol per csw value
  csw <- unique(Pdat$csw)
  syms <- c()
  sym.cur <- 0
  for(i in csw){
    syms <- c(syms,rep(sym.cur, times=length(which(Pdat$csw==i))) )
    sym.cur <- sym.cur + 1
  }
  
  # one colour per mass value
  mu <- unique(Pdat$mu)
  clrs <- c()
  clr.idx <- 1
  clr.pal <- brewer.pal(n=9,name="Set1")
  for(i in mu){
    clrs <- c(clrs,rep(clr.pal[clr.idx], times=length(which(Pdat$mu==i))) )
    clr.idx <- clr.idx+1
    # we need to wrap around if there are too many
    if(clr.idx==(length(clr.pal)+1) ) clr.idx <- 1
  }

  tikzfiles <- tikz.init("csw_beta_P",width=3.5,height=3.5,lwdUnit=0.7)
   
  # the uncertainty in x comes from the slopes only
  # but we don't take it into account here
  coef.deriv <- matrix(nrow=2,ncol=n.predict)
  coef.deriv[1,] <- newdata$csw
  coef.deriv[2,] <- newdata$beta
  xvar <- t(coef.deriv) %*% coefs.cov[2:3,2:3] %*% coef.deriv
  dx <- sqrt(diag(xvar))

  poly.x <- (coefs[2]*newdata$csw+coefs[3]*newdata$beta)
  
  poly.x <- c(poly.x-dx,rev(poly.x)+dx)

  # the uncertainty in y we can take from the simulated confidence interval
  # plus the confidence interval of the model (which includes factors such as the removal of points)
  
  dyp <- sqrt( (P.confidence[,3]-P.confidence[,2])^2 + 
               (sim.confidence[,3]-sim.confidence[,2])^2 )
  dym <- sqrt( (P.confidence[,2]-P.confidence[,1])^2 + 
               (sim.confidence[,2]-sim.confidence[,1])^2 ) 
   
  poly.y <- c(P.confidence[,2]+dyp,rev(P.confidence[,2]-dym))

    
  # and finally account for the uncertainty of the x position of the
  # data points (since it depends on the model)
  coef.deriv <- matrix(nrow=2,ncol=length(Pdat$csw))
  coef.deriv[1,] <- Pdat$csw
  coef.deriv[2,] <- Pdat$beta
  xvar <- t(coef.deriv) %*% coefs.cov[2:3,2:3] %*% coef.deriv
  dx <- sqrt(diag(xvar))
  xpts <- (coefs[2]*Pdat$csw+coefs[3]*Pdat$beta)

  # compute plot limits based on the supplied beta and csw ranges
  xlims <- c(coefs[2]*cswrange[1]+coefs[3]*betarange[1], coefs[2]*cswrange[2]+coefs[3]*betarange[2] )

  # prepare plot area
  plot(x=xpts,y=Pdat$P,#log='y',
       xlim=xlims, ylim=Prange,
       xlab=sprintf("$ %.4f~c_\\mathrm{sw} + %.4f~\\beta $",coefs[2],coefs[3]),
       ylab="$ \\left\\langle P \\right\\rangle $", 
       pch=syms, col=clrs, type='n', las=1)
   
  # error band and fit line
  polygon(x=poly.x,y=poly.y,col="grey65", border=NA)
  lines(x=coefs[2]*newdata$csw+coefs[3]*newdata$beta,y=P.confidence[,2],lty=3)

  # now we load amodel and invert it for a (fm)
  #load("a_csw_beta.model.Rdata")
  load("a_csw_beta.sim.model.Rdata")

  #a.cfs <- amodel$coefficients
  #print(summary(amodel))
  
  #test.csw <- (1/a.cfs[2])*(log(a)-a.cfs[1]-a.cfs[3]*newdata$beta)
  
  # now we load amodel and invert it for a (fm)
  match.csw <- lapply(X=sim.amodel,
                     FUN=function(x) {
                       cfs <- x$coefficients
                       (1/cfs[2])*(log(a)-cfs[1]-cfs[3]*newdata$beta)
                     }
                     )

  match.csw <- t(array(unlist(match.csw),dim=c(n.predict,n.sim)))
  match.csw <- apply(X=match.csw,MARGIN=2,FUN=quantile,probs=c(0.1573,0.5,0.8427)) 

  match.x <- coefs[2]*match.csw[2,]+coefs[3]*newdata$beta
  
  # the plaquette value for tadpole improved perturbation theory as a function of csw and beta
  match.P <- (0.113*6)/(match.csw[2,]*newdata$beta-newdata$beta)

  # this can become singular, so let's filter it so it fits into the Prange
  div.idx.match.P <- union(which(match.P > 1),which(match.P < 0))
  match.csw.restrict <- match.csw[,-div.idx.match.P]
  match.x <- match.x[-div.idx.match.P]
  match.P <- match.P[-div.idx.match.P]
  poly.x <- c(match.x,rev(match.x))

  # consider the uncertainty on 0.113(3) and the one on match.csw coming from the ln(a) fit
  # noting that when match.csw fluctuates up, match.P fluctuates down
  match.mdP <- sqrt( (match.P*0.003/0.113)^2 
                     + (match.P/(match.csw.restrict[2,]-1))^2*(match.csw.restrict[3,]-match.csw.restrict[2,])^2
                   )
  match.pdP <- sqrt( (match.P*0.003/0.113)^2 
                     + (match.P/(match.csw.restrict[2,]-1))^2*(match.csw.restrict[2,]-match.csw.restrict[1,])^2
                   )
  poly.y <- c(match.P+match.pdP,rev(match.P-match.mdP))
  #polygon(x=poly.x,y=poly.y,col=a.col, border=NA)
  #lines(x=match.x,y=match.P,lty=3,col="forestgreen")

  # add points on top
  plotwitherror(x=xpts,dx=dx,
                y=Pdat$P,dy=Pdat$dP,
                pch=syms, col=clrs, rep=TRUE )
  
  csw.legend <- sprintf("$ c_\\mathrm{sw} = %s $", csw)
  legend(x="bottomright",legend=csw.legend,pch=unique(syms), bty='n')

  mu.legend <- c(sprintf("$ a\\mu = %s $", mu ),"$ \\langle P \\rangle=\\frac{6\\cdot0.113}{\\beta (\\tilde c_\\mathrm{sw} - 1)} $",sprintf("$\\tilde c_\\mathrm{sw}(a=%s~\\mathrm{fm},\\beta)$",a))
  #legend(x="topleft",legend=mu.legend,pch=c(15,15,15,NA),col=c(unique(clrs),a.col,NA),lty=c(NA,NA,3,NA), 
  #       bty='n',pt.cex=1.5,cex=0.8)
  # removed the green line here
  legend(x="topleft",legend=mu.legend[1:2],pch=c(15,15),col=c(unique(clrs)), 
         bty='n',pt.cex=1.5,cex=0.8) 

  
  plot(x=newdata$beta,y=match.csw[2,],
       xlab="$\\beta$",ylab="$c_\\mathrm{sw}$",
       xlim=betarange,ylim=cswrange,las=1,type='n')

  poly.x <- c(newdata$beta,rev(newdata$beta))
  poly.y <- c(match.csw[3,],rev(match.csw[1,]))
  polygon(x=poly.x,y=poly.y,col=a.col,border=NA)
  lines(x=newdata$beta,y=match.csw[2,],col="forestgreen",lty=3)

  #
  P.csw <- lapply(
                  X=sim.Pmodel,
                  FUN=function(x){
                        coefs <- x$coefficients
                        q.a <- coefs[2]*newdata$beta
                        q.b <- (newdata$beta*coefs[1]+coefs[3]*newdata$beta^2-coefs[2]*newdata$beta)
                        q.c <- -(coefs[1]*newdata$beta+coefs[3]*newdata$beta^2+6*0.113)

                        Re( (-q.b+sqrt(q.b^2-4*q.a*q.c))/(2*q.a) )
                  }
        )
  #


  q.a <- coefs[2]*newdata$beta
  q.b <- (newdata$beta*coefs[1]+coefs[3]*newdata$beta^2-coefs[2]*newdata$beta)
  q.c <- -(coefs[1]*newdata$beta+coefs[3]*newdata$beta^2+6*0.113)

  P.dcsw <- 6/sqrt(q.b^2-4*q.a*q.c)

  P.csw <- t(array(unlist(P.csw),dim=c(length(newdata$beta),n.sim)))
  P.csw <- apply(X=P.csw,MARGIN=2,FUN=quantile,probs=c(0.1573,0.5,0.8427))
 
  P.csw.df <- data.frame(beta=newdata$beta, csw=P.csw[2,], match.csw=match.csw[2,])
  print(P.csw.df)

  
  P.pdcsw <- sqrt( (P.csw[3,]-P.csw[2,])^2
                  + P.dcsw^2*0.003^2 
                 ) 
  P.mdcsw <- sqrt( (P.csw[2,]-P.csw[1,])^2
                  + P.dcsw^2*0.003^2 
                 ) 
  #
  #print(data.frame(P.mdcsw,P.csw[2,],P.pdcsw))


  poly.y <- c(P.csw[2,]+P.pdcsw,rev(P.csw[2,]-P.mdcsw))

  polygon(x=poly.x,y=poly.y,col=rgb(0.0,0.0,0.0,0.3),border=NA)
  lines(y=P.csw[2,],x=newdata$beta)
  # indicate tree level value for c_sw
  abline(col="grey45",h=1.0,lty=2)


  lg <- c(sprintf("$c_\\mathrm{sw}(a=%s~\\mathrm{fm},\\beta)$",a),"$c_\\mathrm{sw}=1+0.113(3)\\frac{6}{\\beta \\langle P \\rangle}$")
  legend(x="topright",lg,pch=c(15,15),col=c(a.col,rgb(0.0,0.0,0.0,0.3)),lty=c(3,1), 
         bty='n',pt.cex=1.5,cex=0.8)
  tikz.finalize(tikzfiles) 
 
  # finally, plot csw(beta) for a larger range of values
  tikzfiles <- tikz.init("csw_beta_wide",width=4,height=4  ,lwdUnit=0.7)
  newdata <- data.frame(csw=seq(0.8,5,length.out=n.predict),
                        beta=seq(0.001,10,length.out=n.predict))
  plot(x=newdata$beta,y=match.csw[2,],
       xlab="$\\beta$",ylab="$c_\\mathrm{sw}$",
       xlim=c(1,3),ylim=c(0.8,3),las=1,type='n')

  #
  P.csw <- lapply(
                  X=sim.Pmodel,
                  FUN=function(x){
                        coefs <- x$coefficients
                        q.a <- coefs[2]*newdata$beta
                        q.b <- (newdata$beta*coefs[1]+coefs[3]*newdata$beta^2-coefs[2]*newdata$beta)
                        q.c <- -(coefs[1]*newdata$beta+coefs[3]*newdata$beta^2+6*0.113)

                        Re( (-q.b+sqrt(q.b^2-4*q.a*q.c))/(2*q.a) )
                  }
        )
  #
  q.a <- coefs[2]*newdata$beta
  q.b <- (newdata$beta*coefs[1]+coefs[3]*newdata$beta^2-coefs[2]*newdata$beta)
  q.c <- -(coefs[1]*newdata$beta+coefs[3]*newdata$beta^2+6*0.113)

  P.dcsw <- 6/sqrt(q.b^2-4*q.a*q.c)

  P.csw <- t(array(unlist(P.csw),dim=c(length(newdata$beta),n.sim)))
  P.csw <- apply(X=P.csw,MARGIN=2,FUN=quantile,probs=c(0.1573,0.5,0.8427))
 
  P.pdcsw <- sqrt( (P.csw[3,]-P.csw[2,])^2
                  + P.dcsw^2*0.003^2 
                 ) 
  P.mdcsw <- sqrt( (P.csw[2,]-P.csw[1,])^2
                  + P.dcsw^2*0.003^2 
                 ) 
  #
  #print(data.frame(P.mdcsw,P.csw[2,],P.pdcsw))

  poly.x <- c(newdata$beta,rev(newdata$beta))
  poly.y <- c(P.csw[2,]+P.pdcsw,rev(P.csw[2,]-P.mdcsw))
  poly.y <- c(P.csw[2,]+P.pdcsw,rev(P.csw[2,]-P.mdcsw))
  polygon(x=poly.x,y=poly.y,col=rgb(0.0,0.0,0.0,0.3),border=NA)
  lines(y=P.csw[2,],x=newdata$beta)
  # indicate tree level value for c_sw
  abline(col="grey45",h=1.0,lty=2)

  tikz.finalize(tikzfiles)

}
  
