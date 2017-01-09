# this function performs a study of all possible fit ranges starting from timeslice 3
# for connected meson correlation functions, it systematically repeats both the direct
# fit to the functional form of the correlator dominated by the ground state
# as well as a constant fit to the effective mass

# it takes as input: hadron cf and effectivemass objects, the kappa parameter of the simulation as well
# as a pair of quark masses (a*mu) which are required for extracting the decay constant

meson_2pt_study_fitrange <- function(cf,effmass,name,kappa,m.sea,q_masses,Tmax,Tmin=5,minrange=5,fps.disprel='continuum',useCov=FALSE,debug=FALSE,boot.fit=FALSE) {
  if(!any(class(cf) == "cf")) {
    stop("study_fitrange requires that 'cf' is of class 'cf'!\n")
  }
  if(!any(class(effmass) == "effectivemass")) {
    stop("study_fitrange requires that 'effmass' is of class 'effectivemass'!\n")
  }
  if(!cf$boot.samples) {
    stop("study_fitrange requires cf to be bootstrapped!\n")
  }

  if(debug) {
    cat("Performing study of fit range dependence\n")
    print(q_masses)
  }
  require("plotrix")
  require("tikzDevice")

  # construct list of fitranges to be studied with placeholders for analysis results
  if(missing(Tmax)){
    Tmax <- (cf$Time/2)-1
  }
  fr <- list()
  for( T1 in Tmin:(Tmax-minrange) ) {
    for( T2 in (T1+minrange):Tmax ) {
      fr[[length(fr)+1]] <- list( t1=T1, t2=T2, boot.fit=boot.fit, m.sea=m.sea, q_masses=q_masses, name=name, Time=cf$Time, kappa=kappa,
                                   useCov=useCov, fps.disprel=fps.disprel,
                                    M=NA,P1=NA,P2=NA,f=NA,Meff=NA
                                  )
    }
  }
  class(fr) <- c(class(fr),"fitrange")
  
  require("parallel")
  ptm <- proc.time()
  rval <- mclapply(X=fr,FUN=do_fit.fitrange,cf=cf,effmass=effmass)
  if(debug){
    cat("Time for fit-range analysis\n")
    print(proc.time()-ptm)
  }

  error <- lapply(X=rval,FUN=function(x) { any(class(x)=="try-error") })
  if(any(unlist(error)))
    stop("try-errors detected in parallel execution of do_fit.fitrange")

  class(rval) <- c(class(rval),"fitrange")
  savename <- sprintf("%s.fitrange",gsub("-","_",name))
  assign(savename,rval)
  save(list=savename,file=sprintf("%s.Rdata",savename),compress=FALSE)
  
#  ptm <- proc.time()
#  rvaldf <- list(M=NULL,P1=NULL ,P2=NULL,f=NULL,Meff=NULL)
#  for( obj in rval ){
#    if(!is.na(obj$M)){
#        rvaldf$M <- cbind(rvaldf$M,obj$M$t)
#    }
#    if(!is.na(obj$f)){
#        rvaldf$f <- cbind(rvaldf$f,obj$f$t)
#    }
#    if(!is.na(obj$P1)){
#        rvaldf$P1 <- cbind(rvaldf$P1,obj$P1$t)
#    }
#    if(!is.na(obj$P2)){
#        rvaldf$P2 <- cbind(rvaldf$P2,obj$P2$t)
#    }
#    if(!is.na(obj$Meff)){
#        rvaldf$Meff <- cbind(rvaldf$Meff,obj$Meff$t)
#    }
#  }
#  for(df in rvaldf){
#    names(df) <- NULL
#  }
#  if(debug){
#    cat("Time for data reformatting\n")
#    print(proc.time()-ptm)
#  }
#
#  class(rvaldf) <- c(class(rvaldf),"fitrangedf")
#  savename <- sprintf("%s.fitrangedf",gsub("-","_",name))
#  assign(savename,rvaldf)
#  save(list=savename,file=sprintf("%s.Rdata",savename),compress=FALSE)
  
  rval
}

do_fit.fitrange <- function(fr,cf,effmass) {
  cat(fr$t1, "to", fr$t2, "\n")
  cf.matrixfit <- try(matrixfit(cf=cf, t1=fr$t1, t2=fr$t2, parlist=array(c(1,1,1,2,2,2), dim=c(2,3)),
                                matrix.size=3, symmetrise=T, boot.R=cf$boot.R, boot.l=cf$boot.l, 
                                useCov=fr$useCov, boot.fit=fr$boot.fit))
 
  if(!any(class(cf.matrixfit)=="try-error")) {
    cf.matrixfit <- computefps( cf.matrixfit, mu1=fr$q_masses$m1, mu2=fr$q_masses$m2, Kappa=fr$kappa, #normalisation="mpsexplicit",
                                disprel=fr$fps.disprel, boot.fit=fr$boot.fit )
  }

  effectivemass.fit <- try(fit.effectivemass(effmass, t1=fr$t1, t2=fr$t2, useCov=fr$useCov, replace.na=TRUE, boot.fit=fr$boot.fit))
  
  # if there were no errors during the fit procedure, propagate the fit results to the return value 
  if(fr$boot.fit){
    if(!any(class(cf.matrixfit)=="try-error")){
      fr$M <- list(t=cf.matrixfit$opt.tsboot[1,],t0=cf.matrixfit$opt.res$par[1],
                   Q=cf.matrixfit$Qval,chisq=cf.matrixfit$opt.res$value,dof=cf.matrixfit$dof )
      fr$P1 <- list(t=cf.matrixfit$opt.tsboot[2,],t0=cf.matrixfit$opt.res$par[2],
                    Q=cf.matrixfit$Qval,chisq=cf.matrixfit$opt.res$value,dof=cf.matrixfit$dof )
      fr$P2 <- list(t=cf.matrixfit$opt.tsboot[3,],t0=cf.matrixfit$opt.res$par[3],
                    Q=cf.matrixfit$Qval,chisq=cf.matrixfit$opt.res$value,dof=cf.matrixfit$dof )
      fr$f <- list(t=cf.matrixfit$fps.tsboot,t0=cf$matrixfit$fps,
                   Q=cf.matrixfit$Qval,chisq=cf.matrixfit$opt.res$value,dof=cf.matrixfit$dof )
    }
    if(!any(class(effectivemass.fit)=="try-error")){
      fr$Meff <- list(t=effectivemass.fit$massfit.tsboot[,1],t0=effectivemass.fit$opt.res$par[1],
                      Q=effectivemass.fit$Qval,chisq=effectivemass.fit$opt.res$value,dof=effectivemass.fit$dof )
    }
  }else{
    if(!any(class(cf.matrixfit)=="try-error")){
      fr$M <- list(t=NA,t0=cf.matrixfit$opt.res$par[1],
                   Q=cf.matrixfit$Qval,chisq=cf.matrixfit$opt.res$value,dof=cf.matrixfit$dof )
      fr$P1 <- list(t=NA,t0=cf.matrixfit$opt.res$par[2],
                    Q=cf.matrixfit$Qval,chisq=cf.matrixfit$opt.res$value,dof=cf.matrixfit$dof )
      fr$P2 <- list(t=NA,t0=cf.matrixfit$opt.res$par[3],
                    Q=cf.matrixfit$Qval,chisq=cf.matrixfit$opt.res$value,dof=cf.matrixfit$dof )
      fr$f <- list(t=NA,t0=cf$matrixfit$fps,
                   Q=cf.matrixfit$Qval,chisq=cf.matrixfit$opt.res$value,dof=cf.matrixfit$dof )
    }
    if(!any(class(effectivemass.fit)=="try-error")){
      fr$Meff <- list(t=NA,t0=effectivemass.fit$opt.res$par[1],
                      Q=effectivemass.fit$Qval,chisq=effectivemass.fit$opt.res$value,dof=effectivemass.fit$dof )
    }
  }
  return(invisible(fr))
}

plot.fitrange <- function(fr) {
  require("plotrix")
  require("tikzDevice")
  if(!any(class(fr)=="fitrange")) {
    stop("object for plotting must be of class 'fitrange'")
  }

#  # we now remove the outliers using the usual idea of computing quartiles and the interquartile range
#  quants <- quantile(fr$M)
#  tshld.hi <- quants[4] + 1.5*IQR(fr$M)
#  tshld.lo <- quants[2] - 1.5*IQR(fr$M)
#
#  outlier.indices <- which( fr$M < tshld.lo | fr$M > tshld.hi )
#
#  if( length(fr[,1]) == length(outlier.indices) ) {
#    warning("For ", fr$name[1], " no entries remain after removal of outliers! Continuing with full set!")
#  } else if( length(outlier.indices) == 0 ) {
#    warning("Interquartile range very large, no outliers found!\n")
#  } else {
#    fr <- fr[-outlier.indices,]
#  }

  # assemble relevant data for convenient plotting

  # first compute the weights that are going to be used for the weighted histogram
  # and the final determination of the error
  wM <- wP1 <- wP2 <- (1-2*abs(fr$matrixfit.Q-0.5))^2
  wMeff <- (1-2*abs(fr$effectivemass.Q-0.5))^2
  if(fr$boot.fit){
    wM <- wM*(min(fr$dM)/fr$dM)^2
    wP1 <- wP1*(min(fr$dP1)/fr$dP1)^2
    wP2 <- wP2*(min(fr$dP2)/fr$dP2)^2
    wMeff <- wMeff*(min(fr$dMeff)/fr$dMeff)^2 
  }

  l.matrixM <- list(df=data.frame( val=fr$M, dval=fr$dM, t1=fr$t1, t2=fr$t2,
                               ChiSqr.ov.dof=(fr$matrixfit.ChiSqr/fr$matrixfit.dof),
                               Q=fr$matrixfit.Q, w=wM),
                               label="$M_{\\mathrm{mtx}}$",
                               name="matrixfit M" )

  l.matrixP1 <- list(df=data.frame( val=fr$P1, dval=fr$dP1, t1=fr$t1, t2=fr$t2,
                               ChiSqr.ov.dof=(fr$matrixfit.ChiSqr/fr$matrixfit.dof),
                               Q=fr$matrixfit.Q, w=wP1),
                               label="$P1_{\\mathrm{mtx}}$",
                               name="matrixfit P1" )
  
  l.matrixP2 <- list(df=data.frame( val=fr$P2, dval=fr$dP2, t1=fr$t1, t2=fr$t2,
                               ChiSqr.ov.dof=(fr$matrixfit.ChiSqr/fr$matrixfit.dof),
                               Q=fr$matrixfit.Q, w=wM),
                               label="$P2_{\\mathrm{mtx}}$",
                               name="matrixfit P2" )

  l.effM <- list(df=data.frame( val=fr$Meff, dval=fr$dMeff,t1=fr$t1, t2=fr$t2,
                                ChiSqr.ov.dof=(fr$effectivemass.ChiSqr/fr$effectivemass.dof),
                                Q=fr$effectivemass.Q, w=wMeff),
                                label="$M_{\\mathrm{eff}}$",
                                name="eff M" )

  # produce a number of plots relating to the fit range analysis
  temp <- sprintf("%s.fitrange.%s",fr$name[1],c("tex","pdf","aux","log"))
  tikzfiles <- list(tex=temp[1],pdf=temp[2],aux=temp[3],log=temp[4])
  rm(temp)
  tikz(tikzfiles$tex, standAlone = TRUE, width=5, height=5)

  # colours to add some timeslice information
  colours <- rainbow(n=fr$Time[1]/2,start=0.5)
  # colours to provide information about Q
  Qcolours <- rainbow(n=100,start=0.5)
  
  for( l in list( l.matrixM, l.matrixP1, l.matrixP2, l.effM ) ) {
    df <- l$df
    qtyname <- l$name
    label <- l$label
    title <- sprintf("%s %s",qtyname,fr$name[1])
    # prefix all underscore with a backslash so that latex understands properly
    # note that we need to escape the backslash twice in order to have two backslashes
    # in the resulting string
    title <- gsub("_","\\\\_",title)
    print(title)

    hist(df$val,breaks=40,main=title,xlab=label)

    #plot(density(df$val),main=title, xlab=label)
    #abline(v=mean(df$val),col='red')
    #abline(v=median(df$val),col='blue')

    plot(y=df$ChiSqr.ov.dof,x=df$val,main=title,
         col=colours[df$t2],ylab="$\\chi^2 / \\mathrm{d.o.f}$",xlab=label
        )

    plot(y=df$w,x=df$val,main=title,col=colours[df$t2],ylab='$w$',xlab=label)
    
    weighted.hist(x=df$val,w=df$w,breaks=40,main=paste("weighted",title),xlab=label)

    # save some more lines by doing two sets of plots in one go
    # plot the weight and quantity as a function of several combinations of the fit range
    for( dat in list( list(qty=df$w,lab='w'), list(qty=df$val,lab=label) ) ) {
      #plot(x=df$t2-df$t1,y=dat$qty,main=title,col=colours[df$t2],xlab="$t_f - t_i$",ylab=dat$lab)
      #plot(x=df$t1+df$t2,y=dat$qty,main=title,col=colours[df$t2],xlab="$t_i + t_f$",ylab=dat$lab)
      #plot(x=df$t1,y=dat$qty,main=title,col=colours[df$t2],xlab="$t_i$",ylab=dat$lab)
      #plot(x=df$t2,y=dat$qty,main=title,col=colours[df$t1],xlab="$t_f$",ylab=dat$lab)
      weighted.hist(x=df$t1,w=df$w,breaks=0:(fr$Time[1]/2),main=paste("weighted",title),xlab="$t_i$")
      weighted.hist(x=df$t2,w=df$w,breaks=0:(fr$Time[1]/2),main=paste("weighted",title),xlab="$t_f$")
      weighted.hist(x=df$t2-df$t1,w=df$w,breaks=0:fr$Time[1],main=paste("weighted",title),xlab="$t_f-t_i$")
      weighted.hist(x=df$t2+df$t1,w=df$w,breaks=0:fr$Time[1],main=paste("weighted",title),xlab="$t_f+t_i$")
    }

    boxplot(df$val,main=title,ylab=label)
  }

  # weighted mean and variance for the matrixfit masses
  matrixfit.mustar <- weighted.mean(l.matrixM$df$val,l.matrixM$df$w)
  matrixfit.varstar <- weighted.variance(l.matrixM$df$val,l.matrixM$df$w)

  # weighted mean and variance for the effective masses
  effective.mustar <- weighted.mean(l.effM$df$val,l.effM$df$w)
  effective.varstar <- weighted.variance(l.effM$df$val,l.effM$df$w)

  cat("Effective mass:", effective.mustar, sqrt(effective.varstar), "\n")
  cat("Matrixfit:", matrixfit.mustar, sqrt(matrixfit.varstar), "\n")

#  effseq <- seq(from=min(l.effmass$df$val),to=max(l.effmass$df$val),length.out=500)
#  effnorm <- dnorm(x=effseq,mean=effective.mustar,sd=sqrt(effective.varstar))
#  matrixseq <- seq(from=min(l.matrix$df$val),to=max(l.matrix$df$val),length.out=500)
#  matrixnorm <- dnorm(x=matrixseq,mean=matrixfit.mustar,sd=sqrt(matrixfit.varstar))
#
#  ylims <- c(0,max( c(effnorm,matrixnorm) ) )
#
#  plot(y=effnorm, x=effseq,lwd=3,type='l',col='blue',
#      xlim=effective.mustar+4*c(-1,1)*sqrt(effective.varstar),ylim=ylims,xlab="mass",ylab="",
#      main=sprintf("weighted gaussians %s", gsub("_","\\\\_",name)))
#
#  lines(y=matrixnorm, x=matrixseq,lwd=3,col='red',type='l')


  # plot effective mass vs. matrixfit mass
  plot(y=l.effM$df$val,x=l.matrixM$df$val, main="effective vs. matrixfit",
       ylab="effective mass", xlab="matrixfit mass",col=colours[l.effM$df$t2])

  # plot chosen fit ranges, indicating the weight by the colour
  plot(x=l.matrixM$df$t1,y=l.matrixM$df$t2,main=sprintf("chosen fitranges %s", gsub("_","\\\\_",fr$name[1])),
       xlab=expression(t[i]),ylab=expression(t[f]),col=Qcolours[as.integer(100*l.matrixM$df$w)])

  dev.off()
  tools::texi2dvi(tikzfiles$tex,pdf=T)                                                                                                                                                                         
  # use pdfcrop tool to remove plot borders
  command <- sprintf("pdfcrop %s %s",tikzfiles$pdf,tikzfiles$pdf)
  system(command)
  # remove temporary files 
  command <- sprintf("rm %s %s %s", tikzfiles$tex, tikzfiles$log, tikzfiles$aux)
  system(command)
  
  return()
}

summary.fitrange <- function(fr) {
  if(!any(class(fr)=="fitrange")) {
    stop("object for plotting must be of class 'fitrange'")
  }
  
#  t1=T1, t2=T2, boot.fit=boot.fit, q_masses=q_masses, name=name, Time=cf$Time, kappa=kappa,

#  sprintf("%4s %3s %3s %5s %5s %7s
#   lapply(X=fitrange,FUN=function(x) { sprintf("t1: 
}
