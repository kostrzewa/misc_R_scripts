## TODO: propagate some relevant paramters or extract from cf/effmass

meson_2pt_study_fitrange <- function(cf,effmass,name,debug=FALSE) {
  if(!any(class(cf) == "cf")) {
    stop("study_fitrange requires that 'cf' is of class 'cf'!\n")
  }
  if(!any(class(effmass) == "effectivemass")) {
    stop("study_fitrange requires that 'effmass' is of class 'effectivemass'!\n")
  }
  
  if(debug) {
    cat("Performing study of fit range dependence\n")
  }
  require("plotrix")

  res <- data.frame(name=c(), t1=c(), t2=c(), M=c(), dM=c(), 
                Meff=c(), dMeff=c(), P1=c(), dP1=c(), 
                P2=c(), dP2=c(), f=c(), df=c(), 
                qm1=c(), qm2=c(), kappa=c(),matrixfit.ChiSqr=c(),
                matrixfit.dChiSqr=c(), matrixfit.dof=c(), effectivemass.ChiSqr=c(), 
                effectivemass.dChiSqr=c(), effectivemass.dof=c() )

  for( T1 in 3:((cf$Time/2)-1-2) ) {
    for( T2 in (T1+2):((cf$Time/2)-1) ) {
      cat(T1, "to", T2, "\n")
      
      cf.matrixfit <- try(matrixfit(cf=cf, t1=T1, t2=T2, parlist=array(c(1,1,1,2,2,2), dim=c(2,3)),
                                matrix.size=3, symmetrise=T, boot.R=boot.R, boot.l=boot.l, useCov=useCov, boot.fit=boot.fit))
      
      # if the fit fails, we go to the next fit range
      if( any(class(cf.matrixfit) == "try-error" )) next
      
      cf.matrixfit <- computefps( cf.matrixfit, mu1=q_masses$m1, mu2=q_masses$m2, Kappa=kappa, disprel='continuum', boot.fit=boot.fit )
  
      effectivemass.fit <- try(fit.effectivemass(effmass, t1=T1, t2=T2, useCov=useCov, replace.na=TRUE, boot.fit=boot.fit))
      if( any(class(effectivemass.fit) == "try-error" )) next

      # if we skip bootstrapping the two fits (boot.fit==FALSE), then we can't compute errors
      dM <- NA
      dMeff <- NA
      dP1 <- NA
      dP2 <- NA
      df <- NA
      matrixfit.dChiSqr <- NA
      effectivemass.dChiSqr <- NA

      if(boot.fit) {
        dM <- sd(cf.matrixfit$opt.tsboot[1,])
        dMeff <- sd(effectivemass.fit$massfit.tsboot[,1])
        dP1 <- sd(cf.matrixfit$opt.tsboot[2,])
        dP2 <- sd(cf.matrixfit$opt.tsboot[3,])
        df <- sd(cf.matrixfit$fps.tsboot)
        matrixfit.dChiSqr <- sd(cf.matrixfit$opt.tsboot[,4])
        effectivemass.dChiSqr <- sd(effectivemass.fit$massfit.tsboot[,2])
      }

      res <- rbind(res,data.frame(name=name, t1=T1, t2=T2, M=cf.matrixfit$opt.res$par[1], dM=dM, 
                Meff=effectivemass.fit$opt.res$par[1], dMeff=dMeff, 
                P1=cf.matrixfit$opt.res$par[2], dP1=dP1, 
                P2=cf.matrixfit$opt.res$par[3], dP2=dP2, 
                f=cf.matrixfit$fps, df=df, 
                qm1=q_masses$m1, qm2=q_masses$m2, kappa=kappa,
                matrixfit.ChiSqr=cf.matrixfit$opt.res$value, matrixfit.dChiSqr=matrixfit.dChiSqr,
                matrixfit.dof=cf.matrixfit$dof,
                effectivemass.ChiSqr=effectivemass.fit$opt.res$value, effectivemass.dChiSqr=effectivemass.dChiSqr,
                effectivemass.dof=effectivemass.fit$dof ))
    }
  }

  # need to rename object before saving to file
  savename <- sprintf("%s.fitrange",name)
  assign(savename,res)
  save(list=savename,file=sprintf("%s.Rdata",savename)) 

  # we now remove the outliers using the usual idea of computing quartiles and the interquartile range

  quants <- quantile(res$M)
  tshld.hi <- quants[4] + 1.5*IQR(res$M)
  tshld.lo <- quants[2] - 1.5*IQR(res$M)

  outlier.indices <- which( res$M < tshld.lo | res$M > tshld.hi )

  res <- res[-outlier.indices,]
  
  # assemble relevant data for convenient plotting
  plotdf.matrix <- data.frame( name="matrixfit", label=expression(M[mtx]),
                               val=res$M, t1=res$t1, t2=res$t2, 
                               ChiSqr.ov.dof=(res$matrixfit.ChiSqr/res$matrixfit.dof), 
                               Q=(1-pchisq(df=res$matrixfit.dof,q=res$matrixfit.ChiSqr)) )
  
  plotdf.effmass <- data.frame( name="effective", label=expression(M[eff]),
                                val=res$Meff, t1=res$t1, t2=res$t2, 
                                ChiSqr.ov.dof=(res$effectivemass.ChiSqr/res$effectivemass.dof), 
                                Q=(1-pchisq(df=res$effectivemass.dof,q=res$effectivemass.ChiSqr)) )

  # produce a number of plots relating to the fit range analysis 
  pdf(sprintf("%s.fitrange.pdf",directory),title=directory)
  
  # colours to add some timeslice information
  colours <- rainbow(n=cf$Time/2)
  for( df in list( plotdf.matrix, plotdf.effmass ) ) {
    title <- paste(df$name,name)
    
    hist(df$val,breaks=40,main=title,xlab=df$label)
    
    plot(density(df$val),main=title, xlab=df$label)
    abline(v=mean(df$val),col='red')
    abline(v=median(df$val),col='blue')
    
    plot(y=df$ChiSqr.ov.dof,x=df$val,main=title,
         col=colours[df$t2],ylab=expression(Chi^2 %/% d.o.f),xlab=df$label
        )
    
    plot(y=df$Q,x=df$val,main=title,col=colours[df$t2],ylab='Q',xlab=df$label)
    
    weighted.hist(x=df$val,w=df$Q,breaks=40,main=paste("weighted",title),xlab=df$label)
    
    # save some more lines by doing two sets of plots in one go
    for( dat in list( data.frame(qty=df$Q,lab='Q'), data.frame(qty=df$val,lab=df$label) ) ) {
      plot(x=df$t2-df$t1,y=dat$qty,main=title,col=colours[df$t2],xlab=expression(t[f] %-% t[i]),ylab=dat$lab)
      plot(x=df$t1+df$t2,y=df$Q,main=title,col=colours[df$t2],xlab=expression(t[i] %+% t[f]),ylab=dat$lab)
      plot(x=df$t1,y=dat$qty,main=title,col=colours[df$t2],xlab=expression(t[i]),ylab=dat$lab)
      plot(x=df$t2,y=dat$qty,main=title,col=colours[df$t1],xlab=expression(t[f]),ylab=dat$lab)
    }
    
    boxplot(df$val,main=title,ylab=df$label)
  }
  
  # weighted mean and variance for the matrixfit masses
  matrixfit.mustar <- weighted.mean(plotdf.matrix$val,plotdf.matrix$Q)
  matrixfit.varstar <- weighted.variance(plotdf.matrix$val,plotdf.matrix$Q)
  
  # weighted mean and variance for the effective masses
  effective.mustar <- weighted.mean(plotdf.effective$val,plotdf.effective$Q)
  effective.varstar <- weighted.variance(effective.masses,effective.Q)
  
  cat("Effective mass:", effective.mustar, sqrt(effective.varstar), "\n")
  cat("Matrixfit:", matrixfit.mustar, sqrt(matrixfit.varstar), "\n")  
  
  effseq <- seq(from=min(effective.masses),to=max(effective.masses),length.out=500)
  effnorm <- dnorm(x=effseq,mean=effective.mustar,sd=sqrt(effective.varstar))
  matrixseq <- seq(from=min(matrixfit.masses),to=max(matrixfit.masses),length.out=500)
  matrixnorm <- dnorm(x=matrixseq,mean=matrixfit.mustar,sd=sqrt(matrixfit.varstar))
  
  ylims <- c(0,max( c(effnorm,matrixnorm) ) )

  plot(y=effnorm, x=effseq,lwd=3,type='l',col='blue',
      xlim=effective.mustar+4*c(-1,1)*sqrt(effective.varstar),ylim=ylims,xlab="mass",ylab="",
      main=sprintf("weighted gaussians %s",name))
  
  lines(y=matrixnorm, x=matrixseq,lwd=3,col='red',type='l')

  plot(y=plotdf.effmass$val,x=plotdf.matrix$val, main="effective vs. matrixfit", 
       ylab="effective mass", xlab="matrixfit mass",col=colours[plotdf.effmass$t2])  

  plot(x=plotdf.matrix$t1,y=plotdf.matrix$t2,main=sprintf("chosen fitranges %s", name),
       xlab=expression(t[i]),ylab=expression(t[f]),col=colurs[plotdf.matrix$t2])

  dev.off()
  
  return(invisible(res))
}