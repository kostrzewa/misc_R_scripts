# this function performs a study of all possible fit ranges starting from timeslice 3
# for connected meson correlation functions, it systematically repeats both the direct
# fit to the functional form of the correlator dominated by the ground state
# as well as a constant fit to the effective mass

# it takes as input: hadron cf and effectivemass objects, the kappa parameter of the simulation as well
# as a pair of quark masses (a*mu) which are required for extracting the decay constant

meson_2pt_study_fitrange <- function(cf,effmass,name,kappa,q_masses,fps.disprel='continuum',useCov=FALSE,debug=FALSE,boot.fit=FALSE) {
  minrange <- 5
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
  }
  require("plotrix")
  require("tikzDevice")

  # extract some stuff from cf

  res <- data.frame(name=c(), t1=c(), t2=c(), M=c(), dM=c(),
                Meff=c(), dMeff=c(), P1=c(), dP1=c(),
                P2=c(), dP2=c(), f=c(), df=c(),
                qm1=c(), qm2=c(), kappa=c(),matrixfit.ChiSqr=c(),
                matrixfit.dChiSqr=c(), matrixfit.dof=c(), effectivemass.ChiSqr=c(),
                effectivemass.dChiSqr=c(), effectivemass.dof=c() )

  for( T1 in 3:((cf$Time/2)-1-minrange) ) {
    for( T2 in (T1+minrange):((cf$Time/2)-1) ) {
      cat(T1, "to", T2, "\n")

      cf.matrixfit <- try(matrixfit(cf=cf, t1=T1, t2=T2, parlist=array(c(1,1,1,2,2,2), dim=c(2,3)),
                                matrix.size=3, symmetrise=T, boot.R=cf$boot.R, boot.l=cf$boot.l, useCov=useCov, boot.fit=boot.fit))

      # if the fit fails, we go to the next fit range
      if( any(class(cf.matrixfit) == "try-error" )) next

      cf.matrixfit <- computefps( cf.matrixfit, mu1=q_masses$m1, mu2=q_masses$m2, Kappa=kappa, #normalisation="mpsexplicit",
       disprel=fps.disprel, boot.fit=boot.fit )

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
                matrixfit.dof=cf.matrixfit$dof, matrixfit.Q=cf.matrixfit$Qval,
                effectivemass.ChiSqr=effectivemass.fit$opt.res$value, effectivemass.dChiSqr=effectivemass.dChiSqr,
                effectivemass.dof=effectivemass.fit$dof, effectivemass.Q=effectivemass.fit$Qval ))
    }
  }

  # need to rename object before saving to file so that it can easily be loaded
  # and referred to in an external script
  # we make sure that no "-" characters remain in the name because they would
  # make it very difficult to use the objects in the future
  savename <- sprintf("%s.fitrange",gsub("-","_",name) )
  assign(savename,res)
  save(list=savename,file=sprintf("%s.Rdata",savename))

  # we now remove the outliers using the usual idea of computing quartiles and the interquartile range
  quants <- quantile(res$M)
  tshld.hi <- quants[4] + 1.5*IQR(res$M)
  tshld.lo <- quants[2] - 1.5*IQR(res$M)

  outlier.indices <- which( res$M < tshld.lo | res$M > tshld.hi )

  if( length(res[,1]) == length(outlier.indices) ) {
    warning("For ", name, " no entries remain after removal of outliers! Continuing with full set!")
  } else if( length(outlier.indices) == 0 ) {
    warning("Interquartile range very large, no outliers found!\n")
  } else {
    res <- res[-outlier.indices,]
  }

  # assemble relevant data for convenient plotting
  l.matrix <- list(df=data.frame( val=res$M, t1=res$t1, t2=res$t2,
                               ChiSqr.ov.dof=(res$matrixfit.ChiSqr/res$matrixfit.dof),
                               Q=res$matrixfit.Q),
                               label="$M_{\\mathrm{mtx}}$",
                               name="matrixfit" )

  l.effmass <- list(df=data.frame( val=res$Meff, t1=res$t1, t2=res$t2,
                                ChiSqr.ov.dof=(res$effectivemass.ChiSqr/res$effectivemass.dof),
                                Q=res$effectivemass.Q),
                                label="$M_{\\mathrm{eff}}$",
                                name="effmass" )


  # produce a number of plots relating to the fit range analysis
  temp <- sprintf("%s.fitrange.%s",name,c("tex","pdf","aux","log"))
  tikzfiles <- list(tex=temp[1],pdf=temp[2],aux=temp[3],log=temp[4])
  rm(temp)
  tikz(tikzfiles$tex, standAlone = TRUE, width=5, height=5)

  # colours to add some timeslice information
  colours <- rainbow(n=cf$Time/2)
  # colours to provide information about Q
  Qcolours <- rainbow(n=100)
  
  for( l in list( l.matrix, l.effmass ) ) {
    df <- l$df
    qtyname <- l$name
    label <- l$label
    title <- sprintf("%s %s",qtyname,name)
    # prefix all underscores with a backslash so that latex understands properly
    # note that we need to escape the backslash twice in order to have two backslashes
    # in the resulting string
    title <- gsub("_","\\\\_",title)
    print(title)

    hist(df$val,breaks=40,main=title,xlab=label)

    plot(density(df$val),main=title, xlab=label)
    abline(v=mean(df$val),col='red')
    abline(v=median(df$val),col='blue')

    plot(y=df$ChiSqr.ov.dof,x=df$val,main=title,
         col=colours[df$t2],ylab="$\\chi^2 / \\mathrm{d.o.f}$",xlab=label
        )

    plot(y=df$Q,x=df$val,main=title,col=colours[df$t2],ylab='$Q$',xlab=label)

    weighted.hist(x=df$val,w=df$Q,breaks=40,main=paste("weighted",title),xlab=label)

    # save some more lines by doing two sets of plots in one go
    for( dat in list( list(qty=df$Q,lab='Q'), list(qty=df$val,lab=label) ) ) {
      plot(x=df$t2-df$t1,y=dat$qty,main=title,col=colours[df$t2],xlab="$t_f - t_i$",ylab=dat$lab)
      plot(x=df$t1+df$t2,y=df$Q,main=title,col=colours[df$t2],xlab="$t_i + t_f$",ylab=dat$lab)
      plot(x=df$t1,y=dat$qty,main=title,col=colours[df$t2],xlab="$t_i$",ylab=dat$lab)
      plot(x=df$t2,y=dat$qty,main=title,col=colours[df$t1],xlab="$t_f$",ylab=dat$lab)
    }

    boxplot(df$val,main=title,ylab=label)
  }

  # weighted mean and variance for the matrixfit masses
  matrixfit.mustar <- weighted.mean(l.matrix$df$val,l.matrix$df$Q)
  matrixfit.varstar <- weighted.variance(l.matrix$df$val,l.matrix$df$Q)

  # weighted mean and variance for the effective masses
  effective.mustar <- weighted.mean(l.effmass$df$val,l.effmass$df$Q)
  effective.varstar <- weighted.variance(l.effmass$df$val,l.effmass$df$Q)

  cat("Effective mass:", effective.mustar, sqrt(effective.varstar), "\n")
  cat("Matrixfit:", matrixfit.mustar, sqrt(matrixfit.varstar), "\n")

  effseq <- seq(from=min(l.effmass$df$val),to=max(l.effmass$df$val),length.out=500)
  effnorm <- dnorm(x=effseq,mean=effective.mustar,sd=sqrt(effective.varstar))
  matrixseq <- seq(from=min(l.matrix$df$val),to=max(l.matrix$df$val),length.out=500)
  matrixnorm <- dnorm(x=matrixseq,mean=matrixfit.mustar,sd=sqrt(matrixfit.varstar))

  ylims <- c(0,max( c(effnorm,matrixnorm) ) )

  plot(y=effnorm, x=effseq,lwd=3,type='l',col='blue',
      xlim=effective.mustar+4*c(-1,1)*sqrt(effective.varstar),ylim=ylims,xlab="mass",ylab="",
      main=sprintf("weighted gaussians %s", gsub("_","\\\\_",name)))

  lines(y=matrixnorm, x=matrixseq,lwd=3,col='red',type='l')

  plot(y=l.effmass$df$val,x=l.matrix$df$val, main="effective vs. matrixfit",
       ylab="effective mass", xlab="matrixfit mass",col=colours[l.effmass$df$t2])

  # plot chosen fit ranges, indicating the Q value with a
  plot(x=l.matrix$df$t1,y=l.matrix$df$t2,main=sprintf("chosen fitranges %s", gsub("_","\\\\_",name)),
       xlab=expression(t[i]),ylab=expression(t[f]),col=Qcolours[as.integer(100*l.matrix$df$Q)])

  dev.off()
  tools::texi2dvi(tikzfiles$tex,pdf=T)                                                                                                                                                                         
  # use pdfcrop tool to remove plot borders
  command <- sprintf("pdfcrop %s %s",tikzfiles$pdf,tikzfiles$pdf)
  system(command)
  # remove temporary files 
  command <- sprintf("rm %s %s %s", tikzfiles$tex, tikzfiles$log, tikzfiles$aux)
  system(command)

  return(invisible(res))
}
