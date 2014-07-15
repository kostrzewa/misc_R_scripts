# function for the automated analysis of connected meson two-point functions
# with the aim of extracting the (ground-state) masses and decay constants

# to adjust it for a given purpose the 
#   *_masses
#   mass_comb
#   dirs
#   analyses
# objects must be edited

analysis_conn_meson_2pt <- function(analyses_to_be_done_input,kappa,boot.R=400,boot.l=20,debug=F,pause=F,skip=0,seed=12345,useCov=F,read.cor=T,study.fitrange=F) {
  # masses to be used in this analysis
#  light_masses <- c(0.0009)
  strange_masses <- c(0.0238,0.0245,0.0252,0.0259)
  charm_masses <- c(0.2822,0.294,0.3058,0.3176)
  light_masses <- c(0.0009)
  
  # for testing purposes
  #strange_masses <- c(0.021,0.023,0.025,0.027,0.031)
  #charm_masses <- c(0.18,0.25,0.27,0.29,0.31,0.36,0.37,0.38)

  # combinations of these masses
  mass_comb <- list( ll=data.frame( m1=light_masses, m2=light_masses ),
                     ls=expand.grid( m1=light_masses, m2=strange_masses),
                     lc=expand.grid( m1=light_masses, m2=charm_masses),
                     sc=expand.grid( m1=strange_masses, m2=charm_masses) )
  
  # (sub-)directories where the correlators will be loaded from
  dirs <- list( ll_c=sprintf("llc_u_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ll_n_ud=sprintf("lln_u_%g-d_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ll_n_du=sprintf("lln_d_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ls_c=sprintf("ls_u_%g-sp_%g",mass_comb$ls$m1,mass_comb$ls$m2),
                lc_c=sprintf("lc_u_%g-cp_%g",mass_comb$lc$m1,mass_comb$lc$m2),
                sc_c=sprintf("sc_sp_%g-cp_%g",mass_comb$sc$m1,mass_comb$sc$m2) )

  analyses <- NULL
  # the "analyses" list drives the automation process, the first four elements are self-explanatory
  #   t1,t2 are the fit range, t1_plot,t2_plot the x plotting limits
  #   basename is the correlator filename "prefix"
  #   observable is a numerical vector identifying which "gamma combinations" will be fitted (in the CMI format)
  #   sign is a numerical vector indicating whether the correlator is of "cosh" (+1) or "sinh" (-1) form
  analyses[[1]] <- list( dirs=dirs$ll_c, name="ll_c", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=10, t2=42, t1_plot=5, t2_plot=48, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[2]] <- list( dirs=dirs$ls_c, name="ls_c", mass_diagonal=F, q_masses=mass_comb$ls,
                         t1=12, t2=40, t1_plot=8, t2_plot=48, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[3]] <- list( dirs=dirs$lc_c, name="lc_c", mass_diagonal=F, q_masses=mass_comb$lc,
                         t1=15, t2=30, t1_plot=8, t2_plot=48, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[4]] <- list( dirs=dirs$sc_c, name="sc_c", mass_diagonal=F, q_masses=mass_comb$sc,
                         t1=16, t2=35, t1_plot=8, t2_plot=48, basename="outprcv.", observables=c(1), sign=c(1) )
  
  analyses[[5]] <- list( dirs=dirs$ll_n_ud, name="ll_n_ud", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=11, t2=47, t1_plot=5, t2_plot=48, basename="outprcvn.", observable=c(5), sign=c(-1) )
  
  analyses[[6]] <- list( dirs=dirs$ll_n_du, name="ll_n_du", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=11, t2=47, t1_plot=5, t2_plot=48, basename="outprcvn.", observable=c(5), sign=c(-1) )

  analyses_to_be_done <- NULL
  if( missing(analyses_to_be_done_input) ) {
    cat("Analyses to be done missing, doing all!\n")
    analyses_to_be_done <- seq(1,length(analyses))
  } else {
    analyses_to_be_done <- analyses_to_be_done_input
  }

  analysis_results <- NULL
  for( ctr_analyses in analyses_to_be_done ) {
    for( ctr_dirs in seq(1,length( analyses[[ctr_analyses]]$dirs ) ) ) {
      if(!file.exists( analyses[[ctr_analyses]]$dirs[ctr_dirs] )) {
        cat("## Skipping", analyses[[ctr_analyses]]$dirs[ctr_dirs], "because it doesn't exist!\n")
        next
      }

      result <- do_meson_analysis(directory=analyses[[ctr_analyses]]$dirs[ctr_dirs], name=analyses[[ctr_analyses]]$name,debug=debug, pause=pause,
                                  basename=analyses[[ctr_analyses]]$basename, t1=analyses[[ctr_analyses]]$t1, t2=analyses[[ctr_analyses]]$t2,
                                  t1_plot=analyses[[ctr_analyses]]$t1_plot,t2_plot=analyses[[ctr_analyses]]$t2_plot,
                                  observable=analyses[[ctr_analyses]]$observable , sign=analyses[[ctr_analyses]]$sign,
                                  skip=skip, kappa=kappa, q_masses=analyses[[ctr_analyses]]$q_masses[ctr_dirs,],
                                  boot.R=boot.R, boot.l=boot.l, seed=seed, useCov=useCov,read.cor=read.cor)
          
      analysis_results <- rbind(analysis_results, result)
    }
  }

  write.csv(analysis_results, file="meson_2pt.csv")
}

do_meson_analysis <- function(directory,name,t1,t2,t1_plot,t2_plot,kappa,q_masses,
                              debug=F,pause=F,basename="outprcv.",observable=c(1),sign=c(+1),skip=0,boot.R=400,boot.l=10,
                              seed=12345,useCov=F,read.cor=T,study.fitrange=F) {

  if(debug) {
    cat(sprintf("Doing %s meson analysis in %s, observable %d\n",name,directory,observable))
  }

  files <- getorderedfilelist(basename=basename,path=directory)
  
  cmicor <- NULL
  if(read.cor) {
    # prepend directory name
    if(debug) {
      cat("Reading", basename, "correlators from",directory,"\n")
    }
    cmicor <- readcmidatafiles(files[skip+1:length(files)])
    save(cmicor,file=sprintf("%s.cor.Rdata",directory))
  } else {
    load(sprintf("%s.cor.Rdata",directory))
  }
  
  if(debug) {
    cat("Performing matrixfit\n")
  }

  meson.cor <- NULL
  meson.cor.effectivemass <- NULL

  if(study.fitrange) {
    require("plotrix")
    # we start off with the minimum possible number of bootstrap samples because in the first run through the data
    # we are only interested in the effect of the choice of fit range
    #
    meson.cor <- extract.obs(cmicor, vec.obs=observable,sign.vec=sign)
    meson.cor <- bootstrap.cf(meson.cor,boot.R=5,boot.l=1,seed=12345)
    meson.cor.effectivemass <- bootstrap.effectivemass(meson.cor, type="solve", boot.R=5, boot.l=1, seed=12345)
    
    # we will scan the possible fit ranges to estimate the effect of the choice of fit range
    # on the value of the mass and decay constant
    # However, we will return the result of the fit for the specified t1 and t2

    results <- list()
    i <- 1
    count <- 1

    for( T1 in 10:((meson.cor$Time/2)-1-5) ) {
      for( T2 in (T1+5):((meson.cor$Time/2)-1) ) {
        cat(T1, "to", T2, "\n")
        count <- count + 1
        
        meson.cor.matrixfit <- try(matrixfit(meson.cor, t1=T1, t2=T2, symmetrise=T, boot.R=boot.R, boot.l=boot.l, useCov=useCov))
        # if the fit fails, we go to the next fit range
        if( any(class(meson.cor.matrixfit) == "try-error" )) next
        meson.cor.matrixfit <- computefps( meson.cor.matrixfit, mu1=q_masses$m1, mu2=q_masses$m2, Kappa=kappa, reduce.latarts=F )
    
        meson.cor.effectivemass <- try(fit.effectivemass(meson.cor.effectivemass, t1=T1, t2=T2, useCov=useCov, replace.na=TRUE))
        if( any(class(meson.cor.effectivemass) == "try-error" )) next
        
        results[[i]] <- data.frame(name=name, t1=T1, t2=T2, M=meson.cor.matrixfit$opt.res$par[1], dM=sd(meson.cor.matrixfit$opt.tsboot[1,]), 
                  Meff=meson.cor.effectivemass$opt.res$par[1], dMeff=sd(meson.cor.effectivemass$massfit.tsboot[,1]), 
                  P1=meson.cor.matrixfit$opt.res$par[2], dP1=sd(meson.cor.matrixfit$opt.tsboot[2,]), 
                  P2=meson.cor.matrixfit$opt.res$par[3], dP2=sd(meson.cor.matrixfit$opt.tsboot[3,]), 
                  f=meson.cor.matrixfit$fps, df=sd(meson.cor.matrixfit$fps.tsboot), 
                  qm1=q_masses$m1, qm2=q_masses$m2, kappa=kappa,
                  matrixfit.ChiSqr=meson.cor.matrixfit$opt.res$value, matrixfit.dChiSqr=sd(meson.cor.matrixfit$opt.tsboot[,4]),
                  matrixfit.dof=meson.cor.matrixfit$dof,
                  effectivemass.ChiSqr=meson.cor.effectivemass$opt.res$value, effectivemass.dChiSqr=sd(meson.cor.effectivemass$massfit.tsboot[,2]),
                  effectivemass.dof=meson.cor.effectivemass$dof )
        
        i <- i+1
      }
    }
    
    rm(i)

    # need to rename object before saving to file
    savename <- sprintf("%s.fitrange",directory)
    assign(savename,results)
    save(list=savename,file=sprintf("$s.Rdata",savename)) 
    
    matrixfit.masses <- NULL
    effective.masses <- NULL
    ChiSqr.ov.dof <- NULL
    effective.ChiSqr.ov.dof <- NULL
    matrixfit.Q <- NULL
    effective.Q <- NULL
    t1s <- NULL
    t2s <- NULL
    deltat <- NULL
    for( i in 1:length(results) ) {
      criterion <- results[[i]]$matrixfit.ChiSqr/results[[i]]$matrixfit.dof
     # if( criterion <= 3 ) {
        t1s <- c(t1s,results[[i]]$t1)
        t2s <- c(t2s,results[[i]]$t2)
        deltat <- c(deltat,results[[i]]$t2-results[[i]]$t1)
        matrixfit.masses <- c(matrixfit.masses,results[[i]]$M)
        effective.masses <- c(effective.masses,results[[i]]$Meff)
        ChiSqr.ov.dof <- c(ChiSqr.ov.dof,criterion)
        #matrixfit.Q <- c(matrixfit.Q,(1-pgamma(shape=(results[[i]]$matrixfit.dof-2)/2,q=results[[i]]$matrixfit.ChiSqr/2)))
        matrixfit.Q <- c(matrixfit.Q,(1-pgamma(shape=(results[[i]]$matrixfit.dof)/2,q=results[[i]]$matrixfit.ChiSqr/2)))
        effective.ChiSqr.ov.dof <- c(effective.ChiSqr.ov.dof,results[[i]]$effectivemass.ChiSqr/results[[i]]$effectivemass.dof)
        #effective.Q <- c(effective.Q,(1-pgamma(shape=(results[[i]]$effectivemass.dof-2)/2,q=results[[i]]$effectivemass.ChiSqr/2)))
        effective.Q <- c(effective.Q,(1-pgamma(shape=(results[[i]]$effectivemass.dof)/2,q=results[[i]]$effectivemass.ChiSqr/2)))
      #}
    }

    # weighted mean and variance for the matrixfit masses
    matrixfit.mustar <- weighted.mean(matrixfit.masses,matrixfit.Q)
    matrixfit.varstar <- weighted.variance(matrixfit.masses,matrixfit.Q)
    
    # weighted mean and variance for the effective masses
    effective.mustar <- weighted.mean(effective.masses,effective.Q)
    effective.varstar <- weighted.variance(effective.masses,effective.Q)
    
    cat("Effective mass:", effective.mustar, sqrt(effective.varstar), "\n")
    cat("Matrixfit:", matrixfit.mustar, sqrt(matrixfit.varstar), "\n")

    # produce a number of plots relating to the fit range analysis 
    pdf(sprintf("%s.fitrange.pdf",directory),title=directory) 

    hist(matrixfit.masses,breaks=40,main=paste("histogram matrixfit",directory))
    plot(density(matrixfit.masses),main=paste("matrixfit",directory))
    abline(v=mean(matrixfit.masses),col='red')
    abline(v=median(matrixfit.masses),col='blue')
    plot(y=ChiSqr.ov.dof,x=matrixfit.masses,main=paste("matrixfit",directory))
    plot(y=matrixfit.Q,x=matrixfit.masses,main=paste("matrixfit",directory))
    weighted.hist(x=matrixfit.masses,w=matrixfit.Q,breaks=40,main=paste("weighted hist. matrixfit",directory))
    plot(x=deltat,y=matrixfit.Q,main=paste("matrixfit",directory))
    plot(x=t1s+t2s,y=matrixfit.Q,main=paste("matrixfit",directory))
    plot(x=t1s,y=matrixfit.Q,main=paste("matrixfit",directory))
    plot(x=t2s,y=matrixfit.Q,main=paste("matrixfit",directory))
    plot(x=deltat,y=matrixfit.masses,main=paste("matrixfit",directory))
    plot(x=t1s+t2s,y=matrixfit.masses,main=paste("matrixfit",directory))
    plot(x=t1s,y=matrixfit.masses,main=paste("matrixfit",directory))
    plot(x=t2s,y=matrixfit.masses,main=paste("matrixfit",directory))
    boxplot(matrixfit.masses,main=paste("matrixfit",directory))
    
    hist(effective.masses,breaks=40,main=paste("histogram effective",directory))
    plot(density(effective.masses),main=paste("effective",directory))
    abline(v=mean(effective.masses),col='red')
    abline(v=median(effective.masses),col='blue')
    plot(y=effective.ChiSqr.ov.dof,x=effective.masses,main=paste("effective",directory))
    plot(y=effective.Q,x=effective.masses,main=paste("effective",directory))
    weighted.hist(x=effective.masses,w=effective.Q,breaks=40,main=paste("weighted hist. effective",directory))
    plot(x=deltat,y=effective.Q,main=paste("effective",directory))
    plot(x=t1s+t2s,y=effective.Q,main=paste("effective",directory))
    plot(x=t1s,y=effective.Q,main=paste("effective",directory))
    plot(x=t2s,y=effective.Q,main=paste("effective",directory))
    plot(x=deltat,y=effective.masses,main=paste("effective",directory))
    plot(x=t1s+t2s,y=effective.masses,main=paste("effective",directory))
    plot(x=t1s,y=effective.masses,main=paste("effective",directory))
    plot(x=t2s,y=effective.masses,main=paste("effective",directory))
    boxplot(effective.masses,main=paste("effective",directory))

    plot(y=dnorm(x=seq(from=min(matrixfit.masses),to=max(matrixfit.masses),length.out=500),mean=effective.mustar,sd=sqrt(effective.varstar)),
         x=seq(from=min(matrixfit.masses),to=max(matrixfit.masses),length.out=500),lwd=3,type='l',col='blue',xlim=c(0.059,0.064))
    lines(y=dnorm(x=seq(from=min(matrixfit.masses),to=max(matrixfit.masses),length.out=500),mean=matrixfit.mustar,sd=sqrt(matrixfit.varstar)),
         x=seq(from=min(matrixfit.masses),to=max(matrixfit.masses),length.out=500),lwd=3,col='red',type='l')
    plot(y=effective.masses,x=matrixfit.masses, main="effective vs. matrixfit", ylab="effective mass", xlab="matrixfit mass") 
    plot(x=t1s,y=t2s,main=directory)

    dev.off()
  }


  # now we redo the bootstrapping with the correct values and do the fit with our best estimate of the optimal fit range
  # which has been given as the parameter
  if(!study.fitrange) {
    meson.cor <- extract.obs(cmicor, vec.obs=observable,sign.vec=sign)
  }
  meson.cor <- bootstrap.cf(meson.cor,boot.R=boot.R,boot.l=boot.l,seed=seed)
  save.matrixfit <- matrixfit(meson.cor, t1=t1, t2=t2, symmetrise=T, boot.R=boot.R, boot.l=boot.l, useCov=useCov)
  save.matrixfit <- computefps( save.matrixfit, mu1=q_masses$m1, mu2=q_masses$m2, Kappa=kappa, reduce.latarts=F )
  save.effectivemass <- bootstrap.effectivemass(meson.cor, type="solve", boot.R=boot.R, boot.l=boot.l, seed=seed)
  if(debug) {
    cat("Performing bootstrapped effective mass fit\n")
  }
  save.effectivemass <- fit.effectivemass(save.effectivemass, t1=t1, t2=t2, useCov=useCov, replace.na=TRUE)
  
  # we need to rename the object before saving it to file so we can 
  # "load" matrixfit objects for multiple mesons at the same time
  savename <- sprintf("%s.m%g.m%g.matrixfit",name,q_masses$m1,q_masses$m2)
  assign(savename,save.matrixfit)
  save(list=savename,file=sprintf("%s.Rdata",savename))
  summary(save.matrixfit)

  # also rename the effectivemass object 
  savename <- sprintf("%s.m%g.m%g.effectivemass",name,q_masses$m1,q_masses$m2)
  assign(savename,save.effectivemass)
  save(list=savename,file=sprintf("%s.Rdata",savename))
  summary(save.effectivemass)

  rval <- data.frame(name=name, t1=t1, t2=t2, M=save.matrixfit$opt.res$par[1], dM=sd(save.matrixfit$opt.tsboot[1,]), 
                Meff=save.effectivemass$opt.res$par[1], dMeff=sd(save.effectivemass$massfit.tsboot[,1]), 
                P1=save.matrixfit$opt.res$par[2], dP1=sd(save.matrixfit$opt.tsboot[2,]), 
                P2=save.matrixfit$opt.res$par[3], dP2=sd(save.matrixfit$opt.tsboot[3,]), 
                f=save.matrixfit$fps, df=sd(save.matrixfit$fps.tsboot), 
                qm1=q_masses$m1, qm2=q_masses$m2, kappa=kappa,
                matrixfit.ChiSqr=save.matrixfit$opt.res$value, matrixfit.dChiSqr=sd(save.matrixfit$opt.tsboot[,4]),
                matrixfit.dof=save.matrixfit$dof,
                effectivemass.ChiSqr=save.effectivemass$opt.res$value, effectivemass.dChiSqr=sd(save.effectivemass$massfit.tsboot[,2]),
                effectivemass.dof=save.effectivemass$dof )

  if(pause) {
    readline(prompt="Press any key to continue.")
  }

  # estimate some boundaries for the effective mass plot
  ymin <- min( save.effectivemass$effMass[t1_plot:t2_plot], na.rm=T )
  ymax <- max( save.effectivemass$effMass[t1_plot:t2_plot], na.rm=T )

  if( ymin < 0 || is.na(ymin) ) {
    ymin <- 0
  }

  if( ymax > 1.5 || is.na(ymax) ) {
    ymax <- 1.5
  }

  plotfilename <- sprintf("plots/%s.pdf",directory)
  pdf(plotfilename,onefile=T,title=paste(name, directory, sep=" "))
  plot(save.effectivemass, xlab=c("t/a"), ylab="aM", main=paste(name,directory,sep=" "), 
       xlim=c(t1_plot,t2_plot), ylim=c(ymin,ymax) )
  plot(save.matrixfit, xlab="t/a", ylab="C(t)",main=paste(name, directory,sep=" "))     

  # add another plot which shows the difference between the fitted correlation function
  # and the original data
  Cor <- function(t) {
    m <- save.matrixfit$opt.res$par[1] 
    a <- 0.5*save.matrixfit$opt.res$par[2]^2
    Time <- meson.cor$Time 
    return( a * ( exp ( -(m*t) ) + exp( -m*(Time-t) ) ) )
  }
 
  trange <- (t1_plot:t2_plot)
  cval <- Cor(0:(meson.cor$Time/2))
  cvalarray <- t(array( cval,dim=c(length(cval),length(meson.cor$cf[,1]))))
  
  plotwitherror(x=trange,y=meson.cor$cf0[ trange+1 ]/Cor(trange),
                dy=apply(X=(meson.cor$cf[,trange+1]/cvalarray[,trange+1]),MARGIN=2,FUN=sd)/sqrt(length(meson.cor$cf[,1])),
                main=expression(C[i]/C(t)),ylab=NA,xlab="t/a")
  abline(h=1,col="red")
  abline(v=c(t1,t2),col="blue")
  
  dev.off();
  
  cat("\n")

  return(invisible( rval ))
}
