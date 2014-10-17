# function for the automated analysis of connected meson two-point functions
# with the aim of extracting the (ground-state) masses and decay constants

# to adjust it for a given purpose the 
#   *_masses
#   mass_comb
#   dirs
#   analyses
# objects must be edited

source("~/code/R/misc_R_scripts/analysis_conn_meson_2pt/meson_2pt_study_fitrange.R")

analysis_conn_meson_2pt <- function(analyses_to_be_done_input,kappa,boot.R=400,boot.l=20,debug=F,pause=F,skip=0,seed=12345,useCov=F,read.cor=T,study.fitrange=F) {
  # masses to be used in this analysis
  strange_masses <- c(0.0224,0.0231,0.0238,0.0245,0.0252,0.0259)
  charm_masses <- c(0.2586,0.2704,0.2822,0.294,0.3058,0.3176)
  light_masses <- c(0.003)
  
  # combinations of these masses
  mass_comb <- list( ll=data.frame( m1=light_masses, m2=light_masses ),
                     ls=expand.grid( m1=light_masses, m2=strange_masses),
                     lc=expand.grid( m1=light_masses, m2=charm_masses),
                     sc=expand.grid( m1=strange_masses, m2=charm_masses) )
  
  correlators_dir <- "correlators/"

  # (sub-)directories where the correlators will be loaded from
  dirs <- list( ll_c=sprintf("conn_llc_u_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ll_n_ud=sprintf("conn_lln_u_%g-d_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ll_n_du=sprintf("conn_lln_d_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ls_c=sprintf("conn_ls_u_%g-sp_%g",mass_comb$ls$m1,mass_comb$ls$m2),
                lc_c=sprintf("conn_lc_u_%g-cp_%g",mass_comb$lc$m1,mass_comb$lc$m2),
                sc_c=sprintf("conn_sc_sp_%g-cp_%g",mass_comb$sc$m1,mass_comb$sc$m2) )

  analyses <- NULL
  # the "analyses" list drives the automation process, the first four elements are self-explanatory
  #   t1,t2 are the fit range, t1_plot,t2_plot the x plotting limits
  #   basename is the correlator filename "prefix"
  #   observable is a numerical vector identifying which "gamma combinations" will be fitted (in the CMI format)
  #   sign is a numerical vector indicating whether the correlator is of "cosh" (+1) or "sinh" (-1) form
  analyses[[1]] <- list( dirs=dirs$ll_c, name="ll_c", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=10, t2=23, t1_plot=5, t2_plot=24, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[2]] <- list( dirs=dirs$ls_c, name="ls_c", mass_diagonal=F, q_masses=mass_comb$ls,
                         t1=12, t2=20, t1_plot=8, t2_plot=24, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[3]] <- list( dirs=dirs$lc_c, name="lc_c", mass_diagonal=F, q_masses=mass_comb$lc,
                         t1=15, t2=20, t1_plot=8, t2_plot=24, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[4]] <- list( dirs=dirs$sc_c, name="sc_c", mass_diagonal=F, q_masses=mass_comb$sc,
                         t1=16, t2=20, t1_plot=8, t2_plot=24, basename="outprcv.", observables=c(1), sign=c(1) )
  
  analyses[[5]] <- list( dirs=dirs$ll_n_ud, name="ll_n_ud", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=11, t2=23, t1_plot=5, t2_plot=24, basename="outprcvn.", observable=c(5), sign=c(-1) )
  
  analyses[[6]] <- list( dirs=dirs$ll_n_du, name="ll_n_du", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=11, t2=23, t1_plot=5, t2_plot=24, basename="outprcvn.", observable=c(5), sign=c(-1) )

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
                                  boot.R=boot.R, boot.l=boot.l, seed=seed, useCov=useCov, read.cor=read.cor, study.fitrange=study.fitrange
                                 )
          
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
    cat("Extracting observable\n")
  }
  meson.cor <- extract.obs(cmicor, vec.obs=observable,sign.vec=sign)

  # average LF and FL
  meson.cor <- avg.ls.cf(meson.cor)
  
  if(debug) {
    cat("Bootsrapping cf and effectivemass\n")
  }
  meson.cor <- bootstrap.cf(meson.cor,boot.R=boot.R,boot.l=boot.l,seed=seed)
  meson.cor.effectivemass <- bootstrap.effectivemass(meson.cor, type="solve", boot.R=boot.R, boot.l=boot.l, seed=seed)
    
  if(study.fitrange) {
    meson_2pt_study_fitrange(cf=meson.cor,effmass=meson.cor.effectivemass,name=directory,
                                 kappa=kappa,useCov=useCov,q_masses=q_masses,boot.fit=F,
                                 debug=debug)
  }

  if(debug) {
    cat("Performing matrixfit\n")
  }
  # now we fully bootstrap the fits with the fit range that we deem optimal (from the inputs t1, t2)
  save.matrixfit <- matrixfit(meson.cor, t1=t1, t2=t2, symmetrise=T, parlist=array(c(1,1,1,2,2,2), dim=c(2,3)),
                              matrix.size=3, boot.R=boot.R, boot.l=boot.l, useCov=useCov, boot.fit=T )

  if(debug) {
    cat("Extracting decay constant\n")
  }
  save.matrixfit <- computefps( save.matrixfit, mu1=q_masses$m1, mu2=q_masses$m2, Kappa=kappa, disprel='continuum', boot.fit=T )
  
  if(debug) {
    cat("Performing bootstrapped effective mass fit\n")
  }
  save.effectivemass <- bootstrap.effectivemass(meson.cor, type="solve", boot.R=boot.R, boot.l=boot.l, seed=seed )
  save.effectivemass <- fit.effectivemass(save.effectivemass, t1=t1, t2=t2, useCov=useCov, replace.na=TRUE, boot.fit=T)
  
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
                matrixfit.ChiSqr=save.matrixfit$opt.res$value, 
                matrixfit.dChiSqr=sd(save.matrixfit$opt.tsboot[length(save.matrixfit$opt.tsboot[,1]),]),
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
       xlim=c(t1_plot,t2_plot), ylim=c(ymin,ymax), pch=4)
  plot(save.matrixfit,plot.errorband=TRUE, xlab="t/a", ylab="C(t)",main=paste(name, directory,sep=" "), pch=4)     

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
  
  ## add another two plots which shows the distribution of ChiSqr over the bootstrap samples
  matrixfit.ChiSqr.tsboot <- save.matrixfit$opt.tsboot[length(save.matrixfit$opt.tsboot[,1]),]
  matrixfit.ChiSqr.x <- seq(min(matrixfit.ChiSqr.tsboot),max(matrixfit.ChiSqr.tsboot),1)
  hist(matrixfit.ChiSqr.tsboot,breaks=40,main=paste("matrixfit.ChiSqr",directory,sep=" "),freq=FALSE)
  lines(x=matrixfit.ChiSqr.x,y=dchisq(x=matrixfit.ChiSqr.x,df=save.matrixfit$dof),lwd=3,col='red')

  effectivemass.ChiSqr.tsboot <- save.effectivemass$massfit.tsboot[,2]
  effectivemass.ChiSqr.x <- seq(min(effectivemass.ChiSqr.tsboot),max(effectivemass.ChiSqr.tsboot),1)
  hist(effectivemass.ChiSqr.tsboot,breaks=40,main=paste("effectivemass.ChiSqr",directory,sep=" "),freq=FALSE)
  lines(x=effectivemass.ChiSqr.x,y=dchisq(x=effectivemass.ChiSqr.x,df=save.effectivemass$dof),lwd=3,col='blue')

  dev.off();
  
  cat("\n")

  return(invisible( rval ))
}
