# function for the automated analysis of connected meson two-point functions
# with the aim of extracting the (ground-state) masses and decay constants

# to adjust it for a given purpose the 
#   *_masses
#   mass_comb
#   dirs
#   analyses
# objects must be edited

analysis_conn_meson_2pt <- function(analyses_to_be_done_input,kappa,boot.R=400,boot.l=10,debug=F,pause=F,skip=0) {
  # masses to be used in this analysis
  light_masses <- c(0.0009)
  strange_masses <- c(0.0238,0.0245,0.0252,0.0259)
  charm_masses <- c(0.2822,0.294,0.3058,0.3176)
  
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
                ll_n=sprintf("lln_u_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ls_c=sprintf("ls_u_%g-sp_%g",mass_comb$ls$m1,mass_comb$ls$m2),
                lc_c=sprintf("lc_u_%g-cp_%g",mass_comb$lc$m1,mass_comb$lc$m2),
                sc_c=sprintf("sc_sp_%g-cp_%g",mass_comb$sc$m1,mass_comb$sc$m2) )

  analyses <- NULL
  # the "analyses" list drives the automation process, the first four elements are self-explanatory
  #   t1,t2 are the fit range, t1_plot,t2_plot the x plotting limits
  #   basename is the correlator filename "prefix"
  #   observable is a numerical vector identifying which "gamma combinations" will be fitted (in the CMI format)
  #   sign is a numerical vector indicating whether the correlator is of "cosh" (+1) or "sinh" (-1) form
  analyses[[1]] <- list( dirs=dirs$ll_c, name="pion_charged", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=12, t2=47, t1_plot=8, t2_plot=48, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[2]] <- list( dirs=dirs$ls_c, name="kaon_charged", mass_diagonal=F, q_masses=mass_comb$ls,
                         t1=12, t2=47, t1_plot=8, t2_plot=48, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[3]] <- list( dirs=dirs$lc_c, name="D_charged", mass_diagonal=F, q_masses=mass_comb$lc,
                         t1=12, t2=47, t1_plot=8, t2_plot=48, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[4]] <- list( dirs=dirs$sc_c, name="Ds_charged", mass_diagonal=F, q_masses=mass_comb$sc,
                         t1=15, t2=47, t1_plot=8, t2_plot=48, basename="outprcv.", observables=c(1), sign=c(1) )

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
                                  boot.R=boot.R, boot.l=boot.l)
      analysis_results <- rbind(analysis_results, result)
    }
  }

  write.csv(analysis_results, file="meson_2pt.csv")
}

do_meson_analysis <- function(directory,name,t1,t2,t1_plot,t2_plot,
                              debug=F,pause=F,basename="outprcv.",observable=c(1),sign=c(+1),skip=0,boot.R=400,boot.l=10,
                              kappa, q_masses) {

  if(debug) {
    cat(sprintf("Doing %s meson analysis in %s, observable %d\n",name,directory,observable))
  }

  files <- getorderedfilelist(basename=basename,path=directory)
# prepend directory name
  if(debug) {
    cat("Reading", basename, "correlators from",directory,"\n")
  }
  cmicor <- readcmidatafiles(files[skip+1:length(files)])
  
  if(debug) {
    cat("Performing matrixfit\n")
  }
  meson.cor <- extract.obs(cmicor, vec.obs=observable,sign.vec=sign)
  meson.cor <- bootstrap.cf(meson.cor,boot.R=boot.R,boot.l=boot.l,seed=232415)
  meson.cor.matrixfit <- matrixfit(meson.cor, t1=t1, t2=t2, symmetrise=T, boot.R=boot.R, boot.l=boot.l, useCov=F)
  meson.cor.matrixfit <- computefps( meson.cor.matrixfit, mu1=q_masses$m1, mu2=q_masses$m2, Kappa=kappa ) 
  summary(meson.cor.matrixfit)

  if(debug) {
    cat("Performing bootstrapped effective mass fit\n")
  }
  meson.cor.effectivemass <- bootstrap.effectivemass(meson.cor, type="acosh", boot.R=boot.R, boot.l=boot.l)
  meson.cor.effectivemass <- fit.effectivemass(meson.cor.effectivemass, t1=t1, t2=t2, useCov=F)
  summary(meson.cor.effectivemass)
 
  rval <- data.frame(name=name, M=meson.cor.matrixfit$opt.res$par[1], dM=sd(meson.cor.matrixfit$opt.tsboot[1,]), 
                Meff=meson.cor.effectivemass$opt.res$par[1], dMeff=sd(meson.cor.effectivemass$massfit.tsboot[,1]), 
                P1=meson.cor.matrixfit$opt.res$par[2], dP1=sd(meson.cor.matrixfit$opt.tsboot[2,]), 
                P2=meson.cor.matrixfit$opt.res$par[3], dP2=sd(meson.cor.matrixfit$opt.tsboot[3,]), 
                f=meson.cor.matrixfit$fps, df=sd(meson.cor.matrixfit$fps.tsboot), 
                qm1=q_masses$m1, qm2=q_masses$m2, kappa=kappa,
                matrixfit.ChiSqr=meson.cor.matrixfit$opt.res$value, matrixfit.dChiSqr=sd(meson.cor.matrixfit$opt.tsboot[,4]),
                matrixfit.dof=meson.cor.matrixfit$dof,
                effectivemass.ChiSqr=meson.cor.effectivemass$opt.res$value, effectivemass.dChiSqr=sd(meson.cor.effectivemass$massfit.tsboot[,2]),
                effectivemass.dof=meson.cor.effectivemass$dof )
 
  if(pause) {
    readline(prompt="Press any key to continue.")
  }

  # estimate some boundaries for the effective mass plot
  ywidth <- 5*max( meson.cor.effectivemass$deffMass[t1:t2], na.rm=T )
  ymid <- mean( meson.cor.effectivemass$effMass[t1:t2], na.rm=T )

  plotfilename <- sprintf("plots/%s_%s.pdf",name,directory)
  pdf(plotfilename,onefile=T,title=paste(name, directory, sep=" "))
  plot(meson.cor.effectivemass, xlab=c("t/a"), ylab="aM", main=paste(name,directory,sep=" "), 
       xlim=c(t1_plot,t2_plot), ylim=c(ymid-ywidth,ymid+ywidth) )
  plot(meson.cor.matrixfit, xlab="t/a", ylab="C(t)",main=paste(name, directory,sep=" "))     
  dev.off();
  
  cat("\n")

  return(invisible( rval ))
}

