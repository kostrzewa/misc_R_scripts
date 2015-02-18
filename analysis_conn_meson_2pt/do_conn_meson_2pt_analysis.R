source("~/code/R/misc_R_scripts/analysis_conn_meson_2pt/meson_2pt_study_fitrange.R")

do_conn_meson_2pt_analysis <- function(directory,name,t1,t2,t1_plot,t2_plot,kappa,q_masses,end=-1,
                              debug=F,pause=F,basename="outprcv.",observable=c(1),sign=c(+1),skip=0,boot.R=400,boot.l=10,
                              seed=12345,useCov=F,read.cor=T,study.fitrange=F,fps.disprel='continuum') {

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
    if(end<0  || end > length(files) ){
      end <- length(files)
    }

    cmicor <- readcmidatafiles(files[(skip+1):end],verbose=debug)
    save(cmicor,file=sprintf("%s.cor.Rdata",directory))
  } else {
    archivename <- sprintf("%s.cor.Rdata",directory)
    wmesg <- sprintf("Warning, archived correlators are being read from: %s\nIf you changed a relevant option like 'skip' or the size of the ensemble, then you should specify 'read.cor=TRUE'!\n",archivename)
    cat(wmesg)
    load(archivename)
  }

  if(debug) {
    cat("Extracting observable\n")
  }
  meson.cor <- extract.obs(cmicor, vec.obs=observable,sign.vec=sign)

  # average LF and FL
  meson.cor <- avg.ls.cf(meson.cor)

  if(debug) {
    cat("Bootsrapping cf\n")
  }
  meson.cor <- bootstrap.cf(meson.cor,boot.R=boot.R,boot.l=boot.l,seed=seed)

  if(study.fitrange) {
    meson.cor.effectivemass <- bootstrap.effectivemass(meson.cor, type="solve", boot.R=boot.R, boot.l=boot.l, seed=seed)
    meson_2pt_study_fitrange(cf=meson.cor,effmass=meson.cor.effectivemass,name=directory,
                             kappa=kappa,useCov=useCov,q_masses=q_masses,boot.fit=F,
                             fps.disprel=fps.disprel,debug=debug)
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
  save.matrixfit <- computefps( save.matrixfit, mu1=q_masses$m1, mu2=q_masses$m2, Kappa=kappa, disprel=fps.disprel, boot.fit=T )

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
                matrixfit.ChiSqr=save.matrixfit$opt.res$value, matrixfit.dChiSqr=sd(save.matrixfit$opt.tsboot[4,]),
                matrixfit.dof=save.matrixfit$dof,
                effectivemass.ChiSqr=save.effectivemass$opt.res$value, effectivemass.dChiSqr=sd(save.effectivemass$massfit.tsboot[,2]),
                effectivemass.dof=save.effectivemass$dof )

  if(pause) {
    readline(prompt="Press any key to continue.")
  }

  # estimate some boundaries for the effective mass plot
  ymin <- save.effectivemass$opt.res$par[1]-20*sd(save.effectivemass$massfit.tsboot[,1])
  ymax <- save.effectivemass$opt.res$par[1]+20*sd(save.effectivemass$massfit.tsboot[,1])

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


  ## add another plot which shows the ratio between the fitted correlation function
  ## and the original data
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
