# convenience function for analyzing online data from a tmLQCD run
### the various parameters are used to build a directory name into which R descends to read
### pionline.dat and output.data
### but a directory can also be provided via rundir with the understanding
### that L, T, kappa and mul must always be provided
# addon can be used to add arbitrary text to the directory name (such as for replicas)
# plaquette and dH control whether these are plotted
# cg_col indicates which column in output.data should be used 

# to use this script, change this to the location where the helper scripts are located!
coderoot <- "~/code/R/misc_R_scripts"

source(paste(coderoot,"/plot_timeseries.R",sep=""))
source(paste(coderoot,"/plot_eigenvalue_timeseries.R",sep=""))

outputonline <- function(L, T, t1, t2, kappa, mul,
  cg_col, evals, rundir, cg.ylim,
  type="", beta=0, csw=0, musigma=0, mudelta=0, muh=0, addon="",
  skip=0,
  plaquette=TRUE, dH=TRUE, acc=TRUE,
  plotsize=5,debug=FALSE,trajlabel=FALSE,title=FALSE,
  pl=FALSE,method="uwerr",fit.routine="optim",oldnorm=FALSE,S=1.5,
  omeas.start=0, omeas.stepsize=1, evals.stepsize=1)
{
  # store analysis results in practical R format, replacing entries as new data is added
  resultsfile <- "omeas.summary.RData"
  
  resultsum <- list()
  if(file.exists(resultsfile)){
    load(resultsfile)
  }

  # vector with NAs to initialise result data frame
  navec <- t(data.frame(val=NA,dval=NA,tauint=NA,dtauint=NA,Wopt=NA,stringsAsFactors=FALSE))
  
  # set up data structure for analysis results 
  result <- list(params=data.frame(L=L,T=T,t1=t1,t2=t2,kappa=kappa,csw=csw,mul=mul,muh=muh,
                                   musigma=musigma,mudelta=mudelta,N.online=0,N.plaq=0,skip=skip,stringsAsFactors=FALSE),
                 obs=data.frame(mpcac_fit=navec, 
                                mpcac_mc=navec, 
                                mpi=navec, 
                                fpi=navec, 
                                P=navec, 
                                dH=navec, 
                                expdH=navec, 
                                mineval=navec, 
                                maxeval=navec, 
                                CG.iter=navec, 
                                accrate=navec, stringsAsFactors=FALSE))
  
  errorband_color <- rgb(0.6,0.0,0.0,0.6)
  
  shift <- 0  
  # something in the skip computation is odd, let's just solve it like this
  if(skip==0){
    shift <- 1
  }
  
  if(missing(rundir)){
    rundir <- construct_rundir(type=type,beta=beta,L=L,T=T,kappa=kappa,mul=mul,
                             csw=csw,musigma=musigma,mudelta=mudelta,muh=muh,addon=addon,
                             debug=debug
                            )
  }
  
  filelabel <- rundir
      
  titletext <- NULL
  if(title) {
    titletext <- rundir
  } else {
    titletext <- ""
  }

  outfile <- sprintf("%s/output.data",rundir)
  
  # read online measurements
  omeas.files <- getorderedfilelist(path=rundir, basename="onlinemeas", last.digits=6)
  omeas.cnums <- getorderedconfignumbers(path=rundir, basename="onlinemeas", last.digits=6)
  omeas.idx   <- c((1+as.integer(ceiling(skip/omeas.stepsize))):length(omeas.files))
  omeas.files <- omeas.files[omeas.idx]
  omeas.cnums <- omeas.cnums[omeas.idx]
  pioncor <- readcmidatafiles( files=omeas.files, skip=0 )
  pioncor <- cbind( pioncor, rep(omeas.cnums, each=3*(T/2+1) ) )

  if(!any(class(pioncor)=='try-error')){
    #if(debug){
    #  pion(pioncor,mu=mul,kappa=kappa,t1=t1,t2=t2,pl=TRUE,skip=skip,matrix.size=1)
    #}

    # the correlation functions have been read externally, taking into account the measurement frequency
    # and possibly missing files. Therefore, skip=0! 
    onlineout <- onlinemeas(pioncor,t1=t1,t2=t2,kappa=kappa,mu=mul,skip=0,method=method,pl=pl,fit.routine=fit.routine,oldnorm=oldnorm,S=S)
    
    result$params$N.online <- onlineout$N

    if(debug){
      print(onlineout)
    }

    if(trajlabel){
      filelabel <- sprintf("%s_traj%d-%d",rundir,skip,(skip-1+length(onlineout$MChist.dpaopp)))
    } else {
      filelabel <- rundir
    }

    dpaopp_filename <- sprintf("01_dpaopp_%s",filelabel)
    result$obs$mpcac_mc <- plot_timeseries(dat=onlineout$MChist.dpaopp,
      trange=skip+c(0,(length(onlineout$MChist.dpaopp)-1)*omeas.stepsize),
      stepsize=omeas.stepsize,
      pdf.filename=dpaopp_filename,
      ylab="$am_\\mathrm{PCAC}$",
      name="am_PCAC (MC history)",
      plotsize=plotsize,
      filelabel=filelabel,
      titletext=titletext,
      errorband_color=errorband_color)
      #ist.by=0.0002))
    
    # adjust autocorrelation times to be in terms of trajectories
    result$obs$mpcac_mc[3] <- result$obs$mpcac_mc[3]*omeas.stepsize
    result$obs$mpcac_mc[4] <- result$obs$mpcac_mc[4]*omeas.stepsize
    result$obs$mpcac_mc[5] <- result$obs$mpcac_mc[5]*omeas.stepsize

    lengthdpaopp <- length(onlineout$MChist.dpaopp)
    mindpaopp <- min(onlineout$MChist.dpaopp)
    maxdpaopp <- max(onlineout$MChist.dpaopp)

    dpaopp_plateau_filename <- sprintf("02_dpaopp_plateau_%s",filelabel)
    tikzfiles <- tikz.init(basename=dpaopp_plateau_filename,width=plotsize,height=plotsize)
    op <- par(family="Palatino",cex.main=0.6,font.main=1)
    par(mgp=c(2,1,0))
    plotwitherror(x=onlineout$dpaopp$t,
      y=onlineout$dpaopp$mass,dy=onlineout$dpaopp$dmass,t='p',
      ylab="$am_\\mathrm{PCAC}$",
      xlab="$t/a$",
      main=titletext)
    rect(xleft=t1,
      xright=t2,
      ytop=onlineout$uwerrresultmpcac$value+onlineout$uwerrresultmpcac$dvalue,
      ybottom=onlineout$uwerrresultmpcac$value-onlineout$uwerrresultmpcac$dvalue,border=FALSE,col=errorband_color)
    abline(h=onlineout$uwerrresultmpcac$value,col="black")
    tikz.finalize(tikzfiles)

    result$obs$mpcac_fit <- t(data.frame(val=(onlineout$fitresult$par[3]*onlineout$fitresult$par[2]/onlineout$fitresult$par[1]/2.),
                                         dval=NA, tauint=NA, dtauint=NA, Wopt=NA, stringsAsFactors=FALSE))
    # no error or tauint from the fit

    mpi_plateau_filename <- sprintf("03_mpi_plateau_%s",filelabel)
    tikzfiles <- tikz.init(mpi_plateau_filename,width=plotsize,height=plotsize)
    op <- par(family="Palatino",cex.main=0.6,font.main=1)
    par(mgp=c(2,1,0))

    ploterror <- try(plotwitherror(x=onlineout$effmass$t,
      y=onlineout$effmass$m,dy=onlineout$effmass$dm,t='p',
      ylab="$aM_\\mathrm{PS}$",
      xlab="$t/a$",
      main=titletext),silent=FALSE)

    if(inherits(ploterror,"try-error")) {
      plot(x=onlineout$effmass$t,y=onlineout$effmass$m)
    }
    rect(xleft=t1,
      xright=t2,
      ytop=onlineout$uwerrresultmps$value+onlineout$uwerrresultmps$dvalue,
      ybottom=onlineout$uwerrresultmps$value-onlineout$uwerrresultmps$dvalue,border=FALSE,col=errorband_color)
    abline(h=onlineout$uwerrresultmps$value,col="black")
    tikz.finalize(tikzfiles)

    result$obs$mpi <- t(data.frame(val=abs(onlineout$fitresultpp$par[2]),
                                   dval=onlineout$uwerrresultmps$dvalue,
                                   tauint=onlineout$uwerrresultmps$tauint*omeas.stepsize,
                                   dtauint=onlineout$uwerrresultmps$dtauint*omeas.stepsize,
                                   Wopt=onlineout$uwerrresultmps$Wopt*omeas.stepsize, stringsAsFactors=FALSE) )

    result$obs$fpi <- t(data.frame(val=2*kappa*2*mul/sqrt(2)*abs(onlineout$fitresultpp$par[1])/sqrt(onlineout$fitresultpp$par[2]^3),
                                   dval=2*kappa*2*mul/sqrt(2)*onlineout$uwerrresultfps$dvalue,
                                   tauint=onlineout$uwerrresultfps$tauint*omeas.stepsize,
                                   dtauint=onlineout$uwerrresultfps$dtauint*omeas.stepsize,
                                   Wopt=onlineout$uwerrresultfps$Wopt*omeas.stepsize, stringsAsFactors=FALSE) )

    # something in the skip computation is odd, let's just solve it like this
    if(skip==0){
      shift <- 1
    } else {
      shift <- 0
    }
  } else { # if(!any(class(pioncor)=='try-error'))
    stop("outputonline: there was an error trying to read the output files (pionionline.dat or output.data)")
  }


  outdat <- NULL
  trange <- NULL
  if( plaquette || dH ) {
    # read output.data
    # determine maximum number of columns in output.data (when the mass preconditioning is changed,
    # the number of columns may change so we need to be able to deal with that)
    no_columns <- max(count.fields(outfile))
    outdat <- read.table(outfile,fill=TRUE,col.names=paste("V",1:no_columns,sep=""))
    trange <- c(skip+shift,length(outdat$V2))
    #result$accrate <- mean(outdat[trange[1]:trange[2],no_columns-2],na.rm=TRUE)
  }

  if(plaquette) {
    plaquette_filename <- sprintf("04_plaquette_%s",filelabel,title=filelabel)
    result$params$N.plaq <- trange[2]-trange[1]
    result$obs$P <- plot_timeseries(dat=outdat$V2[trange[1]:trange[2]],
      trange=trange,
      pdf.filename=plaquette_filename,
      ylab="$ \\langle P \\rangle$" ,
      name="plaquette",
      plotsize=plotsize,
      filelabel=filelabel,
      titletext=titletext,
      errorband_color=errorband_color)
      #ist.by=0.00002))
  }
  if(dH) {
    dH_filename <- sprintf("05_dH_%s",filelabel)
    result$obs$dH <- plot_timeseries(dat=outdat$V3[trange[1]:trange[2]],
      trange=trange,
      pdf.filename=dH_filename,
      ylab="$ \\delta H $",
      name="dH",
      plotsize=plotsize,
      filelabel=filelabel,
      titletext=titletext,
      errorband_color=errorband_color,
      ylim=c(-2,3))
      #ist.by=0.2))

    expdH_filename <- sprintf("06_expdH_%s",filelabel)
    result$obs$expdH <- plot_timeseries(dat=outdat$V4[trange[1]:trange[2]],
      trange=trange,
      pdf.filename=expdH_filename,
      ylab="$ \\exp(-\\delta H) $",
      name="exp(-dH)",
      plotsize=plotsize,
      filelabel=filelabel,
      titletext=titletext,
      errorband_color=errorband_color,
      hist.xlim=c(-2,4),
      ylim=c(-0,6))
      #ist.by=0.2))
  }
  if( !missing("cg_col") ) {
    cg_filename <- sprintf("07_cg_iter_%s", filelabel)
    result$obs$CG.iter <- plot_timeseries(dat=outdat[trange[1]:trange[2],cg_col],
      trange=trange,
      pdf.filename=cg_filename,
      ylab="$N^\\mathrm{iter}_\\mathrm{CG}$",
      name="CG iterations",
      plotsize=plotsize,
      filelabel=filelabel,
      titletext=titletext,
      errorband_color=errorband_color)
      #ist.by=5))
  }
  # TODO: if eigenvalue measurement is not every trajectory, need to correctly handle 
  # evals.stepsize and trange!!!
  if( !missing("evals") ) {
    ev_pdf_filename <- sprintf("08_evals_%02d_%s", evals, filelabel )
    ev_filename <- sprintf("%s/monomial-%02d.data", rundir, evals )
    
    evaldata <- tryCatch(read.table(ev_filename,stringsAsFactors=FALSE), 
                         error=function(e){ stop(sprintf("Reading of %s failed!",ev_filename)) } )

    temp <- plot_eigenvalue_timeseries(dat=evaldata[trange[1]:trange[2],],
        trange=trange,
        stepsize=evals.stepsize,
        pdf.filename = ev_pdf_filename,
        ylab = "eigenvalue",
        plotsize=plotsize,
        filelabel=filelabel,
        titletext=titletext,
        errorband_color=errorband_color )
    result$obs$mineval <- temp$mineval
    result$obs$maxeval <- temp$maxeval
  }
  if( acc == TRUE ){
    # finally add acceptance rate
    accrate_filename <- sprintf("09_accrate_%s",filelabel,title=filelabel)
    result$obs$accrate <- plot_timeseries(dat=outdat[trange[1]:trange[2],no_columns-2],
      trange=trange,
      pdf.filename=accrate_filename,
      ylab="$ \\langle P_\\mathrm{acc} \\rangle$" ,
      name="accrate",
      plotsize=plotsize,
      filelabel=filelabel,
      titletext=titletext,
      errorband_color=errorband_color,
      hist.by=0.5)
  }

  print(result$params)
  print(t(result$obs))

  # if the script "pdfcat" exists, concatenate the plots and remove the individual files
  if( Sys.which("pdfcat") != "" ) {
    commands <- c(sprintf("pdfcat analysis_%s.pdf 0?*_%s.pdf",filelabel,filelabel),
                  sprintf("rm -f 0?_*_%s.pdf",filelabel) )
    for( command in commands ) {
      print(paste("calling",command))
      system(command=command)
    } 
  } else {
    print("pdfcat not found, not concatenating plots!")
  }
  if( Sys.which("pdfcrop") != "" ) {
    commands <- c(sprintf("pdfcrop analysis_%s.pdf analysis_%s.pdf",filelabel,filelabel))
    for( command in commands ) {
      print(paste("calling",command))
      system(command=command)
    } 
  } else {
    print("pdfcrop not found, not cropping plots!")
  }

  resultsum[[rundir]] <- result
  save(resultsum,file=resultsfile)

  return(invisible(result))
}

construct_rundir <- function(type,beta,L,T,kappa=0,mul=0,csw=0,musigma=0,mudelta=0,muh=0,addon="",debug=FALSE) {
  rundir <- NULL
  rundir <- sprintf("%s_b%s-L%dT%d",type,beta,L,T)

  if(csw!=0) {
    rundir <- sprintf("%s-csw%s",rundir,csw)
  }

  # mul=0 && kappa=0 means pure gauge, which is perhaps a bit misleading...
  if( mul != 0 && kappa != 0 ) {
    # silly sprintf prints numbers smaller than 0.001 in scientific notation unless g is used
    # on the other hand, for larger numbers, trailing zeroes are added...
    # so we use %g only in the former case and pass mu as a string otherwise!
    if(mul >= 1e-3) {
      rundir <- sprintf("%s-k%s-mul%s",rundir,kappa,mul)
    } else {
      rundir <- sprintf("%s-k%s-mul%g",rundir,kappa,mul)
    }
  }

  if(muh != 0) {
    rundir <- sprintf("%s-muh%s",rundir,muh)
  }

  if(musigma!=0){
    rundir <- sprintf("%s-musigma%s",rundir,musigma)
  }
  if(mudelta!=0){
    rundir <- sprintf("%s-mudelta%s",rundir,mudelta)
  }
  if(addon!=""){
    rundir <- sprintf("%s_%s",rundir,addon)
  }

  if(debug) {
    cat("Trying to read from directory:", rundir,"\n")
  }

  rundir
}
