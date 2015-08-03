# to use this script, change this to the location where the helper scripts are located!
coderoot <- "~/code/R/misc_R_scripts"

# this function is used by xyplot to split the panels and give titles as given in strip.labels


source(paste(coderoot,"/plot_timeseries.R",sep=""))
source(paste(coderoot,"/plot_eigenvalue_timeseries.R",sep=""))

# convenience function for analyzing online data from a tmLQCD run
### the various parameters are used to build a directory name into which R descends to read
### pionline.dat and output.data
# addon can be used to add arbitrary text to the directory name (such as for replicas)
# plaquette and dH control whether these are plotted
# cg_col indicates which column in output.data should be used 

outputonline <- function(type,beta,L,T,t1,t2,skip,rundir,
  cg_col, evals, kappa=0, mul=0,
  csw=0,musigma=0,mudelta=0,muh=0,addon="",
  plaquette=TRUE, dH=TRUE,
  plotsize=5,debug=FALSE,trajlabel=FALSE,title=TRUE,
  pl=FALSE,method="uwerr",fit.routine="optim",oldnorm=F,S=3)
{
  navec <- c(val=NA,dval=NA,tauint=NA,dtauint=NA,Wopt=NA)
  result <- list(params=data.frame(L=L,T=T,kappa=kappa,csw=csw,mul=mul,muh=muh,musigma=musigma,mudelta=mudelta,N.online=0,N.plaq=0,skip=skip),
                 obs=data.frame(mpcac_fit=navec, mpcac_mc=navec, mpi=navec, fpi=navec, P=navec, dH=navec, expdH=navec, mineval=navec, maxeval=navec, CG.iter=navec))

  result$obs$mpcac_fit[1] <- 0.1
  
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

  filename <- sprintf("%s/piononline.dat",rundir)
  outfile <- sprintf("%s/output.data",rundir)
  pioncor <- try(readcmicor(filename))
  
  if(!any(class(pioncor)=='try-error')){
    if(debug){
      pion(pioncor,mu=mul,kappa=kappa,t1=t1,t2=t2,pl=TRUE,skip=skip,matrix.size=1)
    }

    onlineout <- onlinemeas(pioncor,t1=t1,t2=t2,kappa=kappa,mu=mul,skip=skip,method=method,pl=pl,fit.routine=fit.routine,oldnorm=oldnorm,S=S)

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
      trange=c(skip+1,(skip+length(onlineout$MChist.dpaopp))),
      pdf.filename=dpaopp_filename,
      ylab="$am_\\mathrm{PCAC}$",
      name="am_PCAC (MC history)",
      plotsize=plotsize,
      filelabel=filelabel,
      titletext=titletext,
      errorband_color=errorband_color,
      hist.by=0.0002)

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

    result$obs["val","mpcac_fit"] <- (onlineout$fitresult$par[3]*onlineout$fitresult$par[2]/onlineout$fitresult$par[1]/2.)
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

    result$obs$mpi[1] <- abs(onlineout$fitresultpp$par[2])
    result$obs$mpi[2] <- onlineout$uwerrresultmps$dvalue
    result$obs$mpi[3] <- onlineout$uwerrresultmps$tauint
    result$obs$mpi[4] <- onlineout$uwerrresultmps$dtauint
    result$obs$mpi[5] <- onlineout$uwerrresultmps$Wopt


    result$obs$fpi[1] <- 2*kappa*2*mul/sqrt(2)*abs(onlineout$fitresultpp$par[1])/sqrt(onlineout$fitresultpp$par[2]^3)
    result$obs$fpi[2] <- 2*kappa*2*mul/sqrt(2)*onlineout$uwerrresultfps$dvalue
    result$obs$fpi[3] <- onlineout$uwerrresultfps$tauint
    result$obs$fpi[4] <- onlineout$uwerrresultfps$dtauint
    result$obs$fpi[5] <- onlineout$uwerrresultfps$Wopt

    # something in the skip computation is odd, let's just solve it like this
    if(skip==0){
      shift <- 1
    } else {
      shift <- 0
    }
  } # if(!any(class(pioncor)=='try-error'))

  outdat <- NULL
  trange <- NULL
  if( plaquette || dH ) {
    # read output.data
    # determine maximum number of columns in output.data (when the mass preconditioning is changed,
    # the number of columns may change so we need to be able to deal with that)
    no_columns <- max(count.fields(outfile))
    outdat <- read.table(outfile,fill=TRUE,col.names=paste("V",1:no_columns,sep=""))
    trange <- c(skip+shift,length(outdat$V2))
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
      errorband_color=errorband_color,
      hist.by=0.00002)
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
      hist.xlim=c(-3,3),
      ylim=c(-2,2),
      hist.by=0.2)

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
      ylim=c(-2,4),
      hist.by=0.2)
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
      errorband_color=errorband_color,
      hist.by=5)
  }
  if( !missing("evals") ) {
    ev_pdf_filename <- sprintf("08_evals_%02d_%s", evals, filelabel )
    ev_filename <- sprintf("%s/monomial-%02d.data", rundir, evals )
    evaldata <- read.table(ev_filename)
    temp <- plot_eigenvalue_timeseries(dat=evaldata[trange[1]:trange[2],],
        trange=trange,
        pdf.filename = ev_pdf_filename,
        ylab = "eigenvalue",
        plotsize=plotsize,
        filelabel=filelabel,
        titletext=titletext,
        errorband_color=errorband_color )
    result$obs$mineval <- temp$mineval
    result$obs$maxeval <- temp$maxeval
  }

  print(result)

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
