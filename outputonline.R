# to use this script, change this to the location where the helper scripts are located!
coderoot <- "~/code/R/misc_R_scripts"

# this function is used by xyplot to split the panels and give titles as given in strip.labels

my.strip <- function(which.given, which.panel, ...) {
  strip.labels <- 
    c(expression(am[PCAC]),"<P>")
    panel.rect(0, 0, 1, 1, col="#aabbff", border=1)
    panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,
    lab=strip.labels[which.panel[which.given]])
}

source(paste(coderoot,"/plot_timeseries.R",sep=""))
source(paste(coderoot,"/plot_eigenvalue_timeseries.R",sep=""))

# convenience function for analyzing online data from a tmLQCD run
### the various parameters are used to build a directory name into which R descends to read
### pionline.dat and output.data
# oneplot=TRUE will produce one lattice plot with CG iterations, plaquette and mpcac histories
### otherwise three separate plots are produced
# addon can be used to add arbitrary text to the directory name (such as for replicas)
# plaquette and dH control whether these are plotted
# cg_col indicates which column in output.data should be used 

outputonline <- function(type,beta,L,T,kappa,mul,t1,t2,skip,
  csw=0,musigma=0,mudelta=0,muh=0,addon="",plaquette=TRUE,
  dH=TRUE,oneplot=FALSE,plotsize=5,debug=FALSE,trajlabel=FALSE,
  title=TRUE,pl=FALSE,method="uwerr",fit.routine="optim",oldnorm=F,cg_col,
  evals)
{
  errorband_color <- rgb(0.6,0.0,0.0,0.6)
  rundir <- NULL
  rundir <- sprintf("%s_b%s-L%dT%d",type,beta,L,T)

  if(csw!=0) {
    rundir <- sprintf("%s-csw%s",rundir,csw)
  }

  # silly sprintf prints numbers smaller than 0.001 in scientific notation unless g is used
  # on the other hand, for larger numbers, trailing zeroes are added...
  # so we use %g only in the former case and pass mu as a string otherwise!
  if(mul >= 1e-3) {
    rundir <- sprintf("%s-k%s-mul%s",rundir,kappa,mul)
  } else {
    rundir <- sprintf("%s-k%s-mul%g",rundir,kappa,mul)
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
  
  titletext <- NULL
  if(title) {
    titletext <- rundir
  } else {
    titletext <- ""
  }

  filename <- sprintf("%s/piononline.dat",rundir)
  outfile <- sprintf("%s/output.data",rundir)
  pioncor <- readcmicor(filename)

  if(debug){
    pion(pioncor,mu=mul,kappa=kappa,t1=t1,t2=t2,pl=TRUE,skip=skip,matrix.size=1)
  }

  onlineout <- onlinemeas(pioncor,t1=t1,t2=t2,kappa=kappa,mu=mul,skip=skip,method=method,pl=pl,fit.routine=fit.routine,oldnorm=oldnorm,S=10)

  print(onlineout)

  filelabel <- NULL
  if(trajlabel){
    filelabel <- sprintf("%s_traj%d-%d",rundir,skip,(skip-1+length(onlineout$MChist.dpaopp)))
  } else {
    filelabel <- rundir
  }

  # we are plotting histories of the plaquette and the PCAC mass 
  # if we want to make one plot we will use the "lattice" package, in particular xyplot
  # in order to feed it with data we need to make a data frame which will
  if(oneplot) {
    print("making one plot with two rows")
    # these are the upper boundaries of the PCAC history and contents of output.data, respectively
    # they are different because skip is passed directly to onlinemeas above while the whole of output.data is read 
    l_max <- length(onlineout$MChist.dpaopp) 
    l_skip_max <- (skip-1+length(onlineout$MChist.dpaopp))

    # read output.data
    # determine maximum number of columns in output.data (when the mass preconditioning is changed,
    # the number of columns may change) 
    no_columns <- max(count.fields(outfile))
    outdat <- as.matrix(read.table(outfile,fill=TRUE,col.names=1:no_columns))

    # create a data frame for xyplot
    combined <- data.frame( 
      t=seq(skip,l_skip_max),
      mpcac=onlineout$MChist.dpaopp[0:l_max],
      plaq=outdat[skip:l_skip_max,2])
    
    summary(combined)
    print(length(combined))

    algo_filename <- sprintf("algo_%s.pdf",filelabel)

    require(lattice)
    pdf(algo_filename,family="Palatino",height=3,width=9,title=filelabel)

    # xyplot is a really annoying function to use! The plot here is based on 
    # http://www.sr.bham.ac.uk/~ajrs/R/gallery/plot_midday_weather_profiles.txt
    # which produces this: http://www.sr.bham.ac.uk/~ajrs/R/gallery/midday_weather_profiles.png

    # the first argument is to be understood as
    # mpcac as a function of t, plaq as a function of t and cg as a function of t
    ### in order: bottom to top (for some insane reason) 
    # the "outer" argument causes the plots to be done in different panels, as described by layout
    # the my.strip function gives the correct titles for the panels and draws the grids etc. in REVERSE ORDER
    pl.obj <- xyplot(mpcac + plaq ~ t,
      data=combined,outer=T,ylab="",xlab=expression(t[HMC]),
      scales=list(y="free",rot=0),
      strip=my.strip,
      layout=c(1,2),
      # here the actual plotting is done for each panel
      # the function is executed line by line just like a normal function elsewhere
      # so for each data set y passed, xyplot, abline and rect will be executed
      panel=function(x, y, ...) { 
        panel.grid(h=-1, v=0) # start by plotting gridlines
        uw <- uwerrprimary(y) # do error analysis for the given data
        panel.xyplot(x, y, ..., type="l", lwd=0.5,col="black")  # plot the raw data
        panel.abline(h=uw$value, lty=2, lwd=1) # expectation value from uwerrprimary
        # rect is used to draw a semi-transparent error band
        panel.rect(xleft=(skip-100),
          xright=(skip+length(y)+100),
          ytop=(uw$value+uw$dvalue),
          ybottom=(uw$value-uw$dvalue),
          border=FALSE,col=errorband_color) } )
    # finally, we need to actually draw the plot
    plot(pl.obj)
    dev.off()


  } else {
    dpaopp_filename <- sprintf("01_dpaopp_%s.pdf",filelabel)
    
    plot_timeseries(dat=onlineout$MChist.dpaopp,
      trange=c(skip+1,(skip+length(onlineout$MChist.dpaopp))),
      pdf.filename=dpaopp_filename,
      ylab=expression(am[PCAC]),
      name="am_PCAC (MC history)",
      plotsize=plotsize,
      filelabel=filelabel,
      titletext=titletext,
      errorband_color=errorband_color)    
  
    lengthdpaopp <- length(onlineout$MChist.dpaopp)
    mindpaopp <- min(onlineout$MChist.dpaopp)
    maxdpaopp <- max(onlineout$MChist.dpaopp)
    
    dpaopp_plateau_filename <- sprintf("02_dpaopp_plateau_%s.pdf",filelabel)
    pdf(dpaopp_plateau_filename,width=plotsize,height=plotsize,title=filelabel)
    op <- par(family="Palatino",cex.main=0.6,font.main=1)
    par(mgp=c(2,1,0))
    plotwitherror(x=onlineout$dpaopp$t,
      y=onlineout$dpaopp$mass,dy=onlineout$dpaopp$dmass,t='p',
      ylab=expression(am[PCAC]),
      xlab=expression(t),
      main=titletext)
    rect(xleft=t1,
      xright=t2,
      ytop=onlineout$uwerrresultmpcac$value+onlineout$uwerrresultmpcac$dvalue,
      ybottom=onlineout$uwerrresultmpcac$value-onlineout$uwerrresultmpcac$dvalue,border=FALSE,col=errorband_color)
    abline(h=onlineout$uwerrresultmpcac$value,col="black")
    dev.off()
    
    mpi_plateau_filename <- sprintf("03_mpi_plateau_%s.pdf",filelabel)
    pdf(mpi_plateau_filename,width=plotsize,height=plotsize,title=filelabel)
    op <- par(family="Palatino",cex.main=0.6,font.main=1)
    par(mgp=c(2,1,0))

    ploterror <- try(plotwitherror(x=onlineout$effmass$t,
      y=onlineout$effmass$m,dy=onlineout$effmass$dm,t='p',
      ylab=expression(am[PS]),
      xlab=expression(t),
      main=titletext),silent=FALSE)

    if(inherits(ploterror,"try-error")) { 
      plot(x=onlineout$effmass$t,y=onlineout$effmass$m)
    }
    rect(xleft=t1,
      xright=t2,
      ytop=onlineout$uwerrresultmps$value+onlineout$uwerrresultmps$dvalue,
      ybottom=onlineout$uwerrresultmps$value-onlineout$uwerrresultmps$dvalue,border=FALSE,col=errorband_color)
    abline(h=onlineout$uwerrresultmps$value,col="black")
    dev.off()
    
    # something in the skip computation is odd, let's just solve it like this
    if(skip==0){
      shift <- 1
    } else {
      shift <- 0
    }
    
    outdat <- NULL
    trange <- NULL
    if( plaquette || dH ) {
      # read output.data
      # determine maximum number of columns in output.data (when the mass preconditioning is changed,
      # the number of columns may change) 
      no_columns <- max(count.fields(outfile))
      outdat <- read.table(outfile,fill=TRUE,col.names=paste("V",1:no_columns,sep=""))
      trange <- c(skip+shift,length(outdat$V2))
    }

    if(plaquette) {
      plaquette_filename <- sprintf("04_plaquette_%s.pdf",filelabel,title=filelabel)
      plot_timeseries(dat=outdat$V2[trange[1]:trange[2]],
        trange=trange,
        pdf.filename=plaquette_filename,
        ylab=expression("<P>"),
        name="plaquette",
        plotsize=plotsize,
        filelabel=filelabel,
        titletext=titletext,
        errorband_color=errorband_color)
    }
    if(dH) {
      dH_filename <- sprintf("05_dH_%s.pdf",filelabel)
      plot_timeseries(dat=outdat$V3[trange[1]:trange[2]],
        trange=trange,
        pdf.filename=dH_filename,
        ylab=expression(paste(delta,"H",sep="")),
        name="dH",
        plotsize=plotsize,
        filelabel=filelabel,
        titletext=titletext,
        errorband_color=errorband_color)
        
      expdH_filename <- sprintf("06_expdH_%s.pdf",filelabel)  
      plot_timeseries(dat=outdat$V4[trange[1]:trange[2]],
        trange=trange,
        pdf.filename=expdH_filename,
        ylab=expression(paste(paste("exp(-",delta),"H)")),
        name="expdH",
        plotsize=plotsize,
        filelabel=filelabel,
        titletext=titletext,
        errorband_color=errorband_color)
    }
    if( !missing("cg_col") ) {
      cg_filename <- sprintf("07_cg_iter_%s.pdf", filelabel)
      plot_timeseries(dat=outdat[trange[1]:trange[2],cg_col],
        trange=trange,
        pdf.filename=cg_filename,
        ylab="CG iterations",
        name="CG iterations",
        plotsize=plotsize,
        filelabel=filelabel,
        titletext=titletext,
        errorband_color=errorband_color)
    }
    if( !missing("evals") ) {
      ev_pdf_filename <- sprintf("08_evals_%02d_%s.pdf", evals, filelabel )
      ev_filename <- sprintf("%s/monomial-%02d.data", rundir, evals )
      evaldata <- read.table(ev_filename)
      plot_eigenvalue_timeseries(dat=evaldata[trange[1]:trange[2],],
         trange=trange,
         pdf.filename = ev_pdf_filename,
         ylab = "eigenvalue",
         plotsize=plotsize,
         filelabel=filelabel,
         titletext=titletext,
         errorband_color=errorband_color )
    }
  }

  if( Sys.which("pdfcat") != "" ) {
    command <- sprintf("pdfcat analysis_%s.pdf 0?*_%s.pdf",filelabel,filelabel)
    print(paste("calling",command))
    system(command=command)
  } else {
    print("pdfcat not found, not concatenating plots!")
  }
}
