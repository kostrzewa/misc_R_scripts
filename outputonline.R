
# this function is used by xyplot to split the panels and give titles as given in strip.labels

my.strip <- function(which.given, which.panel, ...) {
  strip.labels <- 
    c("CG iterations (lightest)",expression(am[PCAC]),"<P>")
    panel.rect(0, 0, 1, 1, col="#aabbff", border=1)
    panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,
    lab=strip.labels[which.panel[which.given]])
}

# convenience function for analyzing online data from a tmLQCD run
### the various parameters are used to build a directory name into which R descends to read
### pionline.dat and output.data
# oneplot=TRUE will produce one lattice plot with CG iterations, plaquette and mpcac histories
### otherwise three separate plots are produced
# addon can be used to add arbitrary text to the directory name (such as for replicas)
# plaquette and dH control whether these are plotted
# cg_col indicates which column in output.data should be used 

outputonline <- function(type,beta,L,T,kappa,mu,t1,t2,skip,csw=0,musigma=0,mudelta=0,addon="",cg_col,plaquette=TRUE,dH=TRUE,oneplot=FALSE,plotsize=5,debug=FALSE,trajlabel=FALSE)
{
  errorband_color <- rgb(0.6,0.0,0.0,0.6)
  rundir <- NULL
  rundir <- sprintf("%s_b%g-L%dT%d",type,beta,L,T)

  if(csw!=0) {
    rundir <- sprintf("%s-csw%g",rundir,csw)
  }

  rundir <- sprintf("%s-k%g-mu%g",rundir,kappa,mu)

  if(musigma!=0){
    rundir <- sprintf("%s-musigma%g",rundir,musigma)
  }
  if(mudelta!=0){
    rundir <- sprintf("%s-mudelta%g",rundir,mudelta)
  }
  if(addon!=""){
    rundir <- sprintf("%s_%s",rundir,addon)
  }

  filename <- sprintf("%s/piononline.dat",rundir)
  outfile <- sprintf("%s/output.data",rundir)
  pioncor <- readcmicor(filename)

  if(debug){
    pion(pioncor,mu=mu,kappa=kappa,t1=t1,t2=t2,pl=TRUE,skip=skip,matrix.size=1)
  }

  onlineout <- onlinemeas(pioncor,t1=t1,t2=t2,kappa=kappa,mu=mu,skip=skip,method="uwerr")

  print(onlineout)

  filelabel <- NULL
  if(trajlabel){
    filelabel <- sprintf("%s_traj%d-%d",rundir,skip,(skip-1+length(onlineout$MChist.dpaopp)))
  } else {
    filelabel <- sprintf("%s",rundir) 
  }

  # we are plotting histories of the plaquette, PCAC mass and the number of acceptance CG iterations 
  # if we want to make one plot we will use the "lattice" package, in particular xyplot
  # in order to feed it with data we need to make a data frame which will
  if(oneplot) {
    print("making one plot with three rows")
    # these are the upper boundaries of the PCAC history and contents of output.data, respectively
    # they are different because skip is passed directly to onlinemeas above while the whole of output.data is read 
    l_max <- length(onlineout$MChist.dpaopp) 
    l_skip_max <- (skip-1+length(onlineout$MChist.dpaopp))

    # read output.data
    outdat <- as.matrix(read.table(outfile))

    # create a data frame for xyplot
    combined <- data.frame( 
      t=seq(skip,l_skip_max),
      mpcac=onlineout$MChist.dpaopp[0:l_max],
      plaq=outdat[skip:l_skip_max,2],
      cg=outdat[skip:l_skip_max,cg_col])
    
    summary(combined)
    print(length(combined))

    algo_filename <- sprintf("algo_%s.pdf",filelabel)

    require(lattice)
    pdf(algo_filename,family="Palatino",height=5,width=8)

    # xyplot is a really annoying function to use! The plot here is based on 
    # http://www.sr.bham.ac.uk/~ajrs/R/gallery/plot_midday_weather_profiles.txt
    # which produces this: http://www.sr.bham.ac.uk/~ajrs/R/gallery/midday_weather_profiles.png

    # the first argument is to be understood as
    # mpcac as a function of t, plaq as a function of t and cg as a function of t
    ### in order: bottom to top (for some insane reason) 
    # the "outer" argument causes the plots to be done in different panels, as described by layout
    # the my.strip function gives the correct titles for the panels and draws the grids etc. in REVERSE ORDER
    pl.obj <- xyplot(cg + mpcac + plaq ~ t,
      data=combined,outer=T,ylab="",xlab=expression(t[HMC]),
      scales=list(y="free",rot=0),
      strip=my.strip,
      layout=c(1,3),
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
    dpaopp_filename <- sprintf("dpaopp_%s.pdf",filelabel)
    pdf(dpaopp_filename,width=plotsize,height=plotsize)
    op <- par(family="Palatino")
    plot(x=seq(skip+1,(skip+length(onlineout$MChist.dpaopp)),1),
      onlineout$MChist.dpaopp,t='l',
      ylab=expression(am[PCAC]),
      xlab=expression(t[HMC]),
      main=rundir)
    
    rect(xleft=(skip-50),
      xright=(skip+length(onlineout$MChist.dpaopp)+50),
      ytop=onlineout$uwerrresultmpcac$value+onlineout$uwerrresultmpcac$dvalue,
      ybottom=onlineout$uwerrresultmpcac$value-onlineout$uwerrresultmpcac$dvalue,border=FALSE,col=errorband_color)
    
    abline(h=onlineout$uwerrresultmpcac$value,col="black")
    
    lengthdpaopp <- length(onlineout$MChist.dpaopp)
    mindpaopp <- min(onlineout$MChist.dpaopp)
    maxdpaopp <- max(onlineout$MChist.dpaopp)
    dev.off()
    
    # something in the skip computation is odd, let's just solve it like this
    if(skip==0){
      shift <- 1
    } else {
      shift <- 0
    }
    
    if(plaquette) {
      outdat <- read.table(outfile)
      plaquette_filename <- sprintf("plaquette_%s.pdf",filelabel)
      pdf(plaquette_filename,width=plotsize,height=plotsize)
      op <- par(family="Palatino")
      plot(x=seq(skip+shift,length(outdat$V2),1),outdat$V2[skip:length(outdat$V2)],t='l',ylab=expression("<P>"),xlab=expression(t[HMC]),main=rundir)
      dev.off()
    }
    if(dH) {
      dH_filename <- sprintf("dH_%s.pdf",filelabel)
      pdf(dH_filename,width=plotsize,height=plotsize)
      op <- par(family="Palatino")
      plot(x=seq(skip+shift,length(outdat$V3),1),outdat$V3[skip:length(outdat$V3)],t='l',ylab=expression(paste(delta,"H")),xlab=expression(t[HMC]),main=rundir) 
      dev.off()

      expdH_filename <- sprintf("expdH_%s.pdf",filelabel)
      pdf(expdH_filename,width=plotsize,height=plotsize)
      op <- par(family="Palatino")
      plot(x=seq(skip+shift,length(outdat$V4),1),outdat$V4[skip:length(outdat$V4)],t='l',ylab=expression(paste(paste("exp(-",delta),"H)")),xlab=expression(t[HMC]),main=rundir) 
      dev.off()
      
      uw.dH <- uwerrprimary(outdat$V3[skip:length(outdat$V3)])
      uw.expdH <- uwerrprimary(outdat$V4[skip:length(outdat$V4)])
      print("uw.dH")
      summary(uw.dH)
      print("uw.exp(-dH)")
      summary(uw.expdH)
    }
  }
}
