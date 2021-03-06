# this function reads and analyses the output files from the ndcloverrat_num_deriv branch of kostrzewa/tmLQCD
# which writes out the numerical and analytical derivatives for a given monomial

# those files contain the acceptance precision, the force precision and the rotation angle 
# for the numerical derivative in the first three doubles and the derivative itself for the
# remainder of the file

# dir      - subdirectory from which the data files should be read
# single   - the arguments filename and pattern are only considered if single is TRUE
# filename - one specific filename to be read (optional)
# pattern  - one specific pattern to be prepended (see below) (optional)
# indices  - for the non-overview plots, this selects one/multiple particular lattice point(s), direction(s) and generator(s) (optional)
#            otherwise all indices are used         

analysis_num_deriv <- function(dir,filename,pattern,indices,numerical=TRUE,volume=2^4,int.steps=101,trajs=1,tau=1,endian="little",single=FALSE,all=FALSE,width=4,height=4,text=FALSE,finite.width=6,finite.fmax=60) {
 
  # we will use list.files with a pattern to read all the data files of interest
  # in the next few lines we construct which files we would like to read 

  # there are essentially two usecases:
  # 1) we have a list of files of the form
  #      epsX_precY_f_*.bin
  #    where epsX specified the rotation angle and precY the acceptance and force precision
  #   
  # 2) we have one pair of files which begin with the name of the monomial
  #      monomial_f_*.bin
  n_pat <- "*_f_numerical.bin"
  a_pat <- "*_f_analytical.bin"

  # we can also read a single pair of files, which is what happens when single is TRUE  
  if(single) {
    # either the default filenames
    a_pat <- "f_analytical.bin"
    n_pat <- "f_numerical.bin"
    # or ones prepended with "pattern"
    if(!missing(pattern)) {
      a_pat <- sprintf("%s_%s",pattern,a_pat)
      n_pat <- sprintf("%s_%s",pattern,n_pat)
    }
    # or an explicit filename, in which case the same file will be read for
    # analytical and numerical data, which is not necessarily desirable 
    if(!missing(filename)) {
      a_pat <- n_pat <- filename
    }
  }
  
  n_files <- list.files(dir,pattern=n_pat,full.names=FALSE)
  a_files <- list.files(dir,pattern=a_pat,full.names=FALSE)
  
  # these lists will store the analytical and numerical derivative 
  df_a <- list()
  df_n <- list()

  # begin by reading the numerical data
  for( fname in n_files ) {
    n_file <- file(fname,"rb")
    # extract the acceptance and force precisions as well as the rotation angle
    # for the numerical derivative 
    ap <- readBin(n_file,double(),n=1,endian=endian)
    fp <- readBin(n_file,double(),n=1,endian=endian)
    eps <- readBin(n_file,double(),n=1,endian=endian)
    
    # the derivative consists of 8 real numbers per lattice direction per lattice point and there is one entry per integration step
    # we store it as an array with int.steps*trajs rows and volume*4*8 columns (note the transposition) 
    df_n[[length(df_n)+1]] <- list(ap=ap,fp=fp,eps=eps,df=t(array(data=readBin(n_file,double(),n=trajs*int.steps*volume*4*8,endian=endian),dim=c(volume*4*8,int.steps*trajs))))

    close(n_file)
  }

  # we want to load only unique force precisions, these are accumulated here (see below)
  precs <- c()
  # now we read the analytical data
  for( fname in a_files ) {
    cat("Reading ", fname, "\n")
    a_file <- file(fname,"rb")
    ap <- readBin(a_file,double(),n=1,endian=endian)
    fp <- readBin(a_file,double(),n=1,endian=endian)
    eps <- readBin(a_file,double(),n=1,endian=endian)
    
    # from what has been described above, when scanning across the precisions, the analytical derivative
    # may be present far too many times because there is one for each epsilon value, here we 
    # precs keeps track of which one we have already read with the fuzzy condition below
    if(!any( abs(precs-fp) < 10^-30)) {
      df_a[[length(df_a)+1]] <- list(ap=ap,fp=fp,eps=eps,df=t(array(data=readBin(a_file,double(),n=trajs*int.steps*volume*4*8,endian=endian),dim=c(volume*4*8,int.steps*trajs))))
      precs <- c(precs,fp)
    }

    close(a_file)
  }


  # we start by plotting the analytical force for all lattice points, all directions and all generators
  if(all){
    pdf.filename <- "df_a.3D.all.pdf"
    if(!missing(single) && !missing(pattern)) {
      pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
    }
    pdf(file=pdf.filename,width=width,height=height,onefile=T)
    library(plot3D)
    rows <- nrow(df_a[[length(df_a)]]$df)
    cols <- ncol(df_a[[length(df_a)]]$df)
    xcoords <- rep(times=cols,seq(1,rows))
    ycoords <- c()
    for( i in 1:cols ) {
      idcs <- seq((i-1)*rows,i*rows)
      ycoords[idcs] <- i
    }
    # first in 3D as one plot
    points3D(x=xcoords,
            y=ycoords,
            z=as.vector(df_a[[length(df_a)]]$df[1:rows,1:cols]),pch=".",
            xlab="t",ylab="idx",zlab="dP")
    dev.off()
    
    tex.basename <- "df_a.all"
    if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s",pattern,tex.basename)
    }
    tikzfiles <- tikz.init(basename=tex.basename,width=width,height=height)
    # and then we iterate over the columns (which correspond to generators, directions and lattice points,
    # with the rightmost iterating slowest)
    for( i in 1:ncol(df_a[[length(df_a)]]$df) ) {
      plot(y=df_a[[length(df_a)]]$df[,i],
           x=seq(1,length(df_a[[length(df_a)]]$df[,i]))*tau/int.steps,
           t='l',lwd=3,
           xlab="t",ylab="dP")
      # we do a fast fourier transform of the trajectory which provides us with a few of the dominant frequenc(y/ies)
      x <- df_a[[length(df_a)]]$df[,i]
      ft <- Mod(fft(x-mean(x)))/sqrt(length(x))
      ft <- ft[2:length(ft)]
      plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1,length(ft)/2)/(tau),log='x',
           xlab="f",ylab="$Re(\\hat{dP}(f))^2$",lwd=3 )
    }
    tikz.finalize(tikzfiles)
  
    if(numerical){
      # and the same game for the numerical derivative
      pdf.filename <- "df_n.3D.all.pdf"
      if(!missing(single) && !missing(pattern)) {
        pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
      }
      pdf(file=pdf.filename,width=width,height=height,onefile=T)
      library(plot3D)
      rows <- nrow(df_n[[length(df_n)]]$df)
      cols <- ncol(df_n[[length(df_n)]]$df)
      xcoords <- rep(times=cols,seq(1,rows))
      ycoords <- c()
      for( i in 1:cols ) {
        idcs <- seq((i-1)*rows,i*rows)
        ycoords[idcs] <- i
      }
      points3D(x=xcoords,
              y=ycoords,
              z=as.vector(df_n[[length(df_n)]]$df[1:rows,1:cols]),pch=".",
              xlab="t",ylab="idx",zlab="dP")
  
      tex.basename <- "df_n.all"
      if(!missing(single) && !missing(pattern)) {
        tex.basename <- sprintf("%s.%s",pattern,tex.basename)
      }
      tikzfiles <- tikz.init(basename=tex.basename,width=width,height=height)
      for( i in 1:ncol(df_n[[length(df_n)]]$df) ) {
        plot(y=df_n[[length(df_n)]]$df[,i],
             x=seq(1,length(df_n[[length(df_n)]]$df[,i]))*tau/int.steps,
             t='l',xlab="t",ylab="dP",lwd=3)
        x <- df_n[[length(df_n)]]$df[,i]
        ft <- Mod(fft(x-mean(x)))/sqrt(length(x))
        ft <- ft[2:length(ft)]
        plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1,length(ft)/2)/(tau),log='x',
        xlab="f",ylab="$Re(\\hat{dP}(f))^2$",lwd=3 )
      }
      tikz.finalize(tikzfiles)
    }
  } else {
    tex.basename <- "df_a"
    if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s.%s",pattern,int.steps,tex.basename)
    }
    tikzfiles <- tikz.init(basename=tex.basename,width=width,height=height)
    par(mgp=c(2,0.2,0))
    tck <- 0.2/height
    irange <- 1:ncol(df_a[[length(df_a)]]$df)
    if(!missing(indices)){
      irange <- indices
    }
    for( i in irange ) {
      #op <- par(mar = c(4,5,3,1) + 0.1,mgp = c(2,0.3,0) )
      plot(y=df_a[[length(df_a)]]$df[,i],
           x=seq(1,length(df_a[[length(df_a)]]$df[,i]))*tau/int.steps,
           type='l',pch='.',lwd=3,tck=tck,las=1,xaxt='n',
           xlab="",ylab="$\\delta P^a_\\mu(x,\\tau)$")
      mtext(side=1,line=1.0,text="$\\tau$")
      if(length(df_a[[length(df_a)]]$df[,i])*tau/int.steps > 5){
        axis(side=1, labels=TRUE, at=seq(0,20,5),tck=tck)
        axis(side=1, labels=FALSE, at=outer(c(1:4),seq(0,20,5), FUN="+"), tck=0.6*tck)
      } else {
        axis(side=1,tck=tck)
      }

      # we do a fast fourier transform of the trajectory which provides us with a few of the dominant frequenc(y/ies) 
      x <- df_a[[length(df_a)]]$df[,i]
      ft <- Mod(fft(x-mean(x)))/sqrt(length(x))
      ft <- ft[2:length(ft)]
      incr <- 20
      if(length(ft)<4000)
        incr <- incr/2
      if(length(ft)<2000)
        incr <- incr/2
      plot(y=ft,type='l',x=seq(1,length(ft))/(tau),xlim=c(1,length(ft)/2)/(tau),log='y',lwd=3,las=1,xaxt='n',yaxt='n',
           xlab="",ylab="",tck=tck)
      axis(side=1, labels=TRUE, at=seq(-20,500,incr),tck=tck)
      if(incr>=10)
        axis(side=1, labels=FALSE, at=outer(c(2,4,6,8),seq(-10,500,10), FUN="+"), tck=0.6*tck)
      else
        axis(side=1, labels=FALSE, at=outer(c(1:4,6:9),seq(-10,500,10), FUN="+"), tck=0.6*tck)
      if(incr==20)
        axis(side=1, labels=FALSE, at=outer(10,seq(-20,500,20), FUN="+"), tck=tck)
      mtext(side=2,line=2.5,text="$\\| \\mathcal{F} \\left[ {\\delta P}^a_\\mu(x) \\right] (f) \\|$")
      axis(side=2, labels=TRUE, at=10^(seq(-20,20,2)), tck=tck,las=1)                                                                                      
      axis(side=2, labels=FALSE, at=outer(2:9,10^(-20:19)),tck=0.6*tck)
      mtext(side=1,line=1.0,text="$f$")

      plot(y=x^2,
           x=seq(1,length(x))*tau/int.steps,
           type='l',pch='.',lwd=3,tck=tck,las=1,xaxt='n',
           xlab="",ylab="$\\left( \\delta P^a_\\mu(x,\\tau) \\right)^2$")
      mtext(side=1,line=1.0,text="$\\tau$")
      axis(side=1, labels=TRUE, at=seq(0,20,5),tck=tck)
      axis(side=1, labels=FALSE, at=outer(c(1:4),seq(0,20,5), FUN="+"), tck=0.6*tck)
    }
    tikz.finalize(tikzfiles)
    
    if(numerical){ 
      tex.basename <- "df_n"
      if(!missing(single) && !missing(pattern)) {
        tex.basename <- sprintf("%s.%s",pattern,tex.basename)
      }
      tikzfiles <- tikz.init(basename=tex.basename,width=width,height=height)
      irange <- 1:ncol(df_n[[length(df_n)]]$df)
      if(!missing(indices)){
        irange <- indices
      }
      for( i in irange ) {
        plot(y=df_n[[length(df_n)]]$df[,i],
             x=seq(1,length(df_n[[length(df_n)]]$df[,i]))*tau/int.steps,
             t='l',lwd=3,
             xlab="$\\tau$",ylab="$\\delta \\mathcal{P}^a_\\mu(x,\\tau)$")
        x <- df_n[[length(df_n)]]$df[,i]
        ft <- Mod(fft(x-mean(x)))/sqrt(length(x))
        plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1/tau,100),log='x',lwd=3,
        xlab="f",ylab="$\\left( \\Re \\left[ \\tilde{\\delta \\mathcal{P}}^a_\\mu(x,f) \\right] \\right)^2$")
      }
      tikz.finalize(tikzfiles)
    } #if(numerical)
  } #else(if(all))


  # this first more informative plot summarizes the differences between the most precise
  # determination of the analytical derivative ( df_a[[length(df_a)]] ) and all the others
  # for all lattice points, directions and generators, unless indices is defined
 
  if(length(df_a)>1){ 
    tex.basename <- "diff_a"
    if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s",pattern,tex.basename)
    }
    tikzfiles <- tikz.init(basename=tex.basename,width=width,height=height)
    for( i in 1:length(df_a) ) {
      irange <- 1:ncol(df_a[[length(df_a)]]$df)
      trange <- (1:nrow(df_a[[length(df_a)]]$df))/nrow(df_a[[length(df_a)]]$df)
      if(!missing(indices))
        irange <- indices
      for( idx in irange ) {
        plot(y=df_a[[length(df_a)]]$df[,idx]-df_a[[i]]$df[,idx],
             x=trange,
             main=sprintf("$\\delta P^a_\\mu \\left(x, \\tau, \\min( \\alpha_f ) \\right) - \\delta P^a_\\mu\\left(x,\\tau, \\alpha_f = 10^{%d} \\right)$",as.integer(log10(df_a[[i]]$fp))),
            xlab="$\\tau$", ylab="", pch=18, cex.main=0.9 )
      }
    }
    tikz.finalize(tikzfiles)
  }

  #stop("comparison reached")

  # and we do the same for the numerical derivative, checking what kind of acceptance
  # precision we need to select to see basically no difference to the most precise
  # one 
  if(numerical){
    if(length(df_n)>1){
      pdf.filename <- "diff_n.pdf"                                                                                                                                          
      if(!missing(single) && !missing(pattern)) {
        pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
      }
      pdf(file=pdf.filename,width=5,height=5,onefile=T) 
      for( i in 1:length(df_n) ) {
        irange <- 1:ncol(df_n[[length(df_n)]]$df)
        if(!missing(indices))
          irange <- indices
        for( idx in irange ) {
          plot(df_n[[length(df_n)]]$df[,idx]-df_n[[i]]$df[,idx],
               ylab=expression(dP[max]-dP[i]),xlab="t",
               main=sprintf("ap=%.1e eps=%.1e", df_n[[i]]$ap, df_n[[i]]$eps) )
        }
      }
      dev.off()
    }
  }

  if(numerical){
    a.idx <- length(df_a)
    diffs <- data.frame(prec=c(),eps=c(),diff=c())
    tex.basename <- "diff_a_n"
    if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s",pattern,tex.basename)
    }
    tikzfiles <- tikz.init(basename=tex.basename,width=width,height=height)
    for( i in 1:length(df_n) ) {
      irange <- 1:ncol(df_n[[length(df_n)]]$df)
      if(!missing(indices))
        irange <- indices
      for( idx in irange ) {
        plot(y=df_a[[a.idx]]$df[,idx]-df_n[[i]]$df[,idx],x=seq(1,length(df_n[[length(df_n)]]$df[,idx]))*tau/int.steps,
             main="$\\delta P^a_\\mu(x,\\tau)-\\delta \\mathcal{P}^a_\\mu(x,\\tau)$",
             ylab="",xlab="$\\tau$",type='l',lwd=3)
      }
      qt <- quantile( as.vector(df_a[[a.idx]]$df-df_n[[i]]$df), probs=c(0.1573,0.5,0.8427) )
      mx <- max(abs(df_a[[a.idx]]$df-df_n[[i]]$df))
      diffs <- rbind(diffs, data.frame(prec=df_n[[i]]$ap, eps=df_n[[i]]$eps, qt1=qt[1], qt2=qt[2], qt3=qt[3], max=mx) )
    }
    tikz.finalize(tikzfiles)
  }

  mean_df_a <- data.frame(df=c(),ddf=c(),sddf=c(),prec=c())
  mediandiff.df_a <- NULL
  diff.df_a <- matrix(nrow=0,ncol=length(df_a[[1]]$df))
  tex.basename <- "df_a_hist"
  if(!missing(single) && !missing(pattern)) {
    tex.basename <- sprintf("%s.%s",pattern,tex.basename)
  }
  tikzfiles <- tikz.init(basename=tex.basename,width=7,height=2.4)
  for( i in 1:length(df_a) ) {
    tdf <- mean(df_a[[i]]$df)
    tsdf <- sd(df_a[[i]]$df)/sqrt(length(df_a[[i]]$df))
    tddf <- df_a[[length(df_a)]]$df-df_a[[i]]$df
    mean_df_a <- rbind(mean_df_a,data.frame(df=tdf,sdf=tsdf,ddf=max(tddf),sddf=sd(tddf), prec=df_a[[i]]$ap))
    mediandiff.df_a <- rbind(mediandiff.df_a,quantile(as.vector(tddf),probs=c(0.1573,0.5,0.8427)))
    
    diff.df_a <- rbind(diff.df_a,as.vector(abs(df_a[[length(df_a)]]$df-df_a[[i]]$df)))
    cat("Number of breaks: ", floor(20*abs(max(df_a[[i]]$df)-min(df_a[[i]]$df))/sd(df_a[[i]]$df)), "\n")
    xlims <- c(mean(df_a[[i]]$df)-2*sd(df_a[[i]]$df),mean(df_a[[i]]$df)+2*sd(df_a[[i]]$df))
    hist(df_a[[i]]$df,breaks=function(x){floor(20*abs(max(x)-min(x))/sd(x))},main="",ylab="",xlim=xlims,freq=FALSE,xlab="$\\delta P^a_\\mu$")
    if(text){
      lims <- par('usr')
      qt <- quantile(df_a[[i]]$df,probs=c(0.01,0.5,0.99)) 
      text(x=lims[1]+0.18*(lims[2]-lims[1]),y=lims[4]-0.1*(lims[4]-lims[3]),labels=sprintf("$%5g~~%5g~~%5g$",qt[1],qt[2],qt[3]))
    }
    xlims <- c(0,3*sd(df_a[[i]]$df^2))
    hist(df_a[[i]]$df^2,breaks=function(x){floor(20*max(x)/sd(x))},main="",ylab="",xlim=xlims,freq=FALSE,xlab="$\\left( \\delta P^a_\\mu \\right)^2$")
  }
  tikz.finalize(tikzfiles)

  if(numerical){
    tex.basename <- "diff_a_n_hist"
    if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s",pattern,tex.basename)
    }
    tikzfiles <- tikz.init(basename=tex.basename,width=7,height=2.5)
    for( i in 1:length(df_n) ) {
      df <- df_a[[length(df_a)]]$df-df_n[[i]]$df
      print(floor(abs(max(df)-min(df))/0.5))
      hist(df,breaks=50,main="",ylab="",xlab="$\\delta P^a_\\mu(x,\\tau) - \\delta \\mathcal{P}^a_\\mu(x,\\tau)$",freq=FALSE)
      if(text){
        lims <- par('usr')
        qt <- quantile(df,probs=c(0.01,0.5,0.99)) 
        text(x=lims[1]+0.3*(lims[2]-lims[1]),y=lims[4]-0.1*(lims[4]-lims[3]),labels=sprintf("$%f %f %f$",qt[1],qt[2],qt[3]))
      }
      #hist(df_a[[i]]$df,breaks=2000,main="",ylab="",xlim=c(-200,200),xlab="$\\delta P^a_\\mu$",freq=FALSE)
    }
    tikz.finalize(tikzfiles)
  }
 
  if(length(df_a)>1){ 
    tex.basename <- "diff.df_a.prec"
    if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s",pattern,tex.basename)
    }
    tikzfiles <- tikz.init(basename=tex.basename,width=width,height=height)
    # set up plot area
    plot(x=log10(mean_df_a$prec),y=mediandiff.df_a[,2],ylab="",
                 main="$\\delta P^a_\\mu(x,\\tau,\\sigma_f=10^{-26}) - \\delta P^a_\\mu(x,\\tau,\\sigma_f)$",
                 xlab="$\\log_{10}(\\sigma_f)$",
                 type='n',ylim=c(-1e-5,1e-5),
                 xlim=c( min(log10(mean_df_a$prec)) , -3 ) )
    poly.col <- rgb(0.0,0.0,1.0,0.7)
    poly.x <- c(log10(mean_df_a$prec),rev(log10(mean_df_a$prec)))
    poly.y <- c(mediandiff.df_a[,1],rev(mediandiff.df_a[,3]))
    polygon(x=poly.x,y=poly.y,col=poly.col)
    tikz.finalize(tikzfiles)
  }
  
  # optim countour is useless, but the code is nice, let's preserve it for the future
  #diffs <- diffs[with(diffs,order(prec,eps)),]
  #diff.matrix <- matrix(data=diffs$diff,nrow=length(unique(diffs$eps)),ncol=length(unique(diffs$prec)))
  #pdf.filename <- "optim_contour.pdf"
  #if(!missing(single) && !missing(pattern)) {
  #  pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
  #}
  #pdf(file=pdf.filename,width=5,height=5,onefile=T) 
  #filled.contour(x=log10(sort(unique(diffs$eps))),y=log10(sort(unique(diffs$prec))),z=diff.matrix,ylab="eps",xlab="ap")
  #dev.off()
  if(numerical){ 
    tex.basename <- "optim"
    if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s",pattern,tex.basename)
    }
    # these need to be a little larger
    tikzfiles <- tikz.init(basename=tex.basename,width=(width+0.5),height=(height+1.2))
    plot(x=log10(diffs$eps),y=diffs$max,log="y",ylab="",xlab="$\\log_{10}(\\epsilon)$",type='n',main="")
    precs <- c(1e-5,1e-10,1e-15,1e-20,1e-25)
    logprecs <- as.integer(log10(precs))
    library(RColorBrewer)
    clr <- brewer.pal(n=length(precs),name="Set1")
    lt <- 1:length(precs)
    for(i in 1:length(precs)){
      df <- diffs[diffs$prec==precs[i],]
      lines(x=log10(df$eps),y=df$max,col=clr[i],lty=lt[i],lwd=3)
    }
    labels <- sprintf("$\\sigma_s=10^{%d}$",logprecs)
    legend(x="topright",legend=labels,lty=lt,col=clr,lwd=3,bty='n') 
    tikz.finalize(tikzfiles)
  }

  # do fourier analysis of complete data set
  ft <- rep(0,times=length(df_a[[1]]$df[,1]))
  for( i in 1:ncol(df_a[[length(df_a)]]$df) ) {
    # we do a fast fourier transform of the trajectory which provides us with a few of the dominant frequenc(y/ies)
    x <- df_a[[length(df_a)]]$df[,i]
    ft <- ft+Mod(fft(x-mean(x)))/sqrt(length(x))
  }
  ft <- ft[2:length(ft)]
  ft <- ft / max(ft)

  tex.basename <- "ft"
  if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s.%s",pattern,int.steps,tex.basename)
    }
  # these need to be a little larger
  tikzfiles <- tikz.init(basename=tex.basename,width=width,height=height)
  par(mgp=c(3,0.5,0))
  tck <- 0.2/height
  yrange <- range(ft)
  plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1,length(ft)/2)/(tau),#log='x',
       xlab="",ylab="$\\|\\mathcal{F} \\left[ \\delta P \\right] (f)\\|_\\mathrm{av} \\cdot N^{-1}$",lwd=3,las=1,yaxt='n',log='y',xaxt='n')
  abline(h=10^(-4),lty=3)
  mtext(side=1,line=1.3,text="$f$")
  incr <- 20
  if(length(ft)<4000)
    incr <- incr/2
  if(length(ft)<2000)
    incr <- incr/2

  axis(side=1, labels=TRUE, at=seq(0,500,incr),tck=-tck)
  if(incr>=10)
    axis(side=1, labels=FALSE, at=outer(c(2,4,6,8),seq(0,500,10), FUN="+"), tck=-0.6*tck)
  else
    axis(side=1, labels=FALSE, at=outer(c(1:4,6:9),seq(0,500,10), FUN="+"), tck=-0.6*tck)
  if(incr==20)
    axis(side=1, labels=FALSE, at=outer(10,seq(0,500,20), FUN="+"), tck=-0.8*tck)
  axis(side=2, labels=TRUE, at=10^(-10:10), tck=tck,las=1)                                                                                      
  axis(side=2, labels=FALSE, at=outer(2:9,10^(-10:9)),tck=0.5*tck)
  tikz.finalize(tikzfiles)
  
  tex.basename <- "ft.finite"
  if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s.%s",pattern,int.steps,tex.basename)
  }
  tikzfiles <- tikz.init(basename=tex.basename,width=3.4,height=2.5)
  par(mgp=c(3,0.3,0))
  yrange <- range(ft)
  plot(y=ft,t='l',x=seq(1,length(ft))/(tau),log='y',
       xlab="",ylab="$\\|\\mathcal{F} \\left[ \\delta P \\right] (f)\\|_\\mathrm{av} \\cdot N^{-1}$",lwd=3,las=1,yaxt='n',xaxt='n',
       xlim=c(1/tau,finite.fmax),
       ylim=10^c(-3,0))
  abline(h=10^(-2),lty=3)
  # indicate where ft crosses 10^-2, note how x and y are reversed for this purpose
  print(approx(y=seq(1,length(ft))/(tau),x=ft,xout=0.01))
  abline(v=approx(y=seq(1,floor(length(ft)/2))/(tau),x=ft[1:floor(length(ft)/2)],xout=0.01)$y,lty=3)
  mtext(side=1,line=1.0,text="$f$")
  axis(side=1, labels=TRUE, at=seq(-5,500,5), tck=0.1)
  axis(side=1, labels=FALSE, at=outer(c(1,2,3,4,6,7,8,9),seq(-20,500,10), FUN="+"), tck=0.05)
  #axis(side=1, labels=FALSE, at=outer(10,seq(-20,500,20), FUN="+"), tck=-0.045)
  axis(side=2, labels=TRUE, at=10^(-10:10), tck=0.05,las=1)                                                                                      
  axis(side=2, labels=FALSE, at=outer(2:9,10^(-10:9)),tck=0.025)
  tikz.finalize(tikzfiles)
  
  tex.basename <- "ft.loglog"
  if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s.%s",pattern,int.steps,tex.basename)
  }
  tikzfiles <- tikz.init(basename=tex.basename,width=3.4,height=2.5)
  par(mgp=c(3,0.3,0))
  yrange <- range(ft)
  plot(y=ft^2,t='l',x=seq(1,length(ft))/(tau),log='xy',
       xlab="",ylab="$\\|\\mathcal{F} \\left[ \\delta P \\right] (f)\\|_\\mathrm{av}^2 \\cdot N^{-1}$",lwd=3,las=1,yaxt='n',xaxt='n',
       xlim=c(1/tau,finite.fmax),
       ylim=10^c(-6,0))
  #abline(h=10^(-2),lty=3)
  # indicate where ft crosses 10^-2, note how x and y are reversed for this purpose
  #print(approx(y=seq(1,length(ft))/(tau),x=ft,xout=0.01))
  #abline(v=approx(y=seq(1,floor(length(ft)/2))/(tau),x=ft[1:floor(length(ft)/2)],xout=0.01)$y,lty=3)
  mtext(side=1,line=1.0,text="$f$")
  axis(side=1, labels=TRUE, at=seq(-5,500,5), tck=0.1)
  axis(side=1, labels=FALSE, at=outer(c(1,2,3,4,6,7,8,9),seq(-20,500,10), FUN="+"), tck=0.05)
  #axis(side=1, labels=FALSE, at=outer(10,seq(-20,500,20), FUN="+"), tck=-0.045)
  axis(side=2, labels=TRUE, at=10^(-10:10), tck=0.05,las=1)                                                                                      
  axis(side=2, labels=FALSE, at=outer(2:9,10^(-10:9)),tck=0.025)
  tikz.finalize(tikzfiles)

}
