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

analysis_num_deriv <- function(dir,filename,pattern,indices,volume=2^4,int.steps=101,trajs=1,tau=1,endian="little",single=FALSE,all=FALSE,width=4,height=4) {
 
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
      ft <- Re(fft(df_a[[length(df_a)]]$df[,i]))^2/length(df_a[[length(df_a)]]$df[,i])^2
      plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1/tau,100),log='x',
           xlab="f",ylab="$Re(dP(f))^2$",lwd=3 )
    }
    tikz.finalize(tikzfiles)
  
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
      ft <- Re(fft(df_n[[length(df_n)]]$df[,i]))^2/length(df_a[[length(df_a)]]$df[,i])^2
      plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1/tau,100),log='x',
      xlab="f",ylab="$Re(dP(f))^2$",lwd=3 )
    }
    tikz.finalize(tikzfiles)
  } else {
    tex.basename <- "df_a"
    if(!missing(single) && !missing(pattern)) {
      tex.basename <- sprintf("%s.%s",pattern,tex.basename)
    }
    tikzfiles <- tikz.init(basename=tex.basename,width=width,height=height)
    irange <- 1:ncol(df_a[[length(df_a)]]$df)
    if(!missing(indices)){
      irange <- indices
    }
    for( i in irange ) {
      plot(y=df_a[[length(df_a)]]$df[,i],
           x=seq(1,length(df_a[[length(df_a)]]$df[,i]))*tau/int.steps,
           t='l',lwd=3,
           xlab="$\\tau$",ylab="$\\delta P^a_\\mu(x,\\tau)$")
      # we do a fast fourier transform of the trajectory which provides us with a few of the dominant frequenc(y/ies) 
      ft <- Re(fft(df_a[[length(df_a)]]$df[,i]))^2/length(df_a[[length(df_a)]]$df[,i])^2
      plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1/tau,100),log='x',lwd=3,
           xlab="f",ylab="$\\left( \\Re \\left[ \\tilde{\\delta P}^a_\\mu(x,f) \\right] \\right)^2$")
    }
    tikz.finalize(tikzfiles)
    
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
      ft <- Re(fft(df_n[[length(df_n)]]$df[,i]))^2/length(df_a[[length(df_a)]]$df[,i])^2
      plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1/tau,100),log='x',lwd=3,
      xlab="f",ylab="$\\left( \\Re \\left[ \\tilde{\\delta \\mathcal{P}}^a_\\mu(x,f) \\right] \\right)^2$")
    }
    tikz.finalize(tikzfiles)
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
    #hist(df_a[[i]]$df,breaks=as.integer(2*max(df_a[[i]]$df)),main="",ylab="",xlab="",xlim=c(-20,20),freq=FALSE)#,xlab="$\\delta P^a_\\mu$")
    hist(df_a[[i]]$df,breaks=50,main="",ylab="",xlab="",freq=TRUE)#,xlab="$\\delta P^a_\\mu$")
    #hist(df_a[[i]]$df,breaks=2000,main="",ylab="",xlim=c(-200,200),xlab="$\\delta P^a_\\mu$",freq=FALSE)
  }
  tikz.finalize(tikzfiles)


  tex.basename <- "diff_a_n_hist"
  if(!missing(single) && !missing(pattern)) {
    tex.basename <- sprintf("%s.%s",pattern,tex.basename)
  }
  tikzfiles <- tikz.init(basename=tex.basename,width=7,height=2.5)
  for( i in 1:length(df_n) ) {
    hist(df_a[[length(df_a)]]$df-df_n[[i]]$df,breaks=50,main="",ylab="",freq=TRUE,xlab="$\\delta P^a_\\mu(x,\\tau) - \\delta \\mathcal{P}^a_\\mu(x,\\tau)$")
    #hist(df_a[[i]]$df,breaks=2000,main="",ylab="",xlim=c(-200,200),xlab="$\\delta P^a_\\mu$",freq=FALSE)
  }
  tikz.finalize(tikzfiles)

  
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
  
  tex.basename <- "optim"
  if(!missing(single) && !missing(pattern)) {
    tex.basename <- sprintf("%s.%s",pattern,tex.basename)
  }
  # these need to be a little larger
  tikzfiles <- tikz.init(basename=tex.basename,width=(width+0.5),height=(height+1.2))
  plot(x=log10(diffs$eps),y=diffs$max,log="y",ylab="",xlab="$\\log(\\epsilon)$",type='n',main="")
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

  #i <- 1
  #clr <- rainbow(n=length(unique(diffs$prec)))
  #for(prec in unique(diffs$prec)){
  #  df <- diffs[diffs$prec==prec,]
  #  lines(x=log10(df$eps),y=df$diff,col=clr[i])
  #  i <- i+1
  #}
  #i <- 1
  #clr <- rainbow(n=length(unique(diffs$eps)))
  #plot(x=log10(diffs$prec),y=diffs$diff,log="y",ylab="max. diff",main="max. diff. btw. num. and analytical deriv.",xlab="log10(acc. precision)",type='n')
  #for(eps in unique(diffs$eps)){
  #  df <- diffs[diffs$eps==eps,]
  #  lines(x=log10(df$prec),y=df$diff,col=clr[i])
  #  i<-i+1
  #}

  #dev.off()
}
