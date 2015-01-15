analysis_num_deriv <- function(dir,filename,pattern,indices,volume=2^4,int.steps=101,trajs=1,tau=1,endian="little",single=FALSE) {
  
  n_pat <- "*_f_numerical.bin"
  a_pat <- "*_f_analytical.bin"
  
  if(single) {
    a_pat <- "f_analytical.bin"
    n_pat <- "f_numerical.bin"
    if(!missing(pattern)) {
      a_pat <- sprintf("%s_%s",pattern,a_pat)
      n_pat <- sprintf("%s_%s",pattern,n_pat)
    }
    if(!missing(filename)) {
      a_pat <- n_pat <- filename
    }
  }
  
  n_files <- list.files(dir,pattern=n_pat,full.names=FALSE)
  a_files <- list.files(dir,pattern=a_pat,full.names=FALSE)
  
  df_a <- list()
  df_n <- list()

  for( filename in n_files ) {
    n_file <- file(filename,"rb")
    ap <- readBin(n_file,double(),n=1,endian=endian)
    fp <- readBin(n_file,double(),n=1,endian=endian)
    eps <- readBin(n_file,double(),n=1,endian=endian)
    
    df_n[[length(df_n)+1]] <- list(ap=ap,fp=fp,eps=eps,df=t(array(data=readBin(n_file,double(),n=trajs*int.steps*volume*4*8,endian=endian),dim=c(volume*4*8,int.steps*trajs))))

    close(n_file)
  }

  precs <- c()

  for( filename in a_files ) {
    a_file <- file(filename,"rb")
    ap <- readBin(a_file,double(),n=1,endian=endian)
    fp <- readBin(a_file,double(),n=1,endian=endian)
    eps <- readBin(a_file,double(),n=1,endian=endian)
    
    if(!any( abs(precs-fp) < 10^-30)) {
      df_a[[length(df_a)+1]] <- list(ap=ap,fp=fp,eps=eps,df=t(array(data=readBin(a_file,double(),n=trajs*int.steps*volume*4*8,endian=endian),dim=c(volume*4*8,int.steps*trajs))))
      precs <- c(precs,fp)
    }

    close(a_file)
  }

  pdf.filename <- "df_a.all.pdf"
  if(!missing(single) && !missing(pattern)) {
    pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
  }

  pdf(file=pdf.filename,width=5,height=5,onefile=T)
  library(plot3D)
  rows <- nrow(df_a[[length(df_a)]]$df)
  cols <- ncol(df_a[[length(df_a)]]$df)
  xcoords <- rep(times=cols,seq(1,rows))
  ycoords <- c()
  for( i in 1:cols ) {
    idcs <- seq((i-1)*rows,i*rows)
    ycoords[idcs] <- i
  }
  points3D(x=xcoords,
          y=ycoords,
          z=as.vector(df_a[[length(df_a)]]$df[1:rows,1:cols]),pch=".",
          xlab="t",ylab="idx",zlab="dP")

  for( i in 1:ncol(df_a[[length(df_a)]]$df) ) {
    plot(y=df_a[[length(df_a)]]$df[,i],
         x=seq(1,length(df_a[[length(df_a)]]$df[,i]))*tau/int.steps,
         t='l',lwd=3,
         xlab="t",ylab="dP")

    ft <- Re(fft(df_a[[length(df_a)]]$df[,i]))^2/length(df_a[[length(df_a)]]$df[,i])^2
    plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1/tau,100),log='x',
    xlab="f",ylab="Re(dP(f))^2",lwd=3 )
  }
  dev.off()


  pdf.filename <- "df_n.all.pdf"
  if(!missing(single) && !missing(pattern)) {
    pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
  }

  pdf(file=pdf.filename,width=5,height=5,onefile=T)

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
  for( i in 1:ncol(df_n[[length(df_n)]]$df) ) {
    plot(y=df_n[[length(df_n)]]$df[,i],
         x=seq(1,length(df_n[[length(df_n)]]$df[,i]))*tau/int.steps,
         t='l',xlab="t",ylab="dP",lwd=3)
    ft <- Re(fft(df_n[[length(df_n)]]$df[,i]))^2/length(df_a[[length(df_a)]]$df[,i])^2
    plot(y=ft,t='l',x=seq(1,length(ft))/(tau),xlim=c(1/tau,100),log='x',
    xlab="f",ylab="Re(dP(f))^2",lwd=3 )
  }
  dev.off()

  pdf.filename <- "diff_a.pdf"
  if(!missing(single) && !missing(pattern)) {
    pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
  }
  pdf(file=pdf.filename,width=5,height=5,onefile=T) 
  for( i in 1:length(df_a) ) {
    irange <- 1:ncol(df_a[[length(df_a)]]$df)
    if(!missing(indices))
      irange <- indices
    for( idx in irange ) {
      plot(df_a[[length(df_a)]]$df[,idx]-df_a[[i]]$df[,idx],
           ylab=expression(dP[max]-dP[i]),xlab="t",
           main=sprintf("fp=%.1e",df_a[[i]]$fp) )
    }
  }
  dev.off()

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
  
  a.idx <- length(df_a)
  diffs <- data.frame(prec=c(),eps=c(),diff=c())
  pdf.filename <- "diff_a_n.pdf"
  if(!missing(single) && !missing(pattern)) {
    pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
  }
  pdf(file=pdf.filename,width=5,height=5,onefile=T) 
  for( i in 1:length(df_n) ) {
    irange <- 1:ncol(df_n[[length(df_n)]]$df)
    if(!missing(indices))
      irange <- indices
    for( idx in irange ) {
      plot(y=df_a[[a.idx]]$df[,idx]-df_n[[i]]$df[,idx],x=seq(1,length(df_n[[length(df_n)]]$df[,idx]))*tau/int.steps,
           main=sprintf("eps=%.1e ap=%.1e fp=%.1e",df_n[[i]]$eps,df_n[[i]]$ap,df_a[[a.idx]]$fp),
           ylab=expression(dP[a]-dP[n]),xlab="t")
    }
    diffs <- rbind(diffs, data.frame(prec=df_n[[i]]$ap,eps=df_n[[i]]$eps,diff=max(abs(df_a[[a.idx]]$df-df_n[[i]]$df))))
  }
  dev.off()

  mean_df_a <- data.frame(df=c(),ddf=c(),sddf=c(),prec=c())
  diff.df_a <- matrix(nrow=0,ncol=length(df_a[[1]]$df))
  pdf.filename <- "df_a_hist.pdf"
  if(!missing(single) && !missing(pattern)) {
    pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
  }
  pdf(file=pdf.filename,width=5,height=5,onefile=T) 
  for( i in 1:length(df_a) ) {
    tdf <- mean(df_a[[i]]$df)
    tsdf <- sd(df_a[[i]]$df)
    tddf <- abs(df_a[[length(df_a)]]$df-df_a[[i]]$df)
    mean_df_a <- rbind(mean_df_a,data.frame(df=tdf,sdf=tsdf,ddf=max(tddf),sddf=sd(tddf), prec=df_a[[i]]$ap))
    diff.df_a <- rbind(diff.df_a,as.vector(abs(df_a[[length(df_a)]]$df-df_a[[i]]$df)))
    hist(df_a[[i]]$df,breaks=36,main=df_a[[i]]$ap)
  }
  dev.off()
 
  print(mean_df_a)

  pdf.filename <- "mean_df_a.pdf"
  if(!missing(single) && !missing(pattern)) {
    pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
  }
  pdf(file=pdf.filename,width=5,height=5,onefile=T) 
  plot(x=log10(mean_df_a$prec),y=mean_df_a$df,ylab="mean(df_a)",main="average analytical derivative",xlab="log10(force precision)")
  plotwitherror(x=log10(mean_df_a$prec),y=mean_df_a$ddf,dy=mean_df_a$sddf,log="y",ylab="max(abs(df_25-df_x))",main="max analytical derivative difference",xlab="log10(force precision)")
  for( i in 1:nrow(diff.df_a) ) {
    hist(diff.df_a[i,],breaks=20)
  }
  dev.off()  

  diffs <- diffs[with(diffs,order(prec,eps)),]
  diff.matrix <- matrix(data=diffs$diff,nrow=length(unique(diffs$eps)),ncol=length(unique(diffs$prec)))
  pdf.filename <- "optim_contour.pdf"
  if(!missing(single) && !missing(pattern)) {
    pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
  }
  pdf(file=pdf.filename,width=5,height=5,onefile=T) 
  filled.contour(x=log10(sort(unique(diffs$eps))),y=log10(sort(unique(diffs$prec))),z=diff.matrix,ylab="eps",xlab="ap")
  dev.off()
  
  pdf.filename <- "optim.pdf"
  if(!missing(single) && !missing(pattern)) {
    pdf.filename <- sprintf("%s.%s",pattern,pdf.filename)
  }
  pdf(file=pdf.filename,width=5,height=5,onefile=T) 
  plot(x=log10(diffs$eps),y=diffs$diff,log="y",ylab="max. diff",main="max. diff. btw. num. and analytical deriv.",xlab="log10(eps)")
  plot(x=log10(diffs$prec),y=diffs$diff,log="y",ylab="max. diff",main="max. diff. btw. num. and analytical deriv.",xlab="log10(acc. precision)")
  dev.off()
}
