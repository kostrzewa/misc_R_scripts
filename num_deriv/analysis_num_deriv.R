analysis_num_deriv <- function(dir,volume=2^4,int.steps=101,trajs=1,filename,endian="little",idx=1,single=FALSE) {
  
  n_pat <- "*_f_numerical.bin"
  a_pat <- "*_f_analytical.bin"
  if(single) {
    a_pat <- "f_numerical.bin"
    n_pat <- "f_analytical.bin"
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

  pdf(file="df_a.pdf",onefile=T,width=6,height=6)
  for( i in 1:length(df_a) ) {
    plot(df_a[[i]]$df[,idx])
  }
  dev.off()

  pdf(file="df_a.all.pdf",width=6,height=6)
  library(plot3D)
  points3D(x=rep(times=ncol(df_a[[length(df_a)]]$df),seq(1,nrow(df_a[[length(df_a)]]$df))),
          y=rep(times=nrow(df_a[[length(df_a)]]$df),seq(1,ncol(df_a[[length(df_a)]]$df))),
          z=as.vector(df_a[[length(df_a)]]$df),pch=".",
          xlab="int. step",ylab="idx")
  for( i in 1:ncol(df_a[[length(df_a)]]$df) ) {
    plot(df_a[[length(df_a)]]$df[,i])
  }
  dev.off()

  pdf(file="df_n.pdf",onefile=T,width=6,height=6)
  for( i in 1:length(df_n) ) {
    plot(df_n[[i]]$df[,idx])
  }
  dev.off()

  pdf(file="diff_a.pdf",onefile=T,width=6,height=6)
  for( i in 1:length(df_a) ) {
    plot(df_a[[length(df_a)]]$df[,idx]-df_a[[i]]$df[,idx])
  }
  dev.off()

  pdf(file="diff_n.pdf",onefile=T,width=6,height=6)
  for( i in 1:length(df_a) ) {
    plot(df_n[[length(df_n)]]$df[,idx]-df_n[[i]]$df[,idx])
  }
  dev.off()
  
  a.idx <- length(df_a)
  diffs <- data.frame(prec=c(),eps=c(),diff=c())
  pdf(file="diff_a_n.pdf",onefile=T,width=6,height=6)
  for( i in 1:length(df_n) ) {
    plot(df_a[[a.idx]]$df[,idx]-df_n[[i]]$df[,idx])
    diffs <- rbind(diffs, data.frame(prec=df_n[[i]]$ap,eps=df_n[[i]]$eps,diff=max(abs(df_a[[a.idx]]$df-df_n[[i]]$df))))
  }
  dev.off()

  mean_df_a <- data.frame(df=c(),ddf=c(),sddf=c(),prec=c())
  diff.df_a <- matrix(nrow=0,ncol=length(df_a[[1]]$df))
  pdf("df_a_hist.pdf",width=6,height=6)
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
  
  pdf(file="mean.df_a.pdf",width=6,height=6)
  plot(x=log10(mean_df_a$prec),y=mean_df_a$df,ylab="mean(df_a)",main="average analytical derivative",xlab="log10(force precision)")
  plotwitherror(x=log10(mean_df_a$prec),y=mean_df_a$ddf,dy=mean_df_a$sddf,log="y",ylab="max(abs(df_25-df_x))",main="max analytical derivative difference",xlab="log10(force precision)")
  for( i in 1:nrow(diff.df_a) ) {
    hist(diff.df_a[i,],breaks=20)
  }
  dev.off()  


  diffs <- diffs[with(diffs,order(prec,eps)),]
  diff.matrix <- matrix(data=diffs$diff,nrow=length(unique(diffs$eps)),ncol=length(unique(diffs$prec)))
  pdf(file="optim_contour.pdf",width=6,height=6)
  filled.contour(x=log10(sort(unique(diffs$eps))),y=log10(sort(unique(diffs$prec))),z=diff.matrix,ylab="eps",xlab="ap")
  dev.off()
  
  pdf(file="optim.pdf",onefile=T,width=6,height=6)
  plot(x=log10(diffs$eps),y=diffs$diff,log="y",ylab="max. diff",main="max. diff. btw. num. and analytical deriv.",xlab="log10(eps)")
  plot(x=log10(diffs$prec),y=diffs$diff,log="y",ylab="max. diff",main="max. diff. btw. num. and analytical deriv.",xlab="log10(acc. precision)")
  dev.off()
}
