barplot_gpuperf <- function(data,width=4,height=4,title="",basename="gpuperf",ngpu=4) {
  require(RColorBrewer)
  cols <- brewer.pal(3,name="Set1")
  perfmtx <- matrix(nrow=2,ncol=length(data))
  for( i in 1:length(data) ){
    # data contains GPU usage in percent, sampled once a second
    # interpret > 50 as in use
    perfmtx[,i] <- c(sum(data[[i]]>=50)/ngpu,sum(data[[i]]<50)/ngpu)
  }
  totals <- apply(X=perfmtx,MARGIN=2,FUN=sum)
  offset <- max(totals)/20

  tikzfiles <- tikz.init(basename=basename,width=width,height=height)
  mids <- barplot(perfmtx,las=1,ylab="runtime [sec]",main=title,col=cols,names.arg=c("original","QUDA solver","QUDA solver \n + threads"))
  par(xpd=NA)
  for( i in 1:length(totals) ){
    text(x=mids[i],y=totals[i]+offset,labels=sprintf("%.0f sec",totals[i]))
  }
  legend(x="topright",legend=c("GPU idle","GPU active"),fill=rev(cols[1:2]))
  tikz.finalize( tikzfiles )
}

