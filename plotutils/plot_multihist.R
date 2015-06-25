require("hadron")
plot_multihist <- function(dat,lg,basename="multihist",factor=8,pts=FALSE,width=6,height=3,...) {
  tikzfiles <- tikz.init(basename,width=width,height=height)
  cols <- c("red","blue")
  if(length(dat)>2){  
    require("RColorBrewer")
    cols <- brewer.pal(n=length(dat),name="Set1")
  }
  pcols <- cols
  cols <- t(col2rgb(cols))
  cols <- cbind(cols,rep(100,times=nrow(cols)))
  cols <- apply(X=cols,MARGIN=1,FUN=function(x){ rgb(x[1],x[2],x[3],x[4],maxColorValue=255) })
  lims <- c(min(unlist(dat)),max(unlist(dat)))
  dev <- sd(unlist(dat))/factor
  breaks <- seq(lims[1]-dev,lims[2]+dev,dev)
  hst <- list()
  ylim <- c(0,0)
  py <- 0
  
  # determine a sensible ylim
  for( i in 1:length(dat) ) {
    hst[[i]] <- hist(dat[[i]],breaks=breaks,plot=FALSE)
    ymax <- max(hst[[i]]$counts)+0.1*max(hst[[i]]$counts)
    if(ymax>ylim[2]) {
      ylim[2] <- ymax
      py <- ymax - 0.05*max(hst[[i]]$counts)
    }
  }
   
  plot(hst[[1]],col=cols[1],xlim=c(lims[1]-dev,lims[2]+dev),ylim=ylim,...)
  for( i in 2:length(dat) ) {
    plot(hst[[i]],add=TRUE,col=cols[i])
  }
  if(pts){
    for( i in 1:length(dat) ){
      uw <- uwerr(data=dat[[i]])
      plotwitherror(x=uw$value,dx=uw$dvalue,y=py,rep=TRUE,pch='.',col=pcols[i])
    }
  }
  if(!missing(lg)){
    legend(x="topleft",legend=lg,pch=15,col=cols,bty='n')
  }
  tikz.finalize(tikzfiles)
}
