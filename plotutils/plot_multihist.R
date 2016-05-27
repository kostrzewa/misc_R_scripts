require("hadron")
plot_multihist <- function(dat,lg,cols,lg.pos,xlim.probs=c(0.01,0.99),label.x="",label.y="",
                           basename="multihist",factor=8,pts=FALSE,width=6,height=3,
                           medians=FALSE,init=TRUE,...) {
  tikzfiles <- NULL
  if(init) tikzfiles <- tikz.init(basename,width=width,height=height)
  
  par(mgp=c(3,0.6,0))
  if(missing(cols)){
    cols <- c("red","blue")
    if(length(dat)>2){  
      require("RColorBrewer")
      cols <- brewer.pal(n=length(dat),name="PuOr")
    }
  }
  pcols <- cols
  cols <- t(col2rgb(cols))
  cols <- cbind(cols,rep(150,times=nrow(cols)))
  cols <- apply(X=cols,MARGIN=1,FUN=function(x){ rgb(x[1],x[2],x[3],x[4],maxColorValue=255) })
  lims <- quantile(unlist(dat),probs=xlim.probs)
  print(lims)
  dev <- sd(unlist(dat))/factor
  med <- median(unlist(dat))
  breaks <- seq(min(unlist(dat))-dev,max(unlist(dat))+dev,dev)
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
   
  plot(hst[[1]],col=cols[1],xlim=c(lims[1],lims[2]),ylim=ylim,xlab="",yaxt='n',...)
  if(length(hst)>1){
    for( i in 2:length(dat) ) {
      plot(hst[[i]],add=TRUE,col=cols[i])
    }
  }
  if(medians){
    for( i in 1:length(dat) )
      abline(v=median(dat[[i]]),col=cols[i],lwd=3)
  }

  axis(side=2,tck=0.05,labels=TRUE,las=1)
  mtext(line=1.4,text=label.x,side=1)
  if(pts){
    for( i in 1:length(dat) ){
      uw <- uwerr(data=dat[[i]])
      plotwitherror(x=uw$value,dx=uw$dvalue,y=py,rep=TRUE,pch='.',col=pcols[i])
    }
  }
  if(!missing(lg)){
    if(missing(lg.pos)) lg.pos <- "topright"
    legend(x=lg.pos,legend=lg,pch=NA,fill=cols,col=cols,bty='n',pt.cex=1.3,)
  }
  if(init) tikz.finalize(tikzfiles)
}
