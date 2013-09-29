plot_oriented_plaquettes <- function(filename, start, end=-1) {
  name <- strsplit(filename,split=c('/','.'))
  i_start <- start+1
  i_end <- end+1
  plaqdat <- read.table(filename)

  if(end == -1) {
    i_end <- length(plaqdat$V1)
  }
  
  orplaqs <- c( 
    mean(plaqdat$V2[i_start:i_end]), 
    mean(plaqdat$V3[i_start:i_end]),
    mean(plaqdat$V4[i_start:i_end]),
    mean(plaqdat$V5[i_start:i_end]),
    mean(plaqdat$V6[i_start:i_end]),
    mean(plaqdat$V7[i_start:i_end]) )

  dorplaqs <- sqrt( 
    c(
      var(plaqdat$V2[i_start:i_end]), 
      var(plaqdat$V3[i_start:i_end]),
      var(plaqdat$V4[i_start:i_end]),
      var(plaqdat$V5[i_start:i_end]),
      var(plaqdat$V6[i_start:i_end]),
      var(plaqdat$V7[i_start:i_end]) 
    )/(i_end-i_start)  )

  par(mfrow=c(1,2),oma=c(0,0,1.5,0))

  labels <- c("TX","TY","TZ","XY","XZ","YZ") 
  plotwitherror(x=seq(1,6,1),y=orplaqs,dy=dorplaqs,ylab="oriented <P>",xlab="hyperplane",xaxt="n")
  axis(1,labels=labels,tick=TRUE,at=seq(1,6,1))

  ymax <- max(plaqdat$V2[i_start:i_end],plaqdat$V3[i_start:i_end],plaqdat$V4[i_start:i_end],plaqdat$V5[i_start:i_end],plaqdat$V6[i_start:i_end],plaqdat$V7[i_start:i_end])
  ymin <- min(plaqdat$V2[i_start:i_end],plaqdat$V3[i_start:i_end],plaqdat$V4[i_start:i_end],plaqdat$V5[i_start:i_end],plaqdat$V6[i_start:i_end],plaqdat$V7[i_start:i_end])
  
  plot(t='l',ylab="oriented <P>",xlab=expression(t[HMC]),
    x=seq(start,i_end-1,1),
    y=plaqdat$V2[i_start:i_end],
    ylim=c(ymin,ymax),
    pch='.')

  lines(x=seq(start,i_end-1,1),y=plaqdat$V3[i_start:i_end],col="red",pch='.')
  lines(x=seq(start,i_end-1,1),y=plaqdat$V4[i_start:i_end],col="green",pch='.')
  lines(x=seq(start,i_end-1,1),y=plaqdat$V5[i_start:i_end],col="blue",pch='.')
  lines(x=seq(start,i_end-1,1),y=plaqdat$V6[i_start:i_end],col="cyan",pch='.')
  lines(x=seq(start,i_end-1,1),y=plaqdat$V7[i_start:i_end],col="magenta",pch='.')

  title <- name[[1]][length(name[[1]])-1]
  print(title)
  mtext(title,outer=TRUE,cex=1.3)
}
