# * create a pretty vertical plot of a set of quantities with an alternating
#   gray and white background, TeX labels and multiple error bars
#   if required
# * input is a vector of tex labels (labels), a vector of central values (x)
#   and a vector or data frame of errors (dx) and optionally a vector
#   or data frame of negative errors (mdx)
# * if the labels are missing, empty strings will be used instead
# * the basename (base of the resulting filename), x, dx as well as labelpos
#    are required arguments
# * additional arguments are passed to plotwitherror
# * it may be interesting to provide x, dx and mdx in relative terms where 
#   all(x==1) and everyhting else is normalized by the central value

plot_list_vertical <- function(basename,x,dx,mdx,labels,height,width,labelpos,ylim,mar,clr,pch,perlines,errband,cex.label=0.7,item.height=0.15,rep=FALSE,finalize=TRUE,...) {
  if(missing(height)) height <- length(x)*item.height+1 # the one is for the bottom margin
  if(missing(width)) width <- 5
  if(missing(mdx)) mdx <- dx
  if(missing(labels)) labels <- rep("",times=length(x))
  if(missing(ylim)) ylim <- c(0.6,length(x)+0.4)
  if(missing(clr)) clr <- "black"
  if(missing(pch)) pch <- 16
  
  tikzfiles <- NULL
  if(!rep){
    tikzfiles <- tikz.init(basename=basename,height=height,width=width,lwdUnit=0.5)
  }

  if(!missing(mar)){
    par(mar=mar)
  }
 
  y <- 1:length(x)
  # set up plot region
  if(!rep){
    par(mgp=c(2,0.75,0))
    plotwitherror(y=y,x=x,dx=dx,mdx=mdx,yaxt='n',ylab="",type='n',ylim=ylim,...)
    #mtext(text=xlabel,side=1,line=1.75)
    lims <- par("usr")
    # draw some gray rectangles to make reading the plot easier
    for( i in y ) {
      if( i %% 2 == 0 ) {
        col <- rgb(0.8,0.8,0.8,0.4)
        rect(lims[1],i-0.5,lims[2],i+0.5,col=col,border=NA)
      }
    }
    if(!missing(perlines)){
      abline(v=c(perlines-0.01,perlines+0.01),lty=2,col=rgb(0.5,0.5,0.5))
      abline(v=c(perlines-0.001,perlines+0.001),lty=3,col=rgb(0.5,0.5,0.5))
    }
  }

  # the errband argument allows nice rectangles to be placed to indicate an additional error
  # this comes in very handy when trying to compare errors and values in a ratio
  # with the error coming from the denominator depicted by the errorband
  if(!missing(errband)){
    rclr <- rgb(red=1.0,green=0.0,blue=0.0,alpha=0.3)
    print(errband)
    rect(xleft=errband$xleft,xright=errband$xright,ybottom=y-item.height,ytop=y+item.height,col=rclr,border=NA)
  }

  # do the actual plotting
  plotwitherror(y=y,x=x,dx=dx,mdx=mdx,rep=TRUE,pch=pch,col=clr,...)
  if(missing(labelpos)){
    lims <- par('usr')
    labelpos <- lims[1]
  }
  
  # if mar has been provided we probably want to put text outside the plot area
  # so we should disable clipping 
  if(!missing(mar)) par(xpd=NA)
  if(!is.list(labels)) {
    text(x=labelpos,y=y,labels=labels,adj=c(0.5,0.5), cex=cex.label)
  } else {
    increment <- 0.5/length(labels)
    for( i in 1:length(labels) ) {
      text(x=labelpos,y=y+increment*2*(i-1)-0.1,adj=c(0.5,0.5),labels=labels[[i]],cex=cex.label)
    }
  }
  if(finalize){
    tikz.finalize(tikzfiles)
  } else {
    tikzfiles
  }
}

