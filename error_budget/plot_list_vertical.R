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

plot_list_vertical <- function(basename,x,dx,mdx,labels,height,width,labelpos,item.height=0.15,...) {
  if(missing(height)) height <- length(x)*item.height
  if(missing(width)) width <- 5
  if(missing(mdx)) mdx <- dx
  if(missing(labels)) labels <- rep("",times=length(x))
  tikzfiles <- tikz.init(basename=basename,height=height,width=width)

  y <- 1:length(x)
  # set up plot region
  plotwitherror(y=y,x=x,dx=dx,mdx=mdx,xlab="",yaxt='n',ylab="",type='n',...)
  lims <- par("usr")
  # draw some gray rectangles to make reading the plot easier
  for( i in y ) {
    if( i %% 2 == 0 ) {
      col <- rgb(0.8,0.8,0.8,0.4)
      rect(lims[1],i-0.5,lims[2],i+0.5,col=col,border=NA)
    }
  }
  # do the actual plotting
  plotwitherror(y=y,x=x,dx=dx,mdx=mdx,xlab="",yaxt='n',rep=T,...)
  text(x=labelpos,y=y,labels=labels)

  tikz.finalize(tikzfiles)
}

