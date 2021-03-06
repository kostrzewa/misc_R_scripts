construct.empty.hadron_obs <- function() {
  hadron_obs <- list()
  attr(hadron_obs,"class") <- c("hadron_obs",class(hadron_obs))
  hadron_obs
}

construct.hadron_obs <- function(name,texlabel,m.sea,m.val,mean,err,boot) {
  hadron_obs <- construct.empty.hadron_obs()
  hadron_obs[[1]] <- list( name=name, texlabel=texlabel, m.sea=m.sea, m.val=m.val,
                          mean=mean, err=err, boot=boot )
}
  
select.hadron_obs <- function(hadron_obs,by,filter) {
#  if(!any(class(hadron_obs) == "hadron_obs")) {
#    stop("select.hadron_obs: hadron_obs argument must be of class hadron_obs!\n")
#  }
  rval <- list()
  for( i in hadron_obs ) {
    if( i[[by]] == filter )
      rval[[(length(rval)+1)]] <- i
  }
  attr(rval,"class") <- class(hadron_obs)
  rval
}

as.data.frame.hadron_obs <- function(hadron_obs) {
  if(!any(class(hadron_obs) == "hadron_obs")) {
   stop("as.data.frame.hadron_obs: hadron_obs argument must be of class hadron_obs!\n")
  }
}

# hadron_obs should already have been filtered to contain
# only one consistent observable
# pred.idx is a list of two index vecors which specify the
# predictor variables to be used for fitting the model later on

# when "fr" is specified, this indicates that the values
# are to be used for the determination of the systematic errror
# due to the choice of fitrange on the entire analysis chain
# starting with the 
# this means that 

extract.for.fes_fit <- function(hadron_obs, pred.idx, fr=FALSE, weighted=FALSE) {
  if(!any(class(hadron_obs) == "hadron_obs")) {
   stop("extract.for.fes_fit: hadron_obs argument must be of class hadron_obs!\n")
  }
  if(fr==TRUE && is.null(hadron_obs[[1]]$fr) ){
    stop("extract.for.fes_fit: extraction of values from fitrange analysis requested, but no 'fr' object found in 'hadron_obs'")
  }
  n.boot <- length(hadron_obs[[1]]$boot)
  # number of predictor variables
  n.pred <- length(pred.idx$m.sea)+length(pred.idx$m.val)
  pred.names <- sprintf("x%d",1:n.pred)
  # number of response points
  n.resp <- length(hadron_obs)
  
  rval <- vector(mode="list", length=n.boot)

  # number of extra columns
  n.cols <- 2
  col.names <- c("y",pred.names,"weight")
#  if(fr) {
#    n.cols <- 3
#    col.names <- c(col.names,"fr.weight")
#  }


  # pre-allocate memory
  for( i in 1:n.boot ) {
    rval[[i]] <- as.data.frame( array(dim=c(n.resp,n.pred+n.cols),
                                      dimnames=list(c(NULL),col.names)))
  }
  
  # construct the list of response variables with their predictors and weights
  # for weighted ChiSqr
  # at some point in the future one could actually consider a full covariance matrix
  # here!
  for( b in 1:n.boot ) {
    for( r in 1:n.resp ) {
      if(fr){
        # select a fitrange randomly and optionally take into account the weights
        n.fr <- ncol(hadron_obs[[r]]$fr$t)
        fitrange <- NULL
        if(weighted) {
          fitrange <- sample(x=1:n.fr,size=1,prob=hadron_obs[[r]]$fr$w)
        } else {
          fitrange <- sample(x=1:n.fr,size=1)
        }

        # select a specific fit-range independently for each response variable
        # (we use the mean of this result)
        # this way we get n.boot combinations of different fitranges,
        rval[[b]][r,] <- c(mean(hadron_obs[[r]]$fr$t[,fitrange]),
                           hadron_obs[[r]]$m.val[pred.idx$m.val],
                           hadron_obs[[r]]$m.sea[pred.idx$m.sea],
                           1/sd(hadron_obs[[r]]$fr$t[,fitrange])^2)
                          # hadron_obs[[r]]$fr$w[fitrange])
      } else {
        rval[[b]][r,] <- c(hadron_obs[[r]]$boot[b],
                           hadron_obs[[r]]$m.val[pred.idx$m.val],
                           hadron_obs[[r]]$m.sea[pred.idx$m.sea],
                           1/hadron_obs[[r]]$err^2)
      }
    }
  }
  
  rval
}

# a bit like above but extracts values and errors only
# assumes that hadron_obs already contains a consistent observable
# x.name specifies the list member which corresponds to the x variable ('m.val' or 'm.sea', currently)
# x.idx specifies the element of that list member which should be used as
# the x variable, this could be '2', say

extract.for.plot <- function(hadron_obs,x.name,x.idx) {
  N <- length(hadron_obs)
  rval <- as.data.frame( array(dim=c(N,4), dimnames=list(c(NULL),c("y","dy","mdy","x"))))
  for( i in 1:N ) {
    # in these plots we show only the statistical error
    rval[i,] <-  c(y=hadron_obs[[i]]$mean,
                   dy=sqrt(hadron_obs[[i]]$err^2),
                   mdy=sqrt(hadron_obs[[i]]$err^2),
                   x=hadron_obs[[i]][[x.name]][x.idx])
    print(rval[i,])
  }
  rval
}

# function which will produce a plot of data extracted from a hadron_obs object
# a number of optional items can be added to the plot
plot.hadron_obs <- function(df,name,pheno,extrapolations,solutions,lg,labelx,labely,debug=TRUE,...) {
  require(tikzDevice)
  temp <- sprintf("%s.%s",name,c("tex","pdf","aux","log"))
  tikzfiles <- list(tex=temp[1],pdf=temp[2],aux=temp[3],log=temp[4])
  rm(temp)
  tikz(tikzfiles$tex, standAlone = TRUE, width=4, height=4)

  if(debug){
    cat("Supllied data points\n")
    print(df)
  }
  
  # set up plot
  plotwitherror( y=df$y, x=df$x, dy=df$dy, mdy=df$mdy, type='n', ... )
  
  # extract plot boundaries
  lims <- par("usr")
  
  # coordinates which are guaranteed to extend beyond plot boundaries
  deltax <- abs(lims[2]-lims[1])
  xleft <- lims[1]-deltax
  xright <- lims[2]+deltax
  
  deltay <- abs(lims[4]-lims[3])
  ybottom <- lims[3]-deltay
  ytop <- lims[4]+deltay
  
  if(!missing(pheno)){
    # remove any alpha value from pheno$col
    opaque.rgb <- col2rgb(pheno$col)/255
    opaquecolor <- rgb(red=opaque.rgb[1],green=opaque.rgb[2],blue=opaque.rgb[3])
    
    if(debug){
      cat("Phenomenological value supplied\n")
      print(pheno)
    }
    
    rect( xleft=xleft, xright=xright, ybottom=pheno$val-pheno$dval,
          ytop=pheno$val+pheno$dval, col=pheno$col, border=opaquecolor )    

    if(!missing(solutions)){
      if(debug){
        cat("Solutions supplied\n")
        print(solutions)
      }
      
      # add band indicating solution      
      rect( xleft=solutions$val-solutions$dval, xright=solutions$val+solutions$dval, ybottom=ybottom,
            ytop=pheno$val+pheno$dval, col=rgb(red=0.0,blue=1.0,green=0.0,alpha=0.2), border='blue' )
    }
  }
  
  # add data points on top of any bands that were drawn
  plotwitherror( y=df$y, x=df$x, dy=df$dy, mdy=df$mdy, rep=TRUE )
  
  if(!missing(extrapolations)){
    if(debug){
      cat("Extrapolated points supplied\n")
      print(extrapolations)
    }  
    colours <- 'red'
    symbols <- 16
    if(!missing(lg)) {
        colours <- lg$col[2:(2+length(extrapolations$val))]
        symbols <- lg$pch[2:(2+length(extrapolations$val))]
    }
    
    plotwitherror( y=extrapolations$val, dy=extrapolations$dval, 
                   x=extrapolations[,extrapolations$plot.x.idx[1]],  dx=extrapolations[,extrapolations$plot.dx.idx[1]], 
                   col=colours, pch=symbols, rep=T )
  }
  
  # since this is the most important point, we draw it last
  if(!missing(solutions)){
    # and add point on top
    colours <- 'blue'
    symbols <- 18
    if(!missing(lg)) {
      colours <- lg$col[length(lg$col)]
      symbols <- lg$pch[length(lg$pch)]
    }
    plotwitherror( y=pheno$val, x=solutions$val, dx=solutions$dval, dy=pheno$dval, rep=TRUE, col=colours, pch=symbols )
  }        
  
  if(!missing(labelx)) {
    txpd <- par()$xpd
    par(xpd=NA)
    text(x=(lims[1]-0.05*deltax),
         y=(lims[3]-0.18*deltay),labels=labelx) 
    par(xpd=txpd)
  }
  
  if(!missing(lg)){
    #legend(x=lims[1],y=lims[4],legend=lg$labels,pch=lg$pch,col=lg$col,bty="n")
    legend(x="topleft",legend=lg$labels,pch=lg$pch,col=lg$col,bty="n")
  }
  
  dev.off()
  tools::texi2dvi(tikzfiles$tex,pdf=T)
  # use pdfcrop tool to remove plot borders
  command <- sprintf("pdfcrop %s %s",tikzfiles$pdf,tikzfiles$pdf)
  system(command)
  # remove temporary files 
  command <- sprintf("rm %s %s %s", tikzfiles$tex, tikzfiles$log, tikzfiles$aux)
  system(command)
}

plot_syserr_dist.hadron_obs <- function(hadron_obs,width=4.5,height=4,breaks=35) {
  require("weights")
  if(is.null(hadron_obs[[1]]$fr)){
    stop("plot_syserr_dist.hadron_obs: the 'hadron_obs' argument must have a non-null $fr member!\n")
  }
  
  # assume that this is a filtered hadron_obs
  basename <- sprintf("fitrange.syserr.dist.%s.raw",hadron_obs[[1]]$name)
  
  tikzfiles <- tikz.init(basename=basename,width=width,height=height)
 
  for( obs in hadron_obs ){ 
    w <- obs$fr$w
    fr <- apply(X=obs$fr$t,MARGIN=2,mean)
    qt <- quantile(x=fr,probs=c(0.1573,0.5,0.8427))
    wt.qt <- weighted.quantile(x=fr,w=w,probs=c(0.1573,0.5,0.8427))
    texlabel <- obs$texlabel
    
    col <- rgb(0.8,0.8,0.8,0.5)

    hst <- hist(x=fr,main="",yaxt='n',ylab="",xlim=c(min(fr),max(fr)),breaks=breaks,xlab="",freq=FALSE)
    lims <- par("usr")
    width <- abs(lims[2]-lims[1])
    height <- abs(lims[4]-lims[3])
    text(x=lims[2]-0.20*width,y=lims[4]-0.20*height,labels=texlabel)
    rect(xleft=qt[1],ybottom=lims[3],xright=qt[3],ytop=lims[4],col=col,border=NA)
    abline(v=qt[2],lwd=4)
    
    w.hst <- wtd.hist(x=fr,w=w,main="",yaxt='n',ylab="",breaks=breaks,xlab="",xlim=c(min(fr),max(fr)),freq=FALSE)
    lims <- par("usr")
    width <- abs(lims[2]-lims[1])
    height <- abs(lims[4]-lims[3])
    text(x=lims[2]-0.20*width,y=lims[4]-0.20*height,labels=texlabel)
    rect(xleft=wt.qt[1],ybottom=lims[3],xright=wt.qt[3],ytop=lims[4],col=col,border=NA)
    abline(v=wt.qt[2],lwd=4)
  }

  tikz.finalize(tikzfiles)
  return(list(hst,w.hst))
}

# return a data frame to summarize the contents of the hadron_obs object, this is of course rather silly
# because the number of columns could be completely variable and depend on the observable... this
# was the actual reason for introducting hadron_obs in the first place...
summary.hadron_obs <- function(hadron_obs) {
  rval <- NULL
  for(obs in hadron_obs){
    rval <- rbind(rval,data.frame(name=obs$name,texlabel=as.character(obs$texlabel), m1=obs$m.sea[1],m2=obs$m.val[1],m3=obs$m.val[2],
                                  mean=obs$mean,err=obs$err,median=obs$median,mserr=obs$serr[1],serr=obs$serr[2],stringsAsFactors=FALSE))
  }
  class(rval) <- c(class(rval),"summary_hadron_obs")
  rval
}

#plot_fitrange_histograms.hadron_obs <- function(hadron_obs,width=5,height=4,...){
#  tikzfiles <- tikz.init(basename="fitrange_histograms",width=width,height=height)
#  for( obs in hadron_obs ){
#    obs.medians <- apply( X=obs$fr$t, MARGIN=2, FUN=median )

plot_fitrange_errors.summary_hadron_obs <- function(summary_hadron_obs,width=6,height=50,...) {
  if(!any(class(summary_hadron_obs)=="summary_hadron_obs")){
    stop("plot.summary_hadron_obs: argument must be of class 'summary_hadron_obs'")
  }

  tikzfiles <- tikz.init("fitrange_summary",width=width,height=height,sanitize=FALSE)

  # we want to plot everything normalized by the mean
  Ncol <- ncol(summary_hadron_obs)
  summary_hadron_obs[,3:Ncol] <- summary_hadron_obs[,3:Ncol] / summary_hadron_obs$mean
  
  y <- 1:nrow(summary_hadron_obs)
  plotwitherror(x=summary_hadron_obs$mean,y=y,
                dx=cbind(summary_hadron_obs$err,summary_hadron_obs$serr),
                mdx=cbind(summary_hadron_obs$err,summary_hadron_obs$mserr),...)
  #points(x=summary_hadron_obs$median,y=y,pch=4)
  text(x=0.955,y=y,labels=summary_hadron_obs$texlabel)

  tikz.finalize(tikzfiles)
}

compare.extrapolations.syserr.hadron_obs <- function(hadron_obs,ext.sys) {
  rval <- NULL
  for( name in unique(summary(hadron_obs)[,1])){
    obs <- select.hadron_obs(hadron_obs,'name',name)
    for( ext in ext.sys ) {
      if( all(name == ext$name) ) {
        rval <- rbind(rval,data.frame(name=name,mserr.ratio=ext$mserr/obs[[1]]$serr[1],serr.ratio=ext$serr/obs[[1]]$serr[2]))
      }
    }
  }
  rval
}



# it is unclear to me whether this kind of function is useful...
# because it would need to be provided with quite a number of arguments
# fes.hadron_obs <- function(name, hadron_obs, ) {
#   pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
#   obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
#   
#   pred.idx <- list(m.val=c(2,3),m.sea=vector())
#   dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
#   fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
#   pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
#   fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)
# 
#   extrapolations <- rbind(extrapolations, 
#                           data.frame(name=name,
#                                     val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
#                                     dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ),
#                                     x=c(mu_s$val,mu_s$val[3]), dx=c(mu_s$dval,mu_s$dval[3])
#                                     )
#                         )
# 
#   df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
#   plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[extrapolations$name==name,],#solutions=solution,
#                   xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)
# }
