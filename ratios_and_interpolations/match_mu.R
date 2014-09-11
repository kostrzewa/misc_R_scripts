match_mu.1D <- function(name,alldat,masses,pheno,mu,lg,lg.coords,xlab="$\\mu$",debug=FALSE,solve=F,...) {
  if(debug) cat(sprintf("Doing 1D matching for %s\n",name))
  
  if(!missing(lg)){
    print(c(lg$col,pheno$col))
  }
  
  edge <- 2*(max(masses)-min(masses))
  phenoband.x <- c(min(masses)-edge,max(masses)+edge)
  
  nmass <- length(masses)
  
  # list of lists of length "number of bootstrap samples"
  dat <- vector( mode="list", length=length(which( alldat$val.tsboot$name == name ))/nmass)
  weights <- vector(mode='numeric', length=nmass)
  
  # pre-allocate memory for list elements
  for( d in 1:length(dat) ) {
    dat[[d]] <- data.frame(y=vector(mode='numeric',length=nmass),x=vector(mode='numeric',length=nmass))
  }
  
  # loop over mass combinations to construct list of data frames
  for(i in 1:nmass) {
    if(debug) print(masses[i])
    ts.indices <- which( alldat$val.tsboot$name == name &
                         alldat$val.tsboot$m12 == masses[i] )
                         
    for( j in 1:length(ts.indices) ) {
      dat[[j]][i,] <- c( y=alldat$val.tsboot$val[ts.indices[j]], x=masses[i] )
    }
    mass_index <- which( alldat$val$name == name &
                         alldat$val$m12 == masses[i] )
    weights[i] <- (1/alldat$val$dval[mass_index])^2 
  }
  
  indices <- which( alldat$val$name == name )

  fit <- fit_linear_1d.new( dat=dat, weights=weights )
  
  require(tikzDevice)
  texfile <- sprintf("%s.tex",name) 
  pdffile <- sprintf("%s.pdf",name)
  tikz(texfile, standAlone = TRUE, width=4, height=4)
  
  # extrapolate to the input masses mu$val
  predval.tsboot <- extrapolate_1d(fit=fit,predx=mu$val,dpredx=mu$dval)
  predval <- data.frame( val=apply(X=predval.tsboot$y,MARGIN=2,mean), 
                      dval=sqrt( apply(X=predval.tsboot$y,MARGIN=2,sd)^2 + apply(X=predval.tsboot$dy,MARGIN=2,mean)^2 ) )
  
  # find the point where the fitted function matches the pheno$val inputs
  predmu <- NULL
  if(!missing(pheno)) {
    predmu <- data.frame( val=rep(NA,length(pheno$val)), dval=rep(NA,length(pheno$val)) )
    
    if(solve) {
      predmu.tsboot <- solve_linear_1d( fit, pheno$val, pheno$dval )
      predmu <- data.frame( val=apply(X=predmu.tsboot$x,MARGIN=2,mean),
                          dval=sqrt( apply(X=predmu.tsboot$x,MARGIN=2,sd)^2 + predmu.tsboot$dx^2 ) )
    }
  }
  
  # set up plot
  plotwitherror( y=alldat$val$val[indices], x=alldat$val$m12[indices], dy=alldat$val$dval[indices], 
                 ylab=alldat$val$texlabel[indices[1]], xlab=xlab, type='n',...)
  # plot phenomenological value
  if(!missing(pheno)) {
    # remove any alpha value from pheno$col
    color.rgb <- col2rgb(pheno$col)/255
    bordercolor <- rgb(red=color.rgb[1],green=color.rgb[2],blue=color.rgb[3])
    rect( xleft=phenoband.x[1], xright=phenoband.x[2], ybottom=pheno$val-pheno$dval,
          ytop=pheno$val+pheno$dval, col=pheno$col, border=bordercolor )
  }
  # add band indicating solution
  if(solve) {
    # add band indicating solution
    rect( xleft=predmu$val-predmu$dval, xright=predmu$val+predmu$dval, ybottom=0,
          ytop=pheno$val+pheno$dval, col=rgb(red=0.0,blue=1.0,green=0.0,alpha=0.2), border='blue' )
  }
  
  # add points on top of bands
  plotwitherror( y=alldat$val$val[indices], x=alldat$val$m12[indices], dy=alldat$val$dval[indices], rep=T )
  
  # add predictions
  cols <- NULL
  if(!missing(lg)) {
    # remove the "data" colour
    cols <- lg$col[-1]
  } else {
    cols <- rainbow(n=length(predval$val))
  }
  plotwitherror( y=predval$val, x=mu$val, dx=mu$dval, dy=predval$dval, rep=T, col=cols, pch=15:17)
  
  # add point matching the pheno value
  if(!missing(pheno)) {
    if(solve) {
      plotwitherror( y=pheno$val, x=predmu$val, dx=predmu$dval, dy=pheno$dval, rep=T, col='blue', pch=18 )
    }
  }
  
  if( !missing(lg) ) {
    legend(x=lg.coords$x,y=lg.coords$y,
          legend=c(lg$labels,pheno$type),
          col=c(lg$col,pheno$col), 
          pch=c(lg$pch,pheno$pch), bg='white' )
  }
  dev.off()
  tools::texi2dvi(texfile,pdf=T)
  command <- sprintf("pdfcrop %s %s",pdffile,pdffile)
  system(command)
  
  if(!missing(pheno)) {
    return( list( predval=data.frame(val=predval$val,dval=predval$dval,mu=mu$val,dmu=mu$dval), 
                  predmu=data.frame(val=pheno$val,dval=pheno$dval,mu=predmu$val,dmu=predmu$dval) ) )
  } else {
    return( list( predval=data.frame(val=predval$val,dval=predval$dval,mu=mu$val,dmu=mu$dval), 
                  predmu=NA ) )
  }
}
