match_mu.1D <- function(name,alldat,masses,pheno,mu,lg.coords,xlab="$\\mu$",debug=FALSE,...) {
    
  edge <- 0.1*(max(masses)-min(masses))
  phenoband.x <- seq(min(masses)-edge,max(masses)+edge,length.out=50)
  
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
  tikz(texfile, standAlone = TRUE, width=5, height=5)
  
  # with error propagation
  pred.tsboot <- extrapolate_1d(fit=fit,predx=mu$val,dpredx=mu$dval)
  pred <- data.frame( y=apply(X=pred.tsboot$y,MARGIN=2,mean), 
                      dy=sqrt( apply(X=pred.tsboot$y,MARGIN=2,sd)^2 + apply(X=pred.tsboot$dy,MARGIN=2,mean)^2 ) )
  # set up plot
  plotwitherror( y=alldat$val$val[indices], x=alldat$val$m12[indices], dy=alldat$val$dval[indices], 
                 ylab=alldat$val$texlabel[indices[1]], xlab=xlab, type='n',...)
  # plot phenomenological value
  plot.confband( y=rep(pheno$val,50), dy=rep(pheno$dval,50), x=phenoband.x, col=rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3),
                 line=F )
  # add points on top of pheno value
  plotwitherror( y=alldat$val$val[indices], x=alldat$val$m12[indices], dy=alldat$val$dval[indices], rep=T )
  plotwitherror( y=pred$y[1], x=mu$val[1], dx=mu$dval[1], dy=pred$dy[1], rep=T, col='red', pch=15)
  plotwitherror( y=pred$y[2], x=mu$val[2], dx=mu$dval[2], dy=pred$dy[2], rep=T, col='blue', pch=17)
  legend(x=lg.coords$x,y=lg.coords$y,
         legend=c("Measurements","Input: $mu_c = 0.0009*27.46(44)*11.85(16)$", "Input: $mu_c$ from $m_K/f_K=3.163(17)$","FLAG value"),
         col=c('black','red','blue',rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3)), 
         pch=c(1,15,17,15))  
  dev.off()
  tools::texi2dvi(texfile,pdf=T)
  command <- sprintf("pdfcrop %s %s",pdffile,pdffile)
  system(command)
  
}
