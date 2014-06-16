mK_ov_fK.demo <- function(results,strange_masses,debug=FALSE) {
  # extrapolate m_K/f_K in mu_s 
  predx = c(0.0009*27.46)
  dpredx = c(0.0004)
  # find mu s.t. m_K/f_K physical
  predy = c(3.163)
  dpredy = c(0.017)
  
  n.tsboot <- length( which( results$val.tsboot$name == "m_K_ov_f_K" ) )/length(strange_masses)
  
  # choose indices for "m_K/f_K"
  indices <- which( results$val$name == "m_K_ov_f_K" )
  ts.indices <- array(dim=c(n.tsboot,length(strange_masses)))

  for( i in 1:length(strange_masses) ) {
    ts.indices[,i] <- which( results$val.tsboot$name == "m_K_ov_f_K" & results$val.tsboot$m12 == strange_masses[i] )
  }

  legendlabels <- c("Measurements","Input: $ a \\mu_s = 27.46(44) a \\mu_l $","Input: $ m_K/f_K = 3.163(17)$","PDG value")
  legendcols <- c('black','red','blue',rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3))
  legendsyms <- c(1,15,17,15)
  
  tsboot.dat <- cbind( results$val.tsboot$val[ts.indices[,1]], 
                results$val.tsboot$val[ts.indices[,2]], 
                results$val.tsboot$val[ts.indices[,3]],
                results$val.tsboot$val[ts.indices[,4]] )

  # compute standard confidence interval
  confband.x = seq(0.021,0.028,length.out=50)
  fit <- fit_linear_1d(y=tsboot.dat, x=strange_masses, type="lm", 
                              weights=1/results$val$dval[indices]^2, debug=debug)
  confband.y <- extrapolate_1d(fit=fit, predx=confband.x)
  
  # repeat fit using nls so that error propagation will work further below
  fit <- fit_linear_1d(y=tsboot.dat, x=strange_masses, type="nls", 
                                weights=1/results$val$dval[indices]^2, debug=debug)

  filebase <- "mK_ov_fK"
  require(tikzDevice)
  texfile <- sprintf("%s.tex",filebase) 
  pdffile <- sprintf("%s.pdf",filebase)
  tikz(texfile, standAlone = TRUE, width=5, height=5)
  
  # with error propagation from dpredx 
  newy <- extrapolate_1d(fit=fit, predx=predx, dpredx=dpredx)
  # set up plot area
  plotwitherror(y=results$val$val[indices], dy=results$val$dval[indices], x=results$val$m12[ indices ],
                xlim=c(0.022,0.0265),main="with error propagation",
                xlab="$a \\mu_s$",ylab=results$val$texlabel[indices[1]],type='n')
  plot.confband(y=apply(X=confband.y$y,MARGIN=2,FUN=mean),dy=apply(X=confband.y$y,MARGIN=2,FUN=sd),x=confband.x,col=rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4))
  # add phenomenological value
  plot.confband(y=rep(predy,50), dy=rep(dpredy,50), x=confband.x, col=rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3),
                 line=F )  
  # add points
  plotwitherror(y=results$val$val[indices], dy=results$val$dval[indices], x=results$val$m12[ indices ],rep=T)
  # our predicted values
  # the first error is from propagating the error in predx, the second error is evaluated on the bootstrap sample
  plotwitherror(y=apply(X=newy$y,MARGIN=2,mean), dy=sqrt(apply(X=newy$dy,MARGIN=2,mean)^2+apply(X=newy$y,MARGIN=2,sd)^2),x=predx,rep=T,col="red",pch=15)
  # error propagation from dpredy
  newx <- solve_linear_1d( fit=fit, predy=predy, dpredy=dpredy )
  plotwitherror(x=apply(X=newx$x,MARGIN=2,FUN=mean),dx=sqrt(apply(X=newx$x,FUN=sd,MARGIN=2)^2+newx$dx^2),y=predy,rep=T,col='blue',pch=17)
  legend(x=0.0219,y=3.209,pch=legendsyms,col=legendcols,legend=legendlabels)
  dev.off()
  tools::texi2dvi(texfile,pdf=T)
  command <- sprintf("pdfcrop %s %s",pdffile,pdffile)
  system(command)
  
  return( data.frame( mu_s=apply(X=newx$x,MARGIN=2,FUN=mean), dmu_s=sqrt(apply(X=newx$x,FUN=sd,MARGIN=2)^2+newx$dx^2) ) )
}
