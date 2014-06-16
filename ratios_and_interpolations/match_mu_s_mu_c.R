match_mu_s_mu_c.2D <- function(name,alldat,mass_comb,pheno,mu_s,mu_c,m11,m12,m21,m22,xval="m11",debug=FALSE) {
    
  phenoband.x <- seq(0.02,0.03,length.out=50)
  nmass <- length(mass_comb$sc[,1])
  indices <- vector(mode='numeric', length=nmass)
  xvals <- vector(mode='numeric', length=nmass)
  
  # list of lists of length "number of bootstrap samples"
  dat <- vector( mode="list", length=length(which( alldat$val.tsboot$name == name ))/nmass )
  weights <- vector(mode='numeric', length=nmass)
  
  # pre-allocate memory for list elements
  for( d in 1:length(dat) ) {
    dat[[d]] <- data.frame(z=vector(mode='numeric',length=nmass),
      x=vector(mode='numeric',length=nmass),
      y=vector(mode='numeric',length=nmass))
  }
  
  # loop over mass combinations to construct list of data frames
  if( missing(m11) || missing(m12) || missing(m21) || missing(m22) ) { 
    for( i in 1:length(mass_comb$sc[,1]) ){
      if(debug) print(mass_comb$sc[i,])
      ts.indices <- which( alldat$val.tsboot$name == name &
                           alldat$val.tsboot$m11 == mass_comb$sc[i,]$m1 &
                           alldat$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
      for( j in 1:length(ts.indices) ) {
        dat[[j]][i,] <- c( z=alldat$val.tsboot$val[ts.indices[j]], x=mass_comb$sc[i,]$m1, y=mass_comb$sc[i,]$m2 )
      }
      indices[i] <- which( alldat$val$name == name &
                         alldat$val$m11 == mass_comb$sc[i,]$m1 &
                         alldat$val$m12 == mass_comb$sc[i,]$m2 )
      if( xval=="m11" ) {
        xvals[i] <- alldat$val$m11[indices[i]]
      } else if ( xval=="m12" ) {
        xvals[i] <- alldat$val$m12[indices[i]]
      } else if ( xval=="m12" ) {
        xvals[i] <- alldat$val$m21[indices[i]]
      } else {
        xvals[i] <- alldat$val$m22[indices[i]]
      }
      weights[i] <- (1/alldat$val$dval[indices[i]])^2
    }  
  } else {
    # loop over all combinations of m11,m12,m21 and m22
    i <- 1
    for( m11_i in m11 ) {
      for( m12_i in m12 ) {
        for( m21_i in m21 ) {
          for( m22_i in m22 ) {
            ts.indices <- which( alldat$val.tsboot$name == name &
                                 alldat$val.tsboot$m11 == m11_i &
                                 alldat$val.tsboot$m12 == m12_i &
                                 alldat$val.tsboot$m21 == m21_i & 
                                 alldat$val.tsboot$m22 == m22_i )
            #FIXME!! x and y choice need to be an option!
            for( j in 1:length(ts.indices) ) {
              dat[[j]][i,] <- c( z=alldat$val.tsboot$val[ts.indices[j]], x=m22_i, y=m12_i )
            }
            indices[i] <- which( alldat$val$name == name &
                               alldat$val$m11 == m11_i &
                               alldat$val$m12 == m12_i &
                               alldat$val$m21 == m21_i & 
                               alldat$val$m22 == m22_i )
            if( xval=="m11" ) {
              xvals[i] <- alldat$val$m11[indices[i]]
            } else if ( xval=="m12" ) {
              xvals[i] <- alldat$val$m12[indices[i]]
            } else if ( xval=="m12" ) {
              xvals[i] <- alldat$val$m21[indices[i]]
            } else {
              xvals[i] <- alldat$val$m22[indices[i]]
            }                          
            weights[i] <- (1/alldat$val$dval[indices[i]])^2 
            i <- i + 1
          }
        }
      }
    }
  }

  fit <- fit_linear_2d( z=dat, weights=weights, debug=debug )
  
  require(tikzDevice)
  texfile <- sprintf("%s.tex",name) 
  pdffile <- sprintf("%s.pdf",name)
  tikz(texfile, standAlone = TRUE, width=5, height=5)
  
  # with error propagation
  pred.tsboot <- extrapolate_2d(fit=fit,predx=mu_s$val,predy=mu_c$val,dpredx=mu_s$dval,dpredy=mu_c$dval)
  pred <- data.frame( z=apply(X=pred.tsboot$z,MARGIN=2,mean), dz=sqrt( apply(X=pred.tsboot$z,MARGIN=2,sd)^2 +
                      apply(X=pred.tsboot$dz,MARGIN=2,mean)^2 ) )
  if(debug) print(pred)
  
  # set up plot
  plotwitherror( y=alldat$val$val[indices], x=xvals, dy=alldat$val$dval[indices], 
                 main="with error propagation", ylab=alldat$val$texlabel[indices[1]], xlab="$\\mu_s$", type='n')
  # plot phenomenological value
  plot.confband( y=rep(pheno$val,50), dy=rep(pheno$dval,50), x=phenoband.x, col=rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3),
                 line=F )
  # add points on top of pheno value
  plotwitherror( y=alldat$val$val[indices], x=xvals, dy=alldat$val$dval[indices], rep=T )
  
  # this needs to be generalized!
    plotwitherror( y=pred$z[1], x=mu_s$val[1], dx=mu_s$dval[1], dy=pred$dz[1], rep=T, col='red', pch=15)
    plotwitherror( y=pred$z[2], x=mu_s$val[2], dx=mu_s$dval[2], dy=pred$dz[2], rep=T, col='blue', pch=17)
  
  legend(x=0.0212,y=1.23,
         legend=c("Measurements","Input: $mu_c = 0.0009*27.46(44)*11.85(16)$", "Input: $mu_c$ from $m_K/f_K=3.163(17)$","FLAG value"),
         col=c('black','red','blue',rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3)), 
         pch=c(1,15,17,15))  
  dev.off()
  tools::texi2dvi(texfile,pdf=T)
  command <- sprintf("pdfcrop %s %s",pdffile,pdffile)
  system(command)
  
}
