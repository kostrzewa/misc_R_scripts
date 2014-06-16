fDs_ov_fD.2D.demo <- function(results,debug=FALSE,mu_s_from_mK) {
  
  name <- "f_Ds_ov_f_D"
  
  # extrapolate f_Ds/f_D in mu_s/mu_c
  mu_s <- c(0.0009*27.46,mu_s_from_mK$mu_s)
  dmu_s <- c(0.44*0.0009,mu_s_from_mK$dmu_s)
  
  mu_c <- mu_s*11.85
  dmu_c <- c( 0.0009*sqrt( (11.85*0.44)^2 + (27.46*0.16)^2 ), 
              sqrt( (mu_s_from_mK$dmu_s*11.85)^2 + (mu_s_from_mK$mu_s*0.16)^2 ) )
  
  # phenomenological values of f_Ds/f_D
  pheno <- c(1.187,1.26)
  dpheno <- c(0.012,0.06) 
  
  phenoband.x <- seq(0.02,0.03,length.out=50)
  
  zeros <- rep(0,2)
  
  # list of lists of length "number of bootstrap samples"
  dat <- vector( mode="list", length=length(which( results$val.tsboot$name == name ))/length(mass_comb$sc[,1]))
  weights <- vector(mode='numeric', length=16)
  
  # pre-allocate memory for list elements
  for( d in 1:length(dat) ) {
    dat[[d]] <- data.frame(z=vector(mode='numeric',length=16),x=vector(mode='numeric',length=16),y=vector(mode='numeric',length=16))
  }
  
  # loop over f_Ds mass combinations to construct list of data frames
  for( i in 1:length(mass_comb$sc[,1]) ){
    if(debug) print(mass_comb$sc[i,])
    ts.indices <- which( results$val.tsboot$name == name &
                         results$val.tsboot$m11 == mass_comb$sc[i,]$m1 &
                         results$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
    for( j in 1:length(ts.indices) ) {
      dat[[j]][i,] <- c( z=results$val.tsboot$val[ts.indices[j]], x=mass_comb$sc[i,]$m1, y=mass_comb$sc[i,]$m2 )
    }
    sc_index <- which( results$val$name == name &
                       results$val$m11 == mass_comb$sc[i,]$m1 &
                       results$val$m12 == mass_comb$sc[i,]$m2 )
    weights[i] <- (1/results$val$dval[sc_index])^2 
  }
  
  indices <- which( results$val$name == name )

  fit <- fit_linear_2d( z=dat, weights=weights )
  
  pdf("fDs_ov_fD.pdf")
  
  # with error propagation
  pred.tsboot <- extrapolate_2d(fit=fit,predx=mu_s,predy=mu_c,dpredx=dmu_s,dpredy=dmu_c)
  pred <- data.frame( z=apply(X=pred.tsboot$z,MARGIN=2,mean), dz=sqrt( apply(X=pred.tsboot$z,MARGIN=2,sd)^2 +
                      apply(X=pred.tsboot$dz,MARGIN=2,mean)^2 ) )
  # set up plot
  plotwitherror( y=results$val$val[indices], x=results$val$m11[indices], dy=results$val$dval[indices], 
                 main="with error propagation", ylab=expression(f[D[s]]/f[D]), xlab=expression(mu[s]), type='n')
  # plot phenomenological value
  plot.confband( y=rep(pheno[1],50), dy=rep(dpheno[1],50), x=phenoband.x, col=rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3),
                 line=F )
  # add points on top of pheno value
  plotwitherror( y=results$val$val[indices], x=results$val$m11[indices], dy=results$val$dval[indices], rep=T )
  plotwitherror( y=pred$z[1], x=mu_s[1], dx=dmu_s[1], dy=pred$dz[1], rep=T, col='red', pch=15)
  plotwitherror( y=pred$z[2], x=mu_s[2], dx=dmu_s[2], dy=pred$dz[2], rep=T, col='blue', pch=17)
  legend(x=0.0212,y=1.23,
         legend=c("Measurements","Input: mu_c = 0.0009*27.46(44)*11.85(16)", "Input: mu_c from m_K/f_K=3.163(17)","FLAG value"),
         col=c('black','red','blue',rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3)), 
         pch=c(1,15,17,15))  
  
  dev.off()
}

fDs_ov_fD.demo <- function(results,debug=FALSE,mu_s_from_mK) {
  # extrapolate f_Ds/f_D in mu_c for mu_s=0.0238 
  predx <- c(0.0009*27.46*11.85,mu_s_from_mK$mu_s*11.85)
  dpredx <- c( 0.0009*sqrt( (11.85*0.16)^2 + (27.46*0.44)^2), 
              sqrt( (mu_s_from_mK$dmu_s*11.85)^2 + (mu_s_from_mK$mu_s*0.16)^2 ) )
  # find mu_s s.t. f_Ds/f_pi = x
  predy <- c(1.187,1.26)
  dpredy <- c(0.012,0.06)
  
  name <- "f_Ds_ov_f_D"
  
  ## TODO: prepare data into correct format with a little bit more logic!
  
  n.tsboot <- length( which( results$val.tsboot$name == name & results$val.tsboot$m11 == strange_masses[1] ) )/length(charm_masses)
  
  # choose indices for "f_Ds/f_D"
  indices <- which( results$val$name == name & results$val$m11 == strange_masses[1] )
  ts.indices <- array(dim=c(n.tsboot,length(charm_masses)))

  for( i in 1:length(charm_masses) ) {
    ts.indices[,i] <- which( results$val.tsboot$name == name & 
                             results$val.tsboot$m11 == strange_masses[1] & 
                             results$val.tsboot$m12 == charm_masses[i] )
  }

  legendlabels <- c("Measurements",expression(Input: mu[c]),expression(Input: f[D[s]]/f[D]))
  legendcols <- c('black','red','blue')
  legendsyms <- c(1,15,17)
  
  tsboot.dat <- cbind( results$val.tsboot$val[ts.indices[,1]], 
                results$val.tsboot$val[ts.indices[,2]], 
                results$val.tsboot$val[ts.indices[,3]],
                results$val.tsboot$val[ts.indices[,4]] )

  # compute standard confidence interval
  confband.x = seq(0.20,0.4,length.out=50)
  fit <- fit_linear_1d(y=tsboot.dat, x=charm_masses, type="lm", 
                              weights=1/results$val$dval[indices]^2, debug=debug)
  confband.y <- extrapolate_1d(fit=fit, predx=confband.x)
  
  # repeat fit using nls so that error propagation will work further below
  fit <- fit_linear_1d(y=tsboot.dat, x=charm_masses, type="nls", 
                                weights=1/results$val$dval[indices]^2, debug=debug)
 
 
  # without error propagation from dpredx
  newy <- extrapolate_1d(fit=fit, predx=predx)
  pdf(file=sprintf("%s_with_without_error.pdf",name))
  par(family="Times")
  plotwitherror(y=results$val$val[indices], dy=results$val$dval[indices], x=results$val$m12[ indices ],
                xlim=c(0.21,0.39),main="without error propagation",xlab=expression(mu[c]),ylab=expression(f[D[s]]/f[D]))
  plot.confband(y=apply(X=confband.y$y,MARGIN=2,FUN=mean),dy=apply(X=confband.y$y,MARGIN=2,FUN=sd),x=confband.x,col=rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4))
  # our predicted values
  plotwitherror(y=apply(X=newy$y,MARGIN=2,mean), dy=apply(X=newy$y,MARGIN=2,sd),x=predx,rep=T,col="red",pch=15)
  # without error propagation from dpredy
  newx <- solve_linear_1d( fit=fit, predy=predy )
  plotwitherror(x=apply(X=newx$x,MARGIN=2,FUN=mean),dx=apply(X=newx$x,FUN=sd,MARGIN=2),y=predy,rep=T,col='blue',pch=17)
  legend(x=0.21,y=1.22,pch=legendsyms,col=legendcols,legend=legendlabels)
  
  # with error propagation from dpredx 
  newy <- extrapolate_1d(fit=fit, predx=predx, dpredx=dpredx)
  plotwitherror(y=results$val$val[indices], dy=results$val$dval[indices], x=results$val$m12[ indices ],
                xlim=c(0.21,0.39),main="with error propagation",xlab=expression(mu[c]),ylab=expression(f[D[s]]/f[D]))
  plot.confband(y=apply(X=confband.y$y,MARGIN=2,FUN=mean),dy=apply(X=confband.y$y,MARGIN=2,FUN=sd),x=confband.x,col=rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4))
  # our predicted values
  # the first error is from propagating the error in predx, the second error is evaluated on the bootstrap sample
  plotwitherror(y=apply(X=newy$y,MARGIN=2,mean), dy=sqrt(apply(X=newy$dy,MARGIN=2,mean)^2+apply(X=newy$y,MARGIN=2,sd)^2),x=predx,dx=dpredx,rep=T,col="red",pch=15)
  # error propagation from dpredy
  newx <- solve_linear_1d( fit=fit, predy=predy, dpredy=dpredy )
  plotwitherror(x=apply(X=newx$x,MARGIN=2,FUN=mean),dx=sqrt(apply(X=newx$x,FUN=sd,MARGIN=2)^2+newx$dx^2),y=predy,dy=dpredy,rep=T,col='blue',pch=17)
  legend(x=0.21,y=1.22,pch=legendsyms,col=legendcols,legend=legendlabels)
  dev.off()
}
