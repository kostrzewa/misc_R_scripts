mDs_ov_fDs.2D.demo <- function(results,debug=FALSE,mu_s_from_mK) {
  
  name <- "m_Ds_ov_f_Ds"
  
  # extrapolate f_Ds/f_D in mu_s/mu_c
  mu_s <- c(0.0009*27.46,mu_s_from_mK$mu_s)
  dmu_s <- c(0.44*0.0009,mu_s_from_mK$dmu_s)
  
  mu_c <- mu_s*11.85
  dmu_c <- c( 0.0009*sqrt( (11.85*0.44)^2 + (27.46*0.16)^2 ), 
              sqrt( (mu_s_from_mK$dmu_s*11.85)^2 + (mu_s_from_mK$mu_s*0.16)^2 ) )
  
  # phenomenological values of m_Ds/f_Ds
  pheno <- c(7.9)
  dpheno <- c(0.2)  
  
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
  
  pdf("mDs_ov_fDs.pdf")
  
  # with error propagation
  pred.tsboot <- extrapolate_2d(fit=fit,predx=mu_s,predy=mu_c,dpredx=dmu_s,dpredy=dmu_c)
  pred <- data.frame( z=apply(X=pred.tsboot$z,MARGIN=2,mean), dz=sqrt( apply(X=pred.tsboot$z,MARGIN=2,sd)^2 +
                      apply(X=pred.tsboot$dz,MARGIN=2,mean)^2 ) )
  # set up plot
  plotwitherror( y=results$val$val[indices], x=results$val$m11[indices], dy=results$val$dval[indices], 
                 main="with error propagation", ylab=expression(m[D[s]]/f[D[s]]), xlab=expression(mu[s]), type='n')
  # plot phenomenological value
  plot.confband( y=rep(pheno[1],50), dy=rep(dpheno[1],50), x=phenoband.x, col=rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3),
                 line=F )
  # add points on top of pheno value
  plotwitherror( y=results$val$val[indices], x=results$val$m11[indices], dy=results$val$dval[indices], rep=T )
  plotwitherror( y=pred$z[1], x=mu_s[1], dx=dmu_s[1], dy=pred$dz[1], rep=T, col='red', pch=15)
  plotwitherror( y=pred$z[2], x=mu_s[2], dx=dmu_s[2], dy=pred$dz[2], rep=T, col='blue', pch=17)
  legend(x=0.0212,y=1.23,
         legend=c("Measurements","Input: mu_c = 0.0009*27.46(44)*11.85(16)", "Input: mu_c via mu_s from m_K/f_K=3.163(17)","FLAG value"),
         col=c('black','red','blue',rgb(red=0.0,blue=0.0,green=1.0,alpha=0.3)), 
         pch=c(1,15,17,15))  
         
  dev.off()
}
