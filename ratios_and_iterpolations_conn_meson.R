library('propagate')

source("/home/bartek/code/R/misc_R_scripts/fit_extrapolate_solve.R")

ratios_and_iterpolations_conn_meson <- function(debug=F,recompute=T) {
  # masses to be used in this analysis
  strange_masses <<- c(0.0238,0.0245,0.0252,0.0259)
  charm_masses <<- c(0.2822,0.294,0.3058,0.3176)
  light_masses <<- c(0.0009)
  
  # combinations of these masses
  mass_comb <<- list( ll=data.frame( m1=light_masses, m2=light_masses ),
                     ls=expand.grid( m1=light_masses, m2=strange_masses),
                     lc=expand.grid( m1=light_masses, m2=charm_masses),
                     sc=expand.grid( m1=strange_masses, m2=charm_masses) )
  
  # File (and object) names of the data to be loaded
  names <- list(ll_c=sprintf("ll_c.m%g.m%g.matrixfit",mass_comb$ll$m1,mass_comb$ll$m2),
                #ll_na=sprintf("lln_u_%g-d_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                #ll_nb=sprintf("lln_d_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ls_c=sprintf("ls_c.m%g.m%g.matrixfit",mass_comb$ls$m1,mass_comb$ls$m2),
                lc_c=sprintf("lc_c.m%g.m%g.matrixfit",mass_comb$lc$m1,mass_comb$lc$m2),
                sc_c=sprintf("sc_c.m%g.m%g.matrixfit",mass_comb$sc$m1,mass_comb$sc$m2) )

  # single quantities
  quants <- NULL
  quants[["m_pi"]] <- list(name="m_pi", m1=light_masses, m2=light_masses, datanames=names$ll_c, type="mps" )
  quants[["m_K"]] <- list(name="m_K", m1=light_masses, m2=strange_masses, datanames=names$ls_c, type="mps" )
  quants[["m_D"]] <- list(name="m_D", m1=light_masses, m2=charm_masses, datanames=names$lc_c, type="mps" )
  quants[["m_Ds"]] <- list(name="m_Ds", m1=strange_masses, m2=charm_masses, datanames=names$sc_c, type="mps" )
  quants[["f_pi"]] <- list(name="f_pi", m1=light_masses, m2=light_masses, datanames=names$ll_c, type="fps" )
  quants[["f_K"]] <- list(name="f_K", m1=light_masses, m2=strange_masses, datanames=names$ls_c, type="fps" )
  quants[["f_D"]] <- list(name="f_D", m1=light_masses, m2=charm_masses, datanames=names$lc_c, type="fps" )
  quants[["f_Ds"]] <- list(name="f_Ds", m1=strange_masses, m2=charm_masses, datanames=names$sc_c, type="fps" )

  # the different ratios that we would like to compute
  ratios <- NULL
  ratios[["m_pi_ov_f_pi"]] <- list( name="m_pi_ov_f_pi", dividend=quants[["m_pi"]], divisor=quants[["f_pi"]])
  ratios[["m_K_ov_f_K"]] <- list( name="m_K_ov_f_K", dividend=quants[["m_K"]], divisor=quants[["f_K"]])
  ratios[["m_D_ov_f_D"]] <- list( name="m_D_ov_f_D", dividend=quants[["m_D"]], divisor=quants[["f_D"]])
  ratios[["m_Ds_ov_f_Ds"]] <- list( name="m_Ds_ov_f_Ds", dividend=quants[["m_Ds"]], divisor=quants[["f_Ds"]])

  ratios[["m_K_ov_f_pi"]] <- list( name="m_K_ov_f_pi", dividend=quants[["m_K"]], divisor=quants[["f_pi"]])
  ratios[["m_D_ov_f_pi"]] <- list( name="m_D_ov_f_pi", dividend=quants[["m_D"]], divisor=quants[["f_pi"]])
  ratios[["m_Ds_ov_f_pi"]] <- list( name="m_Ds_ov_f_pi", dividend=quants[["m_Ds"]], divisor=quants[["f_pi"]])

  ratios[["f_K_ov_f_pi"]] <- list( name="f_K_ov_f_pi", dividend=quants[["f_K"]], divisor=quants[["f_pi"]])
  ratios[["f_D_ov_f_pi"]] <- list( name="f_D_ov_f_pi", dividend=quants[["f_D"]], divisor=quants[["f_pi"]])
  ratios[["f_Ds_ov_f_pi"]] <- list( name="f_Ds_ov_f_pi", dividend=quants[["f_Ds"]], divisor=quants[["f_pi"]])

  ratios[["f_D_ov_f_K"]] <- list( name="f_D_ov_f_K", dividend=quants[["f_D"]], divisor=quants[["f_K"]])
  ratios[["f_Ds_ov_f_K"]] <- list( name="f_Ds_ov_f_K", dividend=quants[["f_Ds"]], divisor=quants[["f_K"]])

  ratios[["f_Ds_ov_f_D"]] <- list( name="f_Ds_ov_f_D", dividend=quants[["f_Ds"]], divisor=quants[["f_D"]])

  if(recompute) {
    # read all the data files
    for( n in names ) {
      for( file in n ) {
        file <- sprintf("%s.Rdata",file)
        if(debug) cat("Loading" ,file,"\n")
        load(file)
      }
    }
  
    results <- list( res = data.frame(), res.tsboot = data.frame() )
    
    for( qty in quants ) {
      for( objname in qty$datanames ) {
        dat <- get.obs.tsboot(obj=get(objname), type=qty$type)
        results <- compute.quant( name=qty$name, dat=dat, res=results )
      }
    }
    
    for( r in ratios ) {
      # there is a special case if the data of dividend and divisor are exactly the same,
      # we don't need x*x cases but only x!
      if(debug) cat("Computing ratio", r$name, "\n")
      if( all(r$dividend$datanames == r$divisor$datanames) ) {
        for( objname in r$dividend$datanames ) {
          dividend <- get.obs.tsboot(obj=get(objname), type=r$dividend$type )
          divisor <- get.obs.tsboot(obj=get(objname), type=r$divisor$type )
          results <- compute.ratio(name=r$name,dividend=dividend,divisor=divisor,res=results)
        }
      } else {

        # when any pair of mass vectors from the two observables are the same
        # we meed to ensure that "off-diagonal" elements, where say the strange
        # mass for the divisor and dividend are different, are skipped!
        m1m1 <- FALSE
        m1m2 <- FALSE
        m2m1 <- FALSE
        m2m2 <- FALSE
        if( all(r$dividend$m1 == r$divisor$m1) ) {
          m1m1 <- TRUE
        } else if( all(r$dividend$m1 == r$divisor$m2) ) {
          m1m2 <- TRUE
        } else if( all(r$dividend$m2 == r$divisor$m1) ) {
          m2m1 <- TRUE
        } else if( all(r$dividend$m2 == r$divisor$m2) ) {
          m2m2 <- TRUE
        }

        for( dividend.objname in r$dividend$datanames ) {
          for( divisor.objname in r$divisor$datanames ) {
            dividend <- get.obs.tsboot(obj=get(dividend.objname), type=r$dividend$type )
            divisor <- get.obs.tsboot(obj=get(divisor.objname), type=r$divisor$type )
            if( ( m1m1 && dividend$m1 == divisor$m1 ) ||
                  ( m1m2 && dividend$m1 == divisor$m2 ) ||
                    ( m2m1 && dividend$m2 == divisor$m1 ) ||
                      ( m2m2 && dividend$m2 == divisor$m2 ) ||
                        !(m1m1 || m1m2 || m2m1 || m2m2) )  {
              results <- compute.ratio(name=r$name,dividend=dividend,divisor=divisor,res=results)
            }
          }
        }
      }
    }
    save(results,file="results.Rdata")
  } else {
    load("results.Rdata")
  }

  # overview plots of data
  pdf(file="data_overview.pdf",onefile=T)
  for(name in names(quants)) {
    if(debug) print(name)
    ts.indices <- which( results$res.tsboot$name == name )
    indices <- which( results$res$name == name )
    plotwitherror(y=results$res$val[indices], dy=results$res$dval[indices], x=results$res$m12[indices],
                  main=name, xlab=expression(mu),ylab=name)
  }
  for(name in names(ratios)) {
    if(debug) print(name)
    ts.indices <- which( results$res.tsboot$name == name )
    indices <- which( results$res$name == name )
    plotwitherror(y=results$res$val[indices], dy=results$res$dval[indices], x=results$res$m12[indices],
                  main=gsub("_ov_","/",name), xlab=expression(mu),ylab=gsub("_ov_","/",name))
  }
  dev.off()

  # demonstrate functionality on mK/fK
  #mK_ov_fK.demo(results=results,debug=FALSE)
  #fDs_ov_fD.demo(results=results,debug=FALSE)
  
  # mDs_ov_fD.2D.demo(results)
  
  indices1 <- which( results$res$name == "m_Ds" )
  indices2 <- which( results$res$name == "f_Ds" )
  print( data.frame(val=results$res$val[indices1]/results$res$val[indices2],m1=results$res$m11[indices],m2=results$res$m12[indices]))
  
  
}

mDs_ov_fD.2D.demo <- function(results,debug=FALSE) {
  
  name <- "m_Ds_ov_f_Ds"
  
  # list of lists of length "number of bootstrap samples"
  dat <- vector( mode="list", length=length(which( results$res.tsboot$name == name ))/length(mass_comb$sc[,1]))
  
  # pre-allocate memory for list elements
  for( d in 1:length(dat) ) {
    dat[[d]] <- data.frame(z=vector(mode='numeric',length=16),x=vector(mode='numeric',length=16),y=vector(mode='numeric',length=16))
  }
  
  # loop over f_Ds mass combinations to construct list of data frames
  for( i in 1:length(mass_comb$sc[,1]) ){
    if(debug) print(mass_comb$sc[i,])
    ts.indices <- which( results$res.tsboot$name == name &
                         results$res.tsboot$m11 == mass_comb$sc[i,]$m1 &
                         results$res.tsboot$m12 == mass_comb$sc[i,]$m2 )
    for( j in 1:length(ts.indices) ) {
      dat[[j]][i,] <- c( z=results$res.tsboot$val[ts.indices[j]], x=mass_comb$sc[i,]$m1, y=mass_comb$sc[i,]$m2 )
    }
  }
  
  indices <- which( results$res$name == name )
  
  test <- vector(mode="list",length=50)
  
  for( i in 1:50 ) {
    test[[i]] <- dat[[i]]
  }

  fit <- fit_linear_2d( z=test, weights=(1/results$res$dval[indices])^2 )
  
  pars <- data.frame(a=vector(mode='numeric',length=length(fit$fit)),
  b=vector(mode='numeric',length=length(fit$fit)),c=vector(mode='numeric',length=length(fit$fit)))
  for( x in 1:length(fit$fit) ) {
    pars[x,] <- fit$fit[[x]]$m$getPars()
  }
  
  mu_s <- seq(0.022,0.03,length.out=50)
  dmu_s <- rep(0.002,50)
  mu_c <- 11.85*mu_s
  dmu_c <- 0.16*mu_s 
  
  pred <- extrapolate_2d(fit=fit,predx=mu_s,predy=mu_c,dpredx=dmu_s,dpredy=dmu_c)
  
  print(pred)
  
  pdf("mDs_ov_fDs_mu_s.pdf")
  plotwitherror( y=results$res$val[indices], x=results$res$m11[indices], dy=results$res$dval[indices], 
                 main=expression(m[D[s]]/f[D[s]]), ylab=expression(m[D[s]]/f[D[s]]), xlab=expression(mu[s]) )
  plot.confband(y=apply(X=pred$z,MARGIN=2,FUN=mean),dy=sqrt(apply(X=pred$z,MARGIN=2,FUN=sd)^2+apply(X=pred$dz,MARGIN=2,FUN=mean)^2),x=mu_s,col=rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4))
  legend(x=0.0212,y=8.64,legend=c("Measurements","mu_c = 11.85(16)*mu_s"),lty=c(0,1),col=c('black','black'),pch=c(1,NA))
  dev.off()
  
  # sqrt( sd^2 + apply(X=pred$dz,MARGIN=2,FUN=mean)^2)
  
  #print( c( mean(pars$a),sd(pars$a),mean(pars$b),sd(pars$b),mean(pars$c),sd(pars$c) ) )
}

fDs_ov_fD.demo <- function(results,debug=FALSE) {
  # extrapolate f_Ds/f_D in mu_c for mu_s=0.0238 
  predx = c(0.27,0.29,0.31,0.33)
  dpredx = rep(0.03,4)
  # find mu_s s.t. f_Ds/f_pi = x
  predy = c(1.17,1.18,1.19)
  dpredy = c(0.005,0.01,0.02)
  
  name <- "f_Ds_ov_f_D"
  
  ## TODO: prepare data into correct format with a little bit more logic!
  
  n.tsboot <- length( which( results$res.tsboot$name == name & results$res.tsboot$m11 == strange_masses[1] ) )/length(charm_masses)
  
  # choose indices for "f_Ds/f_D"
  indices <- which( results$res$name == name & results$res$m11 == strange_masses[1] )
  ts.indices <- array(dim=c(n.tsboot,length(charm_masses)))

  for( i in 1:length(charm_masses) ) {
    ts.indices[,i] <- which( results$res.tsboot$name == name & 
                             results$res.tsboot$m11 == strange_masses[1] & 
                             results$res.tsboot$m12 == charm_masses[i] )
  }

  legendlabels <- c("Measurements",expression(Input: mu[c]),expression(Input: f[D[s]]/f[D]))
  legendcols <- c('black','red','blue')
  legendsyms <- c(1,15,17)
  
  tsboot.dat <- cbind( results$res.tsboot$val[ts.indices[,1]], 
                results$res.tsboot$val[ts.indices[,2]], 
                results$res.tsboot$val[ts.indices[,3]],
                results$res.tsboot$val[ts.indices[,4]] )

  # compute standard confidence interval
  confband.x = seq(0.20,0.4,length.out=50)
  fit <- fit_linear_1d(y=tsboot.dat, x=charm_masses, type="lm", 
                              weights=1/results$res$dval[indices]^2, debug=debug)
  confband.y <- extrapolate_1d(fit=fit, predx=confband.x)
  
  # repeat fit using nls so that error propagation will work further below
  fit <- fit_linear_1d(y=tsboot.dat, x=charm_masses, type="nls", 
                                weights=1/results$res$dval[indices]^2, debug=debug)
 
 
  # without error propagation from dpredx
  newy <- extrapolate_1d(fit=fit, predx=predx)
  pdf(file=sprintf("%s_with_without_error.pdf",name))
  par(family="Times")
  plotwitherror(y=results$res$val[indices], dy=results$res$dval[indices], x=results$res$m12[ indices ],
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
  plotwitherror(y=results$res$val[indices], dy=results$res$dval[indices], x=results$res$m12[ indices ],
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

mK_ov_fK.demo <- function(results,debug=FALSE) {
  # extrapolate m_K/f_K in mu_s 
  predx = c(0.024,0.025,0.0255,0.0261)
  dpredx = rep(0.0003,4)
  # find mu_s s.t. m_K/f_K =  c(3.13,3.14,3.15)
  predy = c(3.12,3.14,3.15)
  dpredy = c(0.005,0.05,0.5)
  
  ## TODO: prepare data into correct format with a little bit more logic!
  
  n.tsboot <- length( which( results$res.tsboot$name == "m_K_ov_f_K" ) )/length(strange_masses)
  
  # choose indices for "m_K/f_K"
  indices <- which( results$res$name == "m_K_ov_f_K" )
  ts.indices <- array(dim=c(n.tsboot,length(strange_masses)))

  for( i in 1:length(strange_masses) ) {
    ts.indices[,i] <- which( results$res.tsboot$name == "m_K_ov_f_K" & results$res.tsboot$m12 == strange_masses[i] )
  }

  legendlabels <- c("Measurements",expression(Input: mu[s]),expression(Input: m[K]/f[K]))
  legendcols <- c('black','red','blue')
  legendsyms <- c(1,15,17)
  
  tsboot.dat <- cbind( results$res.tsboot$val[ts.indices[,1]], 
                results$res.tsboot$val[ts.indices[,2]], 
                results$res.tsboot$val[ts.indices[,3]],
                results$res.tsboot$val[ts.indices[,4]] )

  # compute standard confidence interval
  confband.x = seq(0.022,0.028,length.out=50)
  fit <- fit_linear_1d(y=tsboot.dat, x=strange_masses, type="lm", 
                              weights=1/results$res$dval[indices]^2, debug=debug)
  confband.y <- extrapolate_1d(fit=fit, predx=confband.x)
  
  # repeat fit using nls so that error propagation will work further below
  fit <- fit_linear_1d(y=tsboot.dat, x=strange_masses, type="nls", 
                                weights=1/results$res$dval[indices]^2, debug=debug)
 
 
  # without error propagation from dpredx
  newy <- extrapolate_1d(fit=fit, predx=predx)
  pdf(file="mKovfK_with_without_error.pdf")
  par(family="Times")
  plotwitherror(y=results$res$val[indices], dy=results$res$dval[indices], x=results$res$m12[ indices ],
                xlim=c(0.022,0.0265),main="without error propagation",xlab=expression(mu[s]),ylab=expression(m[K]/f[K]))
  plot.confband(y=apply(X=confband.y$y,MARGIN=2,FUN=mean),dy=apply(X=confband.y$y,MARGIN=2,FUN=sd),x=confband.x,col=rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4))
  # our predicted values
  plotwitherror(y=apply(X=newy$y,MARGIN=2,mean), dy=apply(X=newy$y,MARGIN=2,sd),x=predx,rep=T,col="red",pch=15)
  # without error propagation from dpredy
  newx <- solve_linear_1d( fit=fit, predy=predy )
  plotwitherror(x=apply(X=newx$x,MARGIN=2,FUN=mean),dx=apply(X=newx$x,FUN=sd,MARGIN=2),y=predy,rep=T,col='blue',pch=17)
  legend(x=0.022,y=3.2,pch=legendsyms,col=legendcols,legend=legendlabels)
  
  # with error propagation from dpredx 
  newy <- extrapolate_1d(fit=fit, predx=predx, dpredx=dpredx)
  plotwitherror(y=results$res$val[indices], dy=results$res$dval[indices], x=results$res$m12[ indices ],
                xlim=c(0.022,0.0265),main="with error propagation",xlab=expression(mu[s]),ylab=expression(m[K]/f[K]))
  plot.confband(y=apply(X=confband.y$y,MARGIN=2,FUN=mean),dy=apply(X=confband.y$y,MARGIN=2,FUN=sd),x=confband.x,col=rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4))
  # our predicted values
  # the first error is from propagating the error in predx, the second error is evaluated on the bootstrap sample
  plotwitherror(y=apply(X=newy$y,MARGIN=2,mean), dy=sqrt(apply(X=newy$dy,MARGIN=2,mean)^2+apply(X=newy$y,MARGIN=2,sd)^2),x=predx,dx=dpredx,rep=T,col="red",pch=15)
  # error propagation from dpredy
  newx <- solve_linear_1d( fit=fit, predy=predy, dpredy=dpredy )
  plotwitherror(x=apply(X=newx$x,MARGIN=2,FUN=mean),dx=sqrt(apply(X=newx$x,FUN=sd,MARGIN=2)^2+newx$dx^2),y=predy,dy=dpredy,rep=T,col='blue',pch=17)
  legend(x=0.022,y=3.2,pch=legendsyms,col=legendcols,legend=legendlabels)
  dev.off()
}

plot.confband <- function(x,y,dy,col) {
  xvals <- c(rev(x),x)
  yvals <- c(rev(y-dy),y+dy)
  polygon(x=xvals,y=yvals,col=col,border=NA)
  lines(x=x,y=y)
}

compute.quant <- function(name,dat,res) {
  n <- length(dat$dat)
  
  rval.tsboot <- data.frame( name=rep(name,n), val=dat$dat,
                m11=rep(dat$m1,n), m12=rep(dat$m2,n), m21=rep(NA,n), m22=rep(NA,n) )
  rval <- data.frame( name=name, val=mean(rval.tsboot$val), dval=sd(rval.tsboot$val),
                      m11=dat$m1, m12=dat$m2, m21=NA, m22=NA )
  
  return( append.result(rval,rval.tsboot,res) )
}

compute.ratio <- function(name,dividend,divisor,res) {
  n <- length(dividend$dat)
  if( n != length(divisor$dat) ) {
    stop("In compute.ratio, dividend and divisor do not have the same length!")
  }
  rval.tsboot <- data.frame( name=rep(name,n), val=dividend$dat/divisor$dat,
                m11=rep(dividend$m1,n), m12=rep(dividend$m2,n), m21=rep(divisor$m1,n), m22=rep(divisor$m2,n) )
  rval <- data.frame( name=name, val=mean(rval.tsboot$val), dval=sd(rval.tsboot$val),
                      m11=dividend$m1, m12=dividend$m2, m21=divisor$m1, m22=divisor$m2 )
                      
  return( append.result(rval,rval.tsboot,res) )
}

append.result <- function(val,val.tsboot,res) {
  if( length(res$res) == 0 ) {
    return( list(res=val,res.tsboot=val.tsboot) )
  } else {
    return( list(res=rbind(res$res,val), res.tsboot=rbind(res$res.tsboot,val.tsboot) ) )
  }
}

get.obs.tsboot <- function(obj,type) {
  rval <- list( m1=obj$mu1, m2=obj$mu2, dat=c() )
  if( type == "fps" ) {
    rval$dat <- obj$fps.tsboot
  } else if ( type == "mps" ) {
    rval$dat <- obj$opt.tsboot[1,]
  } else {
    stop("In 'load.obs.tsboot': type '", type,"' is not known. Exiting!")
  }
  return(rval)
}
