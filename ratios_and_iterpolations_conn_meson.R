library('propagate')

ratios_and_iterpolations_conn_meson <- function(debug=F) {
  # masses to be used in this analysis
  strange_masses <- c(0.0238,0.0245,0.0252,0.0259)
  charm_masses <- c(0.2822,0.294,0.3058,0.3176)
  light_masses <- c(0.0009)
  
  # combinations of these masses
  mass_comb <- list( ll=data.frame( m1=light_masses, m2=light_masses ),
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
 
  # read all the data files
  for( n in names ) {
    for( file in n ) {
      file <- sprintf("%s.Rdata",file)
      if(debug) cat("Loading" ,file,"\n")
      load(file)
    }
  }

  # single quantities
  quants <- NULL
  quants[["m_pi"]] <- list(name="m_pi", m1=light_masses, m2=light_masses, datanames=names$ll_c )
  quants[["m_K"]] <- list(name="m_K", m1=light_masses, m2=strange_masses, datanames=names$ls_c )
  quants[["m_D"]] <- list(name="m_D", m1=light_masses, m2=charm_masses, datanames=names$lc_c )
  quants[["m_Ds"]] <- list(name="m_Ds", m1=strange_masses, m2=charm_masses, datanames=names$sc_c )
  quants[["f_pi"]] <- list(name="f_pi", m1=light_masses, m2=light_masses, datanames=names$ll_c )
  quants[["f_K"]] <- list(name="f_K", m1=light_masses, m2=strange_masses, datanames=names$ls_c )
  quants[["f_D"]] <- list(name="f_D", m1=light_masses, m2=charm_masses, datanames=names$lc_c )
  quants[["f_Ds"]] <- list(name="f_Ds", m1=strange_masses, m2=charm_masses, datanames=names$sc_c )

  # the different ratios that we would like to compute
  ratios <- NULL
  ratios[["m_pi/f_pi"]] <- list( name="m_pi/f_pi", dividend=quants[["m_pi"]], divisor=quants[["f_pi"]],
                               dividend.type="mps", divisor.type="fps" )
  ratios[["m_K/f_K"]] <- list( name="m_K/f_K", dividend=quants[["m_K"]], divisor=quants[["f_K"]],
                               dividend.type="mps", divisor.type="fps" )
  ratios[["m_D/f_D"]] <- list( name="m_D/f_D", dividend=quants[["m_D"]], divisor=quants[["f_D"]],
                               dividend.type="mps", divisor.type="fps" )
  ratios[["m_Ds/f_Ds"]] <- list( name="m_Ds/f_Ds", dividend=quants[["m_Ds"]], divisor=quants[["f_Ds"]],
                               dividend.type="mps", divisor.type="fps" )

  ratios[["m_K/f_pi"]] <- list( name="m_K/f_pi", dividend=quants[["m_K"]], divisor=quants[["f_pi"]],
                               dividend.type="mps", divisor.type="fps" )
  ratios[["m_D/f_pi"]] <- list( name="m_D/f_pi", dividend=quants[["m_D"]], divisor=quants[["f_pi"]],
                               dividend.type="mps", divisor.type="fps" )
  ratios[["m_Ds/f_pi"]] <- list( name="m_Ds/f_pi", dividend=quants[["m_Ds"]], divisor=quants[["f_pi"]],
                               dividend.type="mps", divisor.type="fps" )

  ratios[["f_K/f_pi"]] <- list( name="f_K/f_pi", dividend=quants[["f_K"]], divisor=quants[["f_pi"]], 
                                dividend.type="fps", divisor.type="fps" )
  ratios[["f_D/f_pi"]] <- list( name="f_D/f_pi", dividend=quants[["f_D"]], divisor=quants[["f_pi"]], 
                                dividend.type="fps", divisor.type="fps" )
  ratios[["f_Ds/f_pi"]] <- list( name="f_Ds/f_pi", dividend=quants[["f_Ds"]], divisor=quants[["f_pi"]], 
                                dividend.type="fps", divisor.type="fps" )

  ratios[["f_D/f_K"]] <- list( name="f_D/f_K", dividend=quants[["f_D"]], divisor=quants[["f_K"]], 
                                dividend.type="fps", divisor.type="fps" )
  ratios[["f_Ds/f_K"]] <- list( name="f_Ds/f_K", dividend=quants[["f_Ds"]], divisor=quants[["f_K"]], 
                                dividend.type="fps", divisor.type="fps" )

  ratios[["f_Ds/f_D"]] <- list( name="f_Ds/f_D", dividend=quants[["f_Ds"]], divisor=quants[["f_D"]], 
                                dividend.type="fps", divisor.type="fps" )

  results <- list( res = data.frame(), res.tsboot = data.frame() )
  for( r in ratios ) {
    # there is a special case if the data of dividend and divisor are exactly the same,
    # we don't need x*x cases but only x!
    if(debug) cat("Computing ratio", r$name, "\n")
    if( all(r$dividend$datanames == r$divisor$datanames) ) {
      for( objname in r$dividend$datanames ) {
        dividend <- load.obs.tsboot(obj=get(objname), type=r$dividend.type )
        divisor <- load.obs.tsboot(obj=get(objname), type=r$divisor.type )
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
          dividend <- load.obs.tsboot(obj=get(dividend.objname), type=r$dividend.type )
          divisor <- load.obs.tsboot(obj=get(divisor.objname), type=r$divisor.type )
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


  # extrapolate f_Ds/f_K in mu_s 
  #predx = seq(0.01,0.04,length.out=200) 
  #ts.indices <- which( results$res.tsboot$name == "f_Ds/f_K" & results$res.tsboot$m12 == charm_masses[1] )
  #indices <- which( results$res$name == "f_Ds/f_K" & results$res$m12 == charm_masses[1] )
  #fit.tsboot <<- fit.1d.tsboot(y=results$res.tsboot$r[ts.indices], x=strange_masses)
  #newy <- predx.1d.tsboot(fit.tsboot=fit.tsboot, predx=predx)
  #plotwitherror(y=results$res$r[indices], dy=results$res$dr[indices], x=results$res$m11[ indices ],xlim=c(0.018,0.042))
  #plot.confband(y=apply(X=newy,MARGIN=2,FUN=mean),dy=apply(X=newy,MARGIN=2,FUN=sd),x=predx,col=rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4))

  # extrapolate m_K/f_K in mu_s 
  predx = seq(0.01,0.04,length.out=200)                                                                                                                      
  ts.indices <- which( results$res.tsboot$name == "m_K/f_K" )
  indices <- which( results$res$name == "m_K/f_K" )
  fit.tsboot <- fit.1d.tsboot(y=results$res.tsboot$r[ts.indices], x=strange_masses)
  newy <- predx.1d.tsboot(fit.tsboot=fit.tsboot, predx=predx)
  plotwitherror(y=results$res$r[indices], dy=results$res$dr[indices], x=results$res$m12[ indices ],xlim=c(0.022,0.0265))
  plot.confband(y=apply(X=newy,MARGIN=2,FUN=mean),dy=apply(X=newy,MARGIN=2,FUN=sd),x=predx,col=rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4))


  # find mu_s s.t. m_K/f_K =  c(3.13,3.14,3.15)
  predy = c(3.13,3.14,3.15)
  newx <- predy.1d.tsboot( fit.tsboot, predy )
  print(newx)

  plotwitherror(x=apply(X=newx,MARGIN=2,FUN=mean),dx=apply(X=newx,FUN=sd,MARGIN=2),y=predy,rep=T)
}

plot.confband <- function(x,y,dy,col) {
  xvals <- c(rev(x),x)
  yvals <- c(rev(y-dy),y+dy)
  polygon(x=xvals,y=yvals,col=col,border=NA)
  lines(x=x,y=y)
}

fit.1d.tsboot <- function(y,x,debug=F) {
  # calculate how many fits we have to carry out
  n <- length(y)/length(x)
  stride <- (0:(length(x)-1))*n
  fit <- NULL
  for( index in 1:n ) {
    indices <- index+stride
    dat <- data.frame(y=y[indices],x=x)
    fit[[index]] <- lm(y~x,data=dat)
  }
  return( list(fit=fit,n=n) )
} 

predx.1d.tsboot <- function(fit.tsboot,predx,dpredx) {
  # we want a prediction for some new values of x 
  newy <- array(dim=c(fit.tsboot$n,length(predx)))
  for(index in 1:fit.tsboot$n) { 
    newy[index,] <- predict(fit.tsboot$fit[[index]],newdata=data.frame(x=predx))
  }
  return(newy)
}

predy.1d.tsboot <- function(fit.tsboot,predy,dpredy) {
  newx <- array(dim=c(fit.tsboot$n,length(predy)))
  for(index in 1:fit.tsboot$n) {
    a <- fit.tsboot$fit[[index]]$coefficients[2]
    b <- fit.tsboot$fit[[index]]$coefficients[1]
    newx[index,] <- sapply(X=predy-b,FUN=solve,a=a) 
  }
  return(newx)
}

compute.ratio <- function(name,dividend,divisor,res) {
  n <- length(dividend$dat)
  if( n != length(divisor$dat) ) {
    stop("In compute.ratio, dividend and divisor do not have the same length!")
  }
  rval.tsboot <- data.frame( name=rep(name,n), r=dividend$dat/divisor$dat,
                m11=rep(dividend$m1,n), m12=rep(dividend$m2,n), m21=rep(divisor$m1,n), m22=rep(divisor$m2,n) )
  rval <- data.frame( name=name, r=mean(rval.tsboot$r), dr=sd(rval.tsboot$r),
                      m11=dividend$m1, m12=dividend$m2, m21=divisor$m1, m22=divisor$m2 )

  if( length(res$res) == 0 ) {
    return( list(res=rval,res.tsboot=rval.tsboot) )
  } else {
    return( list(res=rbind(res$res,rval), res.tsboot=rbind(res$res.tsboot,rval.tsboot) ) )
  }
}

append.result <- function(res.dframe,result) {
  if( length(res.dframe) == 0 ) {
    res.dframe <- result
  } else {
    res.dframe <- rbind(res.dframe,result)
  }
  return(res.dframe)
}

load.obs.tsboot <- function(obj,type) {
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
