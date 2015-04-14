compute.quant <- function(name,texlabel,dat) {
  list( name=name, texlabel=texlabel, m.sea=sort(unique(dat$m.sea)), m.val=sort(unique(dat$m.val)),
        mean=mean(dat$boot), err=sd(dat$boot), serr=rep(0,times=2), boot=dat$boot )
}

#TODO: reimagine computation of expressions
compute.expression <- function(expr,envir,m1,m2){
  return( list( dat=eval(expr=expr,envir=envir), m1=m1, m2=m2 ) )
}

compute.ratio <- function(name,texlabel,dividend,divisor) {
  n <- length(dividend$boot)
  if( n != length(divisor$boot) ) {
    stop(sprintf("In compute.ratio for %s, dividend and divisor do not have the same number of bootstrap samples!\n",name))
  }
  
  ratio <- dividend$boot/divisor$boot
  
  list( name=name, texlabel=texlabel, m.sea=sort(unique(c(dividend$m.sea),c(divisor$m.sea))),
        m.val=sort(unique(c(dividend$m.val,divisor$m.val))), mean=mean(ratio), serr=rep(0,times=2),
        err=sd(ratio), boot=ratio )
}

Qw <- function(Q) {
  (1-2*abs(Q-0.5))
}

# this function is the evil child of unholy black magic, look at it only with appropriate protection
# and keep in mind that something like it should never be designed again
# it does the job, but there is nothing about it which is in any shape or form
# easy to understand or to maintain

compute.fitrange_systematic <- function(quants,ratios,m.sea,debug=FALSE,mc=TRUE,probs=c(0,0.1573,0.5,0.8427,1)) {
  lapply.overload <- lapply
  if(mc){
    require("parallel")
    lapply.overload <- mclapply
  }
  # load all data into the environment local to this function (this doesn't work using mclapply unfortunately...)
  env <- environment(NULL)
  temp <- lapply(X=quants,FUN=function(x) {
    for( dataname in x$datanames )  {
      if(!exists(dataname,envir=env)){
        file <- sprintf("%s.Rdata",dataname)
        if(debug) cat("Loading" ,file,"\n")
        load(file,envir=env)
      }
    }
  } )
  
  rval <- list()
  for( qty in quants ) {
    if(debug) cat("Computing fitrange systematic for", qty$name, "\n")
    # loop "over mass combinations" (effectively)
    temp <- list()
    for( index in 1:length(qty$datanames) ) {
      obj <- get(qty$datanames[index],env)
      m.val <- sort(unique(c(obj[[1]]$q_masses$m1,obj[[1]]$q_masses$m2)))
      mclist <- lapply.overload(X=obj,FUN=function(x) {
        if(!(is.na(x$M) || is.na(x$f))){
          val <- NULL
          dval <- NULL
          Q <- NULL
          if(qty$type=="fps"){
            val <- mean(x$f$t)
            dval <- sd(x$f$t)
            Q <- x$f$Q
          }else if(qty$type=="mps"){
            val <- mean(x$M$t)
            dval <- sd(x$M$t)
            Q <- x$M$Q
          }
          return(data.frame(t1=x$t1, t2=x$t2, val=val, dval=dval, w=(Qw(Q)*(1/dval))^2 ))
        } else {
          return(NA)
        }
      } )
      # reformat mclist into some kind of data frame
      mcdf <- NULL
      for( df in mclist ){
        if( !any(class(df)=="try-error") && !any(is.na(df)) ){
          mcdf <- rbind(mcdf,df)
        } else {
          cat("compute.fitrange_systematic: qty try-error\n")
        }
      }
      # store in a sensible format
      rval[[length(rval)+1]] <- list(name=qty$name, m.val=m.val, m.sea=m.sea, df=mcdf, fr.sys=weighted.quantile(x=mcdf$val,w=mcdf$w,na.rm=TRUE,probs=probs))
    }
  }

  for( r in ratios ) {
   # there is a special case if the data of dividend and divisor are exactly the same,
   # we don't need x*y cases but only x! This would be the case for a meson mass
   # and the related decay constant, for example
   # in this case we don't need to care about the masses
    if(debug) cat("Computing fitrange systematic for", r$name, "\n")
    temp <- list()
    if( all(r$dividend$datanames == r$divisor$datanames) ) {
      # loop over mass combinations (effectively)
      for( index in 1:length(r$dividend$datanames) ) {
        obj <- get(r$dividend$datanames[index],env)
        m.val <- sort(unique(c(obj[[1]]$q_masses$m1,obj[[1]]$q_masses$m2)))
        mclist <- lapply.overload(X=obj,FUN=function(x){
          val <- NA
          dval <- NA
          w <- NA 
          if(!is.na(x$M) && !is.na(x$f)){
            if(r$dividend$type == "mps" && r$divisor$type == "fps"){
              val <- mean(x$M$t/x$f$t)
              dval <- sd(x$M$t/x$f$t)
              w <- (Qw(x$M$Q)*(1/dval))^2
            } else {
              val <- mean(x$f$t/x$M$t)
              dval <- sd(x$f$t/x$M$t)
              w <- (Qw(x$M$Q)*(1/dval))^2
            }
          }
        data.frame(t1=x$t1, t2=x$t2, val=val, dval=dval, w=w )
        } )
      
        # reformat mclist into some kind of data frame
        mcdf <- NULL
        for( df in mclist ){
          if( !any(class(df)=="try-error") && !any(is.na(df)) ){
            mcdf <- rbind(mcdf,df)
          } else {
            cat("compute.fitrange_systematic: simple ratio try-error\n")
          }
        }
        # store in a sensible format
        rval[[length(rval)+1]] <- list(name=r$name, m.val=m.val, m.sea=m.sea, df=mcdf, fr.sys=weighted.quantile(x=mcdf$val,w=mcdf$w,na.rm=TRUE,probs=probs))
      } 
    } else {
      # when any pair of mass vectors from the two observables are the same
      # we need to ensure that "off-diagonal" elements, where say the strange
      # mass for the divisor and dividend are different, are skipped!
      # this would be the case, for instance for a ratio involving the kaon
      # and the D_s meson
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
          dividend.list <- get(dividend.objname,env)
          divisor.list <- get(divisor.objname,env)

          dividend.m.val <- sort(unique(c(dividend.list[[1]]$q_masses$m1,dividend.list[[1]]$q_masses$m2)))
          divisor.m.val <- sort(unique(c(divisor.list[[1]]$q_masses$m1,divisor.list[[1]]$q_masses$m2)))

          # note the second condition: if no two mass vectors are the same
          # TODO: solve the first if condition...
          if( any(as.vector(outer(dividend.m.val,divisor.m.val,'=='))) ||
              !(m1m1 || m1m2 || m2m1 || m2m2) )  
          {
            if(length(dividend.list)!=length(divisor.list)){
              msg <- sprintf("Length mismatch between dividend %s and divisor %s!\n",dividend$name,divisor$name)
              stop(msg)
            }

            m.val <- sort(unique(c(dividend.list[[1]]$q_masses$m1,dividend.list[[1]]$q_masses$m2,divisor.list[[1]]$q_masses$m1,divisor.list[[1]]$q_masses$m2)))
            mcdf <- NULL
            # cannot use mclapply here unfortunately....
            for(i in 1:length(dividend.list)){
              t1 <- dividend.list[[i]]$t1
              t2 <- dividend.list[[i]]$t2
              divisor <- NULL
              dividend <- NULL
              if(r$dividend$type == "mps" && r$divisor$type == "mps"){
                dividend <- dividend.list[[i]]$M
                divisor <- divisor.list[[i]]$M
              } else if (r$dividend$type == "fps" && r$divisor$type == "mps"){
                dividend <- dividend.list[[i]]$f
                divisor <- divisor.list[[i]]$M
              } else if (r$dividend$type == "mps" && r$divisor$type == "fps"){
                dividend <- dividend.list[[i]]$M
                divisor <- divisor.list[[i]]$f
              } else {
                dividend <- dividend.list[[i]]$f
                divisor <- divisor.list[[i]]$f
              }
              #if(debug){ print(summary(dividend)); readline("key") }
              if(!(is.na(dividend) || is.na(divisor))){
                val <- mean(dividend$t/divisor$t)
                dval <- sd(dividend$t/divisor$t)
                w <- Qw(dividend$Q)*Qw(divisor$Q)*(1/dval)^2
                mcdf <- rbind(mcdf,data.frame(t1=t1, t2=t2, val=val, dval=dval, w=w ))
              }
            }
            # store in a sensible format
            rval[[length(rval)+1]] <- list(name=r$name, m.val=m.val, m.sea=m.sea, df=mcdf, fr.sys=weighted.quantile(x=mcdf$val,w=mcdf$w,na.rm=TRUE,probs=probs))
          }
        }
      }
    }
  }
  class(rval) <- c(class(rval),"fitrange.serr")
  rval
}

# loop through quants, ratios and expressions and compute expectation values and errors from
# the bootstrap samples which are stored in objects in "envir"
compute.hadron_obs <- function(envir,quants,ratios,expressions,m.sea,debug=F) {
  hadron_obs <- construct.empty.hadron_obs()

  for( qty in quants ) {
    for( objname in qty$datanames ) {
      dat <- get.obs.tsboot(obj=get(objname,envir), m.sea=m.sea, type=qty$type)
      hadron_obs[[length(hadron_obs)+1]] <- compute.quant( name=qty$name, texlabel=qty$texlabel, dat=dat )
    }
  }

  if(!missing(ratios)) {
    for( r in ratios ) {
      # there is a special case if the data of dividend and divisor are exactly the same,
      # we don't need x*x cases but only x! This would be the case for a meson mass
      # and its decay constant, for example
      if(debug) cat("Computing ratio", r$name, "\n")
      if( all(r$dividend$datanames == r$divisor$datanames) ) {
        for( objname in r$dividend$datanames ) {
          dividend <- get.obs.tsboot(obj=get(objname,envir), type=r$dividend$type, m.sea=m.sea )
          divisor <- get.obs.tsboot(obj=get(objname,envir), type=r$divisor$type, m.sea=m.sea )
          hadron_obs[[length(hadron_obs)+1]] <- compute.ratio(name=r$name,texlabel=r$texlabel,dividend=dividend,divisor=divisor)
        }
      } else {
        # when any pair of mass vectors from the two observables are the same
        # we need to ensure that "off-diagonal" elements, where say the strange
        # mass for the divisor and dividend are different, are skipped!
        # this would be the case, for instance for a ratio involving the kaon
        # and the D_s meson
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
            dividend <- get.obs.tsboot(obj=get(dividend.objname,envir), type=r$dividend$type, m.sea=m.sea )
            divisor <- get.obs.tsboot(obj=get(divisor.objname,envir), type=r$divisor$type, m.sea=m.sea )

            # note the second condition: if no two mass vectors are the same
            if( any(as.vector(outer(divisor$m.val,dividend$m.val,'=='))) ||
                !(m1m1 || m1m2 || m2m1 || m2m2) )  {
              hadron_obs[[length(hadron_obs)+1]] <- compute.ratio(name=r$name,texlabel=r$texlabel,dividend=dividend,divisor=divisor)
            }
          }
        }
      }
    }
  }

  if(!missing(expressions)) {
    stop("compute.hadron_obs: Computation of expressions has not been implemented yet!")
    # the expressions are not quite ready yet,,,
#    # ratios from expressions
#    cat("Computing ratio f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD\n")
#    for( i in 1:length(mass_comb$sc[,1]) ) {
#      fDs_indices <- which(hadron_obs$val.tsboot$name == "f_Ds" &
#                          hadron_obs$val.tsboot$m11 == mass_comb$sc[i,]$m1 &
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      mDs_indices <- which(hadron_obs$val.tsboot$name == "m_Ds" &
#                          hadron_obs$val.tsboot$m11 == mass_comb$sc[i,]$m1 &
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      fD_indices <- which(hadron_obs$val.tsboot$name == "f_D" &
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      mD_indices <- which(hadron_obs$val.tsboot$name == "m_D" &
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#
#      dividend <- compute.expression(expr=expression(a*sqrt(b)),
#                    envir=data.frame( a=hadron_obs$val.tsboot$val[fDs_indices], b=hadron_obs$val.tsboot$val[mDs_indices] ),
#                    m1=mass_comb$sc[i,]$m1, m2=mass_comb$sc[i,]$m2 )
#      divisor <- compute.expression(expr=expression(a*sqrt(b)),
#                    envir=data.frame( a=hadron_obs$val.tsboot$val[fD_indices], b=hadron_obs$val.tsboot$val[mD_indices] ),
#                    m1=0.0009, m2=mass_comb$sc[i,]$m2 )
#
#      hadron_obs <- compute.ratio(name="f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD",texlabel="$f_{D_s}\\sqrt{m_{D_s}}/f_D\\sqrt{m_D}$",
#                               dividend=dividend, divisor=divisor, res=hadron_obs )
#    }

#    cat("Computing ratio f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd\n")
#    for( i in 1:length(mass_comb$sc[,1]) ) {
#      fDs_indices <- which(hadron_obs$val.tsboot$name == "f_Ds" &
#                          hadron_obs$val.tsboot$m11 == mass_comb$sc[i,]$m1 &
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      mDs_indices <- which(hadron_obs$val.tsboot$name == "m_Ds" &
#                          hadron_obs$val.tsboot$m11 == mass_comb$sc[i,]$m1 &
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      fD_indices <- which(hadron_obs$val.tsboot$name == "f_D" &
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      mD_indices <- which(hadron_obs$val.tsboot$name == "m_D" &
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#
#      dividend <- compute.expression(expr=expression(a*b^2),
#                    envir=data.frame( a=hadron_obs$val.tsboot$val[fDs_indices], b=hadron_obs$val.tsboot$val[mDs_indices] ),
#                    m1=mass_comb$sc[i,]$m1, m2=mass_comb$sc[i,]$m2 )
#      divisor <- compute.expression(expr=expression(a*b^2),
#                    envir=data.frame( a=hadron_obs$val.tsboot$val[fD_indices], b=hadron_obs$val.tsboot$val[mD_indices] ),
#                    m1=0.0009, m2=mass_comb$sc[i,]$m2 )
#
#      hadron_obs <- compute.ratio(name="f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd",texlabel="$f_{D_s}m_{D_s}^2/f_D{m_D}^2$",
#                               dividend=dividend, divisor=divisor, res=hadron_obs )
#    }
  }
  
  hadron_obs
}


# funtion to extract bootstrap samples from observables created by
# hadron functions where 'obj' is an ojbect returned by 'matrixfit'
# and similar functions

# the 'type' argument determines which list element the bootrap samples
# are extracted from

# the m.sea argument adds some more information about the sea quark mass
# dependence
get.obs.tsboot <- function(obj, m.sea, type) {
  rval <- list( m.val=sort(unique(c(obj$mu1,obj$mu2))), m.sea=m.sea, boot=c() )
  if( type == "fps" ) {
    rval$boot <- obj$fps.tsboot
  } else if ( type == "mps" ) {
    rval$boot <- obj$opt.tsboot[1,]
  } else {
    stop("In 'get.obs.tsboot': type '", type,"' is not known. Exiting!")
  }
  return(rval)
}

expand.grid.unique <- function(a,b) {
  # a and b are vectors or data frames
  if( any( c(length(dim(a)),length(dim(b)))>2 ) ) {
    stop("expand.grid.unique: one of the arguments has more than 2 dimensions!")
  }
  
  dat <- list()
  dat[[1]] <- a
  dat[[2]] <- b
  
  # we deconstruct a and b by column and only keep unique entries
  # in addition we sort the resulting vectors and keep them in a numbered
  # list
  temp <- list()
  for( i in 1:2 ) {
    if( length(dim(dat[[i]])) == 2 ) {
      for ( c.i in 1:(dim(dat[[i]])[2]) ) {
        temp[[length(temp)+1]] <- sort(unique( dat[[i]][,c.i] ))
      }  
    } else {
      temp[[length(temp)+1]] <- sort(unique(dat[[i]]))
    }
  }
  
  # we do a double loop through the list to find any duplicated vectors
  dups.idx <- c()
  for( i in 1:length(temp) ) {
    if( (i+1) < length(temp) ) { 
      for( j in (i+1):length(temp) ) {
        if( all( temp[[i]] == temp[[j]] ) ) {
          dups.idx <- c(dups.idx,j)
        } else if ( any( temp[[i]] == temp[[j]] ) ) {
          stop("expand.grid.unique: after separating and sorting, a and b have common entries but no common columns!\n")
        }
      }
    }
  }
  
  unique.cols <- list()
  for( i in 1:length(temp) ) {
    if( !any( dups.idx == i ) ) {
      unique.cols[[length(unique.cols)+1]] <- temp[[i]]
    }
  }

  expand.grid(unique.cols)  
}

add.fitrange.serr.hadron_obs <- function(hadron_obs,fitrange.serr){
  if(!any(class(hadron_obs)=="hadron_obs")){
    stop("add.fitrange.serr.hadron_obs: hadron_obs argument must be of class hadron_obs!")
  }
  if(!any(class(fitrange.serr)=="fitrange.serr")){
    stop("add.fitrange.serr.hadron_obs: fitrange.serr argument must be of class fitrange.serr!")
  }
  if(length(hadron_obs)!=length(fitrange.serr)){
    stop("add fitrange.serr.hadron_obs: fitrange.serr and hadron_obs do not have the same length!")
  }
  for(i in 1:length(hadron_obs)){
    if(hadron_obs[[i]]$m.val != fitrange.serr[[i]]$m.val){
      msg <- sprintf("add.fitrange.serr.hadron_obs: hadron_obs and fitrange.serr %d seem to have different m.val: %s and %s",i,hadron_obs[[i]]$m.val,fitrange.serr[[i]]$m.val)
      stop(msg)
    }
    m.serr <- abs(fitrange.serr[[i]]$fr.sys[2]-fitrange.serr[[i]]$fr.sys[3])
    p.serr <- abs(fitrange.serr[[i]]$fr.sys[4]-fitrange.serr[[i]]$fr.sys[3])
    hadron_obs[[i]]$serr <- c(m.serr,p.serr)
    hadron_obs[[i]]$median <- fitrange.serr[[i]]$fr.sys[3]
  }
  hadron_obs
}
