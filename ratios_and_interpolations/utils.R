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
      # this gets an object of class "fitrange" as produced by meson_2pt_fitrange_analysis
      obj <- get(qty$datanames[index],env)
      m.val <- sort(unique(c(obj[[1]]$q_masses$m1,obj[[1]]$q_masses$m2)))
      cat("m.val:",m.val,"\n")
      # this will be a data frame containing bootstrap samples in rows and fit ranges in columns
      t <- NULL
      # a data frame containing start and end points of the fitranges
      fr <- NULL
      # a vector containing the weights of the given fitrange
      w <- NULL
      # loop through the fitranges
      my.ptm <- proc.time()
      mclist <- lapply.overload(X=obj,FUN=function(x) {
        Q <- NULL
        bsamples <- NULL
        weight <- NULL
        if(qty$type=="fps" && !is.na(x$f)){
          bsamples <- x$f$t
          Q <- x$f$Q
        } else if (qty$type=="mps" && !is.na(x$M)){
          bsamples <- x$M$t
          Q <- x$M$Q
        }
        tfr <- data.frame(t1=x$t1,t2=x$t2)
        if(!is.null(bsamples)){
          weight <- Qw(Q)^2/sd(bsamples)^2
        }
        list(fr=tfr,w=weight,t=bsamples)
      } )
      if(debug) {
        cat("Time for data reshuffling:\n")
        print(proc.time()-my.ptm)
      }
      # extract data from list
      for(l in mclist){
        summary(l)
        if(is.null(l$t) || is.null(l$w)) next
        t <- cbind(t,l$t)
        fr <- rbind(fr,l$fr)
        w <- c(w,l$w)
      }
      rval[[length(rval)+1]] <- list(name=qty$name,m.val=m.val,m.sea=m.sea,fr=fr,t=t,w=w)
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
        t <- NULL
        fr <- NULL
        w <- NULL
        my.ptm <- proc.time()
        mclist <- lapply.overload(X=obj,FUN=function(x){
          Q <- NULL
          bsamples <- NULL
          weight  <- NULL
          if(!is.na(x$M) && !is.na(x$f)){
            if(r$dividend$type == "mps" && r$divisor$type == "fps"){
               bsamples <- x$M$t/x$f$t
               weight <- (Qw(x$M$Q)/sd(bsamples))^2
            } else {
               bsamples <- x$f$t/x$M$t
               weight <- (Qw(x$M$Q)/sd(bsamples))^2
            }
          }  
          tfr <- data.frame(t1=x$t1,t2=x$t2)
          list(fr=tfr,w=weight,t=bsamples)
        } )
        if(debug) {
          cat("Time for data reshuffling:\n")
          print(proc.time()-my.ptm)
        }
        # extract data from list
        for(l in mclist){
          summary(l)
          if(is.null(l$t) || is.null(l$w)) next
          t <- cbind(t,l$t)
          fr <- rbind(fr,l$fr)
          w <- c(w,l$w)
        }
        rval[[length(rval)+1]] <- list(name=r$name,m.val=m.val,m.sea=m.sea,fr=fr,t=t,w=w)
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
            m.val <- sort(unique(c(dividend.list[[1]]$q_masses$m1,dividend.list[[1]]$q_masses$m2,divisor.list[[1]]$q_masses$m1,divisor.list[[1]]$q_masses$m2)))
            t <- NULL
            w <- NULL
            fr <- NULL
            # we use mclapply to loop over the fitranges of the dividend while an internal for loop
            # will loop over those of the divisor to ensure that only matching fitranges are considered
            my.ptm <- proc.time()
            mclist <- lapply.overload(X=dividend.list,FUN=function(x){
              weight <- NULL
              bsamples <- NULL
              for(obj in divisor.list){
                dividend <- NULL
                divisor <- NULL
                if( (x$t1 != obj$t1) || (x$t2 != obj$t2) ) {
                  next
                } else {
                  if(r$dividend$type == "mps" && r$divisor$type == "mps"){
                    dividend <- x$M
                    divisor <- obj$M
                  } else if (r$dividend$type == "fps" && r$divisor$type == "mps"){
                    dividend <- x$f
                    divisor <- obj$M
                  } else if (r$dividend$type == "mps" && r$divisor$type == "fps"){
                    dividend <- x$M
                    divisor <- obj$f
                  } else {
                    dividend <- x$f
                    divisor <- obj$f
                  }
                  #if(debug){ print(summary(dividend)); readline("key") }
                  if(!(is.na(dividend) || is.na(divisor))){
                    bsamples <- dividend$t/divisor$t
                    weight <- Qw(dividend$Q)*Qw(divisor$Q)*(1/sd(bsamples))^2
                  }
                  # no need to complete this loop since we found the correct fitrange!
                  break
                }
              }
              list(t=bsamples,w=weight,fr=data.frame(t1=x$t1,t2=x$t2))
            } )
            if(debug) {
              cat("Time for data reshuffling:\n")
              print(proc.time()-my.ptm)
            }
            # extract data from list
            for(l in mclist){
              summary(l)
              if(is.null(l$t) || is.null(l$w)) next
              t <- cbind(t,l$t)
              fr <- rbind(fr,l$fr)
              w <- c(w,l$w)
            }
            rval[[length(rval)+1]] <- list(name=r$name,m.val=m.val,m.sea=m.sea,fr=fr,t=t,w=w)
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

add.fitrange.systematic.hadron_obs <- function(hadron_obs,fitrange.serr){
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
    w.quantiles <- weighted.quantile(x=apply(X=fitrange.serr[[i]]$t,MARGIN=2,FUN=mean,na.rm=TRUE),w=fitrange.serr[[i]]$w,probs=c(0,0.1573,0.5,0.8427,1),na.rm=TRUE)
    m.serr <- abs(w.quantiles[2]-w.quantiles[3])
    p.serr <- abs(w.quantiles[4]-w.quantiles[3])
    # THIS IS A TOTAL HACK, we replace the median with the mean, and store the median over the fitranges
    # as the mean because otherwise the changes would be too numerous 
    hadron_obs[[i]]$median <- hadron_obs[[i]]$mean
    hadron_obs[[i]]$serr <- c(m.serr,p.serr)
    hadron_obs[[i]]$mean <- w.quantiles[3]
    hadron_obs[[i]]$boot <- apply(X=fitrange.serr[[i]]$t,MARGIN=1,FUN=weighted.median,w=fitrange.serr[[i]]$w,na.rm=TRUE) 
    hadron_obs[[i]]$err <- sd(hadron_obs[[i]]$boot)
    hadron_obs[[i]]$fr <- fitrange.serr[[i]]
  }
  hadron_obs
}
