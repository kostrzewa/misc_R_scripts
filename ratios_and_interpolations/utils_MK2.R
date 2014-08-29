compute.quant_mk2 <- function(name,texlabel,dat) {
  list( name=name, texlabel=texlabel, m.sea=sort(unique(dat$m.sea)), m.val=sort(unique(dat$m.val)),
        mean=mean(dat$boot), err=sd(dat$boot), boot=dat$boot )
}

#TODO: reimagine computation of expressions

compute.expression_mk2 <- function(expr,envir,m1,m2){
  return( list( dat=eval(expr=expr,envir=envir), m1=m1, m2=m2 ) )
}

compute.ratio_mk2 <- function(name,texlabel,dividend,divisor) {
  n <- length(dividend$boot)
  if( n != length(divisor$boot) ) {
    stop(sprintf("In compute.ratio for %s, dividend and divisor do not have the same number of bootstrap samples!\n",name))
  }
  
  ratio <- dividend$boot/divisor$boot
  
  list( name=name, texlabel=texlabel, m.sea=sort(unique(c(dividend$m.sea),c(divisor$m.sea))),
        m.val=sort(unique(c(dividend$m.val,divisor$m.val))), mean=mean(ratio),
        err=sd(ratio), boot=ratio )
}


# funtion to extract bootstrap samples from observables created by
# hadron functions where 'obj' is an ojbect returned by 'matrixfit'
# and similar functions

# the 'type' argument determines which list element the bootrap samples
# are extracted from

# the m.sea argument adds some more information about the sea quark mass
# dependence
get.obs.tsboot_mk2 <- function(obj, m.sea, type) {
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
        
    
