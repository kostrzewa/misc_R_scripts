construct.empty.meson_obs <- function() {
  meson_obs <- list()
  meson_obs[[1]] <- list( name=vector(mode="character",length=0), texlabel=vector(mode="character", length=0),
                     m.val=vector(mode="numeric",length=0), m.sea=vector(mode="numeric",length=0),
                     mean=vector(mode="numeric", length=0), err=vector(mode="numeric", length=0),
                     boot=vector(mode="numeric", length=0) )
  attr(meson_obs,"class") <- c("meson_obs",class(meson_obs))
  meson_obs
}

select.meson_obs <- function(meson_obs,by,filter) {
  if(!any(class(meson_obs) == "meson_obs") {
    stop("select.meson_obs: meson_obs argument must be of class meson_obs!\n")
  rval <- list()
  for( i in l ) {
    if( i[[by]] == filter )
      rval[[(length(rval)+1)]] <- i
  }
  attr(rval,"class") <- class(meson_obs)
  rval
}

as.data.frame.meson_obs <- function(meson_obs) {
  if(!any(class(meson_obs) {
   stop("as.data.frame.meson_obs: meson_obs argument must be of class meson_obs!\n")
  }
}
