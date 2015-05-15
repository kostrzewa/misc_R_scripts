# given extrapolations and extrapolations.sys from the ratios_and_interpolations
# analysis, extract the various values from the last row for each observable
# as this corresponds to the situation of interest

extract_extrapolations_listplot <- function(ext, ext.sys) {
  if( length(ext) != length(ext.sys) ){
    stop("extract_interpolations_listplot: ext and ext.sys should be of the same length!\n")
  }

  x <- NULL
  dx <- NULL
  mdx <- NULL
  labels <- NULL
  for(i in 1:length(ext)){
    nr <- nrow(ext[[i]])
    labels <- c(labels,ext[[i]]$texlabel[nr])
    x <- c(x,ext[[i]]$val[nr])
    dx <- rbind(dx,cbind(ext[[i]]$dval[nr],ext.sys[[i]]$serr[nr]))
    mdx <- rbind(mdx,cbind(ext[[i]]$dval[nr],ext.sys[[i]]$mserr[nr]))
  }
  list(labels=labels,x=x,dx=dx,mdx=mdx)
}


