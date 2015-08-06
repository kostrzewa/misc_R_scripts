source("~/code/R/misc_R_scripts/plotutils/plot_multihist.R")

plot_force_histograms <- function(datfiles,lg,basename="forces_hist_comparison",width=5,height=3) {
  fdat <- list()
  mon <- list()
  fdat[[1]] <- read.table(header=TRUE,file=datfiles[1],stringsAsFactors=FALSE)
  mon <- unique(fdat[[1]]$mon)
  for( i in 2:length(datfiles) ){
    fdat[[length(fdat)+1]] <- read.table(header=TRUE,file=datfiles[i],stringsAsFactors=FALSE)
    tmon <- unique(fdat[[i]]$mon)
    if(any(mon!=tmon))
      stop("plot_force_histograms: monomials inconsistent!!")
  }
  
  for( m in mon ){
    dat <- list()
    for( f in fdat ){
      dat[[length(dat)+1]] <- f$max[f$mon==m]
    }
    if(!missing(lg)){
      plot_multihist(dat=dat,basename=sprintf("%s.%s",basename,m),lg=lg,width=width,height=height,main="",xlab="$\\max(F^2)$")
    } else {
      plot_multihist(dat=dat,basename=sprintf("%s.%s",basename,m),width=width,height=height,main="",xlab="$\\max(F^2)$")
    }
  }

}


