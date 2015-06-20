source("~/code/R/misc_R_scripts/plotutils/plot_list_vertical.R")

plot_highstat_list_vertical <- function(datfile,basename,labelpos,width=6,...) {
  dat <- list()
  tikzfiles <- NULL
  finalize <- TRUE
  if(length(datfile)>1) finalize <- FALSE
  for( i in 1:length(datfile) ){
    dat <- read.table(header=TRUE,file=datfile[i],stringsAsFactors=FALSE)
    labels <- list()
    if(i<2){
      labels[[1]] <- dat$label
      labels[[2]] <- sprintf("$\\tau_\\mathrm{int}=%.3f, P_\\mathrm{acc}=%.3f$",dat$P.tauint,dat$ar)
    } else {
      labels[[1]] <- sprintf("$P_\\mathrm{acc}=%.3f$",dat$ar)
      labels[[2]] <- sprintf("$\\tau_\\mathrm{int}=%.3f$",dat$P.tauint) 
    }
    rep <- FALSE
    if(i>1) rep <- TRUE
    if(missing(labelpos)){
      temp <- plot_list_vertical(basename=basename,
                         x=dat$P.value,dx=dat$P.dvalue,
                         labels=labels,
                         pch=dat$pch,clr=dat$clr,
                         item.height=0.3, rep=rep,
                         xlab="$\\left\\langle P \\right\\rangle$",finalize=finalize,...)
    } else {
      temp <- plot_list_vertical(basename=basename,
                         x=dat$P.value,dx=dat$P.dvalue,
                         labels=labels,labelpos=labelpos[i],
                         pch=dat$pch,clr=dat$clr,
                         item.height=0.3, rep=rep,
                         xlab="$\\left\\langle P \\right\\rangle$",finalize=finalize,...)
    }
    if(i<2) tikzfiles <- temp
  }
  if(length(datfile)>1){
    tikz.finalize(tikzfiles)
  }
}
