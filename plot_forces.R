plot_forces <- function(datfile,basename="forces",width=10,height=6,ylims=c(1e-5,1e3)) {
  fdat <- read.table(header=TRUE,file=datfile,stringsAsFactors=FALSE)
  mon <- unique(fdat$mon)
  
  nsteps <- length(which(fdat$mon==mon[1]))

  require("RColorBrewer")
  clr <- brewer.pal(name="Paired",n=12)

  tikzfiles <- tikz.init(basename=basename,width=width,height=height)
  # prepare plot area
  plot(NULL,
       xlab="$t_\\mathrm{MD}$",ylab="$F^2$",main="",las=1,tck=0.02,
       ylim=ylims,xlim=c(1,nsteps)
      )
  axis(side=2, labels=FALSE, at=outer(2:9,ylims))
  tikz.finalize(tikzfiles)
}


