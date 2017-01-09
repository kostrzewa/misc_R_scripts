source("~/code/R/misc_R_scripts/prop_error/binary_operations.R")

plot_nd_kaon <- function(datfile) {
  require("RColorBrewer")
  clrs <- brewer.pal(n=9,name="Set1")
  nd_kaon <- read.table(header=TRUE,file=datfile,stringsAsFactors=FALSE)
  print(nd_kaon)
  mpi <- data.frame(val=nd_kaon$mpi,dval=nd_kaon$dmpi)
  mK <- data.frame(val=nd_kaon$mK,dval=nd_kaon$dmK)
  w0 <- data.frame(val=nd_kaon$w0,dval=nd_kaon$dw0)
  gmor <- data.frame(val=nd_kaon$gmor,dval=nd_kaon$dgmor)
 
  mKw0 <- compute_product(mK,w0)
  mKw0sq <- compute_product(mKw0,mKw0)

  w0sq <- compute_product(w0,w0)

  gmorw0sq <- compute_product(gmor,w0sq)

  mpisqw0sq <- compute_product(compute_product(mpi,mpi),w0sq)

  tikzfiles <- tikz.init(basename="nd_kaon_csw_beta",width=2.7,height=2.7, sanitize=FALSE)
  plotwitherror(y=2*gmorw0sq$val, dy=2*gmorw0sq$dval, x=1:length(mKw0$val),
                las=1, ylab="$w_0^2 (2M_K^2 - M_\\pi^2)$", xlab="", main="", xaxt='n',
                col=clrs[1], pch=15 , ylim=c(0.31,0.44),cex=1.0)
  axis(side=1,at=1:length(mKw0$val),labels=nd_kaon$label)
  plotwitherror(y=mpisqw0sq$val,dy=mpisqw0sq$dval,x=1:length(mKw0$val),
                las=1,ylab="$(w_0 M_\\pi)^2$", pch=15 , col=clrs[2], xlab="", main="", xaxt='n',ylim=c(0.0,0.13),cex=1.0)
  axis(side=1,at=1:length(mKw0$val),labels=nd_kaon$label)
  tikz.finalize(tikzfiles)
}
