source("~/code/R/misc_R_scripts/prop_error/binary_operations.R")

plot_pion_splitting <- function(datfile,lg,basename="pion_splitting_csw_beta",width=3.5,height=3.5,...) {
  require("RColorBrewer")
  clrs <- brewer.pal(n=9,name="Set1")
  ps <- read.table(header=TRUE,file=datfile,stringsAsFactors=FALSE,fill=TRUE)
 
  w0 <- data.frame(val=ps$w0,dval=ps$dw0)
  w0sq <- compute_product(w0,w0)
  w0ft <- compute_product(w0sq,w0sq)

  mpisqdiffw0sq <- compute_product(w0sq, data.frame(val=ps$mpi0csq_m_mpisq,dval=ps$dmpi0csq_m_mpisq))
  mpi0sqdiffw0sq <- compute_product(w0sq, data.frame(val=ps$mpisq_m_mpi0dsq,dval=ps$dmpisq_m_mpi0dsq))
  
  mpisqdiffw0ft <- compute_product(w0ft, data.frame(val=ps$mpi0csq_m_mpisq,dval=ps$dmpi0csq_m_mpisq))
  mpi0sqdiffw0ft <- compute_product(w0ft, data.frame(val=ps$mpi0csq_m_mpi0dsq/2,dval=ps$dmpi0csq_m_mpi0dsq/2))

  tikzfiles <- tikz.init(basename=basename,width=width,height=height,lwdUnit=0.8)
  plotwitherror(y=mpisqdiffw0sq$val,dy=mpisqdiffw0sq$dval,x=1:length(ps$L),
                pch=ps$pch,col=ps$clr,xaxt='n',ylab="$w_0^2( M_{\\pi^{(0,c)}}^2 - M_{\\pi^\\pm}^2 )$",xlab="",las=1, 
                ylim=c(0.035,0.095),xlim=c(1,length(ps$L)),... )
  if(!missing(lg)) {
    legend(x="topright",legend=lg$labels,pch=lg$pch,col=lg$col,cex=0.8,pt.cex=1.0,bty='n')
  }
  axis(side=1,at=1:length(ps$L),labels=ps$label)
  
  plotwitherror(y=mpi0sqdiffw0sq$val,dy=mpi0sqdiffw0sq$dval,x=1:length(ps$L),
                pch=ps$pch,col=ps$clr2,xaxt='n',ylab="$w_0^2( M_{\\pi^\\pm}^2 - M_{\\pi^0}^2 )$",xlab="",las=1,
                ylim=c(0.005,0.060),xlim=c(1,length(ps$L)),... )
  if(!missing(lg)) {
    legend(x="topleft",legend=lg$labels,pch=lg$pch,col=lg$col,bty='n',cex=0.8,pt.cex=1.0)
  }
  axis(side=1,at=1:length(ps$L),labels=ps$label)
  
  plotwitherror(y=ps$mpi0c_m_mpi_rel,dy=ps$dmpi0c_m_mpi_rel,x=1:length(ps$L),
                pch=ps$pch,col=ps$clr,xaxt='n',ylab="$\\frac{M_{\\pi^{(0,c)}} - M_{\\pi^\\pm}}{M_{\\pi^\\pm}}$",
                xlab="",las=1,ylim=c(0.1,0.8),xlim=c(1,length(ps$L)),...)
  if(!missing(lg)) {
    legend(x="bottomleft",legend=lg$labels,pch=lg$pch,col=lg$col,bty='n',cex=0.8,pt.cex=1.0)
  }
  axis(side=1,at=1:length(ps$L),labels=ps$label)
  
  plotwitherror(y=ps$mpi_m_mpi0d_rel,dy=ps$dmpi_m_mpi0d_rel,x=1:length(ps$L),
                pch=ps$pch,col=ps$clr2,xaxt='n',ylab="$\\frac{M_{\\pi^\\pm} - M_{\\pi^0}}{M_{\\pi^\\pm}}$",
                xlab="",las=1,ylim=c(0.0,0.6),xlim=c(1,length(ps$L)),... )
  if(!missing(lg)) {
    legend(x="topright",legend=lg$labels,pch=lg$pch,col=lg$col,cex=0.8,pt.cex=1.0,bty='n')
  }
  axis(side=1,at=1:length(ps$L),labels=ps$label)

  tikz.finalize(tikzfiles)
}
