source("~/code/R/misc_R_scripts/prop_error/unary_operations.R")
source("~/code/R/misc_R_scripts/prop_error/binary_operations.R")

plot_mpi0_deltam <- function(mpi0file,r0file) {
  pt.col <- "darkblue"
  ln.col <- "grey45"

  mpi0dat <- read.table(mpi0file,header=TRUE)
  r0dat <- read.table(r0file,header=T)

  mpi0c <- mpi0dat[mpi0dat$name=="m_pi0c",]
  mpi <- mpi0dat[mpi0dat$name=="m_pi",]
  mpisq <- compute_square(mpi)
  mpi0csq <- compute_square(mpi0c)
  mpidiff <- compute_difference(mpisq,mpi0csq)

  r0sq <- compute_square(r0dat)
  r04 <- compute_square(r0sq)

  mpidiffr04 <- compute_product(mpidiff,r04) 
  mpisqr0sq <- compute_product(mpisq,r0sq)

  tikzfiles <- tikz.init(basename="mpi_mpi0_diff",width=3.5,height=3.7,lwdUnit=0.6)

  plotwitherror(x=mpisqr0sq$val,dx=mpisqr0sq$dval,y=mpidiffr04$val,dy=mpidiffr04$dval,
                main="",xlab="",ylab="",pch=mpi0c$pch,
                xlim=c(0,max(mpisqr0sq$val)),
                ylim=c(min(mpidiffr04$val)-2,2.2),las=1,
                col=pt.col)
  
  lims <- par('usr')
  width <- abs(lims[2]-lims[1])
  height <- abs(lims[4]-lims[3])
  text(x=-0.01,y=lims[4]-0.10*height,labels="$\\left( \\frac{M^2_{\\pi^\\pm} - M^2_{\\pi^{(0,c)}}}{a^2} \\right) r_0^4$",adj=0)
  text(x=lims[2],y=lims[3]+0.08*height,labels="$r_0^2 M^2_{\\pi^\\pm}$",adj=1)
  abline(h=0,lty=3,col=ln.col)
  legend(x=lims[2]-0.3*width,y=lims[4]-0.3*height,bty='n',legend=c("$L=24$","$L=32$","$L=48$"),pch=unique(mpi0c$pch),adj=0,col=pt.col)


  mpi0dat <- mpi0dat[mpi0dat$L %in% c(24,32),]
  mpi0c <- mpi0dat[mpi0dat$name=="m_pi0c",]
  mpi0 <- mpi0dat[mpi0dat$name=="m_pi0",]

  mpi0csq <- compute_square(mpi0c)
  mpi0sq <- compute_square(mpi0)

  mpi <- mpi0dat[mpi0dat$name=="m_pi",]
  mpisq <- compute_square(mpi)
  mpisqr0sq <- compute_product(mpisq,r0sq)

  mpidiff <- compute_difference(mpi0csq,mpi0sq)
  mpidiff <- compute_ratio(mpidiff,list(val=2,dval=0))
  mpidiffr04 <- compute_product(mpidiff,r04)

  plotwitherror(x=mpisqr0sq$val,dx=mpisqr0sq$dval,y=mpidiffr04$val,dy=mpidiffr04$dval,
                main="",xlab="",ylab="",pch=mpi0c$pch,
                xlim=c(0,max(mpisqr0sq$val)),
                ylim=c(-2,max(mpidiffr04$val)+2),las=1,
                col=pt.col)

  lims <- par('usr')
  width <- abs(lims[2]-lims[1])
  height <- abs(lims[4]-lims[3])
  text(x=-0.01,y=lims[4]-0.10*height,labels="$\\left( \\frac{M^2_{\\pi^{(0,c)}} - M^2_{\\pi^0}}{2a^2} \\right) r_0^4$",adj=0)
  text(x=lims[2],y=lims[3]+0.08*height,labels="$r_0^2 M^2_{\\pi^\\pm}$",adj=1)
  abline(h=0,lty=3,col=ln.col)
  legend(x=lims[1],y=lims[3]+0.4*height,bty='n',legend=c("$L=24$","$L=32$"),pch=unique(mpi0c$pch),adj=0,col=pt.col)

  tikz.finalize(tikzfiles)
}
