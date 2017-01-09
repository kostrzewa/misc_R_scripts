source("~/code/R/misc_R_scripts/prop_error/unary_operations.R")
source("~/code/R/misc_R_scripts/prop_error/binary_operations.R")

plot_mpi0_deltam <- function(mpi0file,r0file,splitfile,n.sim=1500,n.predict=500) {
  pt.col <- "darkblue"
  ln.col <- "grey45"

  mpi0dat <- read.table(mpi0file,header=TRUE)
  r0dat <- read.table(r0file,header=T)
  splitdat <- read.table(splitfile,header=TRUE)

  mpi0c <- mpi0dat[mpi0dat$name=="m_pi0c",]
  mpi <- mpi0dat[mpi0dat$name=="m_pi",]
  mpisq <- compute_square(mpi)
  mpi0csq <- compute_square(mpi0c)
  mpidiff <- compute_difference(mpisq,mpi0csq)

  r0sq <- compute_square(r0dat)
  r04 <- compute_square(r0sq)

  mpidiffr04 <- compute_product(mpidiff,r04) 
  mpisqr0sq <- compute_product(mpisq,r0sq)

  tikzfiles <- tikz.init(basename="mpi_mpi0_diff",width=3.5,height=3.7,lwdUnit=0.7)

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

  
  # chiral extrapolation of (m_pi^2-m_pi0^2)
  simdata <- data.frame(matrix(ncol=nrow(splitdat),nrow=n.sim))
  for(k in 1:nrow(splitdat)){
    simdata[,k] <- rnorm(mean=splitdat$val[k],sd=splitdat$dval[k],n=n.sim)
  }
  
  # constant extrapolation 
  sim.cfit <- apply(X=simdata,MARGIN=1,
                    FUN=function(x){
                      weighted.mean(x=x,w=1/splitdat$dval^2)
                    } )
 

  sim.lfit <- apply(X=simdata,MARGIN=1,
                    FUN=function(x){
                      lm(val~mu,data=data.frame(mu=splitdat$mu,val=x),weights=1/splitdat$dval^2)
                    } )

  newdat <- data.frame(mu=c(0.0009,seq(-0.002,0.008,length.out=n.predict-1)))

  pred.lfit <- lapply(X=sim.lfit,
                      FUN=function(x){
                        predict(x,newdata=newdat)
                      } )
  pred.lfit <- t(array(unlist(pred.lfit),dim=c(n.predict,n.sim)))
  pred.lfit <- t(apply(X=pred.lfit,MARGIN=2,FUN=quantile,probs=c(0.1573,0.5,0.8427)))



  # prepare plot for mass splitting 
  plot(x=splitdat$mu,y=splitdat$val,
       main="",xlab="",ylab="",yaxt='n',
       xlim=c(0,max(splitdat$mu)),
       ylim=c(-0.003,max(splitdat$val+splitdat$dval)),las=1,type='n')
  axis(side=2,at=seq(-0.003,0.013,0.0015),las=1)

  lims <- par('usr')
  width <- abs(lims[2]-lims[1])
  height <- abs(lims[4]-lims[3])
  
  #rect(xleft=lims[1],xright=lims[2],ytop=mean(sim.cfit)+sd(sim.cfit),ybottom=mean(sim.cfit)-sd(sim.cfit),col="#e41a1c66",border=NA)
  abline(h=mean(sim.cfit),lty=2)

  poly.x <- c(newdat$mu[-1],rev(newdat$mu[-1]))
  poly.y <- c(pred.lfit[-1,3],rev(pred.lfit[-1,1]))
  polygon(x=poly.x,y=poly.y,border=NA,col="#377eb844")
  lines(x=newdat$mu[-1],pred.lfit[-1,2],lty=3)
  plotwitherror(rep=TRUE,x=0.0009,y=mean(sim.cfit),dy=sd(sim.cfit),pch=15,col="red")
  plotwitherror(rep=TRUE,x=0.0009,y=pred.lfit[1,2],dy=pred.lfit[1,3]-pred.lfit[1,2],mdy=pred.lfit[1,2]-pred.lfit[1,1],pch=15,col="blue")
     
  plotwitherror(x=splitdat$mu,y=splitdat$val,dy=splitdat$dval,
                pch=splitdat$pch,col=pt.col,rep=TRUE)

  text(x=lims[1]+0.0003,y=lims[4]-0.10*height,labels="$a^2\\left( M^2_{\\pi^\\pm} - M^2_{\\pi^0} \\right)$",adj=0)
  text(x=lims[2],y=lims[3]+0.08*height,labels="$a\\mu_\\ell$",adj=1)
  abline(h=0,lty=3,col=ln.col)
  legend(x=lims[1],y=lims[3]+0.22*height,bty='n',legend=c("$L=24$","$L=32$"),pch=unique(mpi0c$pch),adj=0,col=pt.col)

  # extract square of charged pion mass with error from data file  
  mpi0dat <- read.table(mpi0file,header=TRUE)
  mpisqphys <- mpi0dat[mpi0dat$name=="m_pi",]
  mpisqphys <- mpisqphys[mpisqphys$mu1==0.0009,]
  mpisqphys <- compute_square(mpisqphys)

  mpi0sqphys.cfit <- compute_difference(mpisqphys,list(val=mean(sim.cfit),dval=sd(sim.cfit)),name="mpi0sqphys.cfit")
  mpi0sqphys.lfit <- compute_difference(mpisqphys,list(val=pred.lfit[1,2],dval=pred.lfit[1,3]-pred.lfit[1,2]),name="mpi0sqphys.lfit")
  
  mpi0phys.cfit <- compute_sqrt(mpi0sqphys.cfit)
  mpi0phys.cfit <- compute_ratio(mpi0phys.cfit,list(val=0.092,dval=0.002))
  mpi0phys.cfit <- compute_product(mpi0phys.cfit,list(val=197.3,dval=0.0),name="mpi0phys.cfit(mev)")

  mpi0phys.lfit <- compute_sqrt(mpi0sqphys.lfit)
  mpi0phys.lfit <- compute_ratio(mpi0phys.lfit,list(val=0.092,dval=0.002))
  mpi0phys.lfit <- compute_product(mpi0phys.lfit,list(val=197.3,dval=0.0),name="mpi0phys.lfit(mev)")

  print(mpi0sqphys.cfit)
  print(mpi0sqphys.lfit)
  print(mpi0phys.cfit)
  print(mpi0phys.lfit)

  tikz.finalize(tikzfiles)
}
