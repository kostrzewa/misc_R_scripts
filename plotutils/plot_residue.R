
plot_residue <- function(datfiles,lg,basename="residues",width=7,height=5,ylims=c(-20,3)) {
  fdat <- list()
  for(f in datfiles){
    fdat[[length(fdat)+1]] <- read.table(header=TRUE,file=f,stringsAsFactors=FALSE)[,5]
  }
  
  maxiter <- 0
  for(rsq in fdat){
    if(length(rsq)>maxiter){
      maxiter <- 1.1*length(rsq)
    }
  }
  
  require("RColorBrewer")
  clr <- rep(brewer.pal(name="Set1",n=9),2)
  
  tikzfiles <- tikz.init(basename=basename,width=width,height=height)
  # prepare plot area
  plot(NULL,
       xlab="",ylab="$\\|r\\|^2$",main="",las=1,tck=0.02,
         ylim=10^ylims,xlim=c(0,maxiter),log='y',yaxt='n',
        )
  # draw major and minor tick marks
  axis(side=2, labels=TRUE, at=10^(ylims[1]:(ylims[2])),tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^(ylims[1]:(ylims[2]-1))),tck=0.01)
  mtext(side=1, line=1.1, text="$N_\\mathrm{iter}$") 
   
  for( i in 1:length(fdat) ) {
    # in case there is some slight disagreement among the number of steps, we determine them for each monomial separately
    lines(fdat[[i]],lty=1,col=clr[i],lwd=3)
  }
  if(!missing(lg)){
    legend(lty=1,lwd=3,col=clr[1:length(fdat)],legend=lg,x=0.8*maxiter,y=10^ylims[2],bty='n')
  }
  
  tikz.finalize(tikzfiles)
}


