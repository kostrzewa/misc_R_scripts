plot_hmc_timings <- function(path=".", ntraj=1, basename="hmc_timings", slims=c(0,8000), flims=c(0,0.3), DDalphaAMG=TRUE){
  datfiles <- c("heatbath_time.dat","accept_time.dat","derivative_time.dat","eigenvalues_time.dat","ddalphaamg_setup_time.dat","qphix_packing_time.dat")
  components <- c("heatbath","accept","derivative","eigenvalues","ddalphaamg_setup","qphix_packing")
 
  #if(!DDalphaAMG){
  #  datfiles <- datfiles[1:(length(datfiles)-1)]
  #  components <- components[1:(length(components)-1)]
  #}

  monomials <- NULL
  {
    for( i_cpt in 1:length(datfiles) ){
      temp.dat <- read.table(sprintf("%s/%s",path,datfiles[i_cpt]),header=TRUE,stringsAsFactors=FALSE)
      monomials <- unique(c(monomials,temp.dat$monomial))
    }
  }
  
  # sort monomials in a sensible order, this results in a lexical ordering
  monomials <- sort(monomials,method="radix")
  
  totals <- array(double(0),dim=c(length(components),length(monomials)))
  
  row.names(totals) <- components
  colnames(totals) <- monomials

  for( i_cpt in 1:length(datfiles) ){
    dat <- read.table(sprintf("%s/%s",path,datfiles[i_cpt]),stringsAsFactors=FALSE,header=TRUE)
    for( i_mon in 1:length(monomials) ){
      totals[ i_cpt, i_mon ] <- sum( dat[ which(dat$monomial == monomials[i_mon]), 2] )
    }
  }

  # for test trajectories, we can run one integration step on the outermost
  # timescale (compared to N integration steps) and pass ntraj = 1/N
  # such that we simply have to scale the *derivative* by this factor
  # to get a good estimate for the total cost (does not include
  # DDalphaAMG setup updates)
  if(ntraj >= 1){
    totals <- totals/ntraj
  } else {
    totals["derivative",] <- totals["derivative",]/ntraj
  }
  print(totals)
  fracs <- totals/sum(totals)
  print(fracs)
  cat(sprintf("\nTotal trajectory time, averaged over %f trajectories: %f hours\n",ntraj,sum(totals/3600)))
  
  require("RColorBrewer")
  clrs <- brewer.pal(n=length(components),name="Dark2")
  
  pdf(sprintf("%s.pdf",basename),width=7,height=4.5)
  
  par(mar=c(5,10,4,2))
  bpos <- barplot(totals,col=clrs, names.arg=rep("",length(monomials)), 
                  horiz=TRUE,beside=FALSE,xlab="time [s]", xlim=slims)
  coords <- par()$usr
  legend(x="topright",bty='n',pch=15,col=clrs,legend=components,pt.cex=2)
  text(x=coords[1]-0.01*slims[2], y=bpos, labels=monomials, srt=0, xpd=TRUE, adj=c(1,0.5))
  
  par(mar=c(5,10,4,2))
  bpos <- barplot(fracs,col=clrs, names.arg=rep("",length(monomials)), 
                  horiz=TRUE,beside=FALSE,xlab="fraction of total time",xlim=flims)
  coords <- par()$usr
  legend(x="topright",bty='n',pch=15,col=clrs,legend=components,pt.cex=2)
  text(x=coords[1]-0.25*0.01, y=bpos, labels=monomials, srt=0, xpd=TRUE, adj=c(1,0.5))
  
  dev.off()

  return( list(totals=totals, fracs=fracs) )
}
