plot_BSM_etascan <- function(datfile,debug=FALSE) {
  # kappa eta conf iters res
  iters <- read.table(datfile,header=TRUE)
  
  residuals <- data.frame()
  kappas <- unique(iters$kappa)
  confs <- unique(iters$conf)
  etas <- unique(iters$eta)
  
  require("RColorBrewer")
  clrs <- brewer.pal(n=length(confs),name="Set1")
  pchs <- 15:(15+length(kappas))
  
  for( kappa in kappas ){
    kapiter <- iters[ iters$kappa == kappa, ]
    for( eta in etas ){
      for( conf in confs ){
        singleconf <- kapiter[ ( kapiter$eta == eta & kapiter$conf == conf ), ]
        if(nrow(singleconf)>=1){
          residuals <- rbind( residuals, data.frame(kappa=kappa,eta=eta,conf=conf,
                                                    clr=clrs[which(confs==conf)],pch=pchs[which(kappas==kappa)],
                                                    rmin=min(singleconf$res),rmean=mean(singleconf$res),rmax=max(singleconf$res),
                                                    imin=min(singleconf$iters),imean=mean(singleconf$iters),imax=max(singleconf$iters),
                                                    stringsAsFactors=FALSE)
                            )
        }
      }
    }
  }
  if(debug) print(residuals)
  
  tikzfiles <- tikz.init(basename=datfile,width=5,height=5)
  positions <- c("bottomleft","topleft")
  for(kappa in kappas){
    kres <- residuals[ residuals$kappa == kappa, ]
    plotwitherror(y=kres$rmean,x=kres$eta,
                  dy=kres$rmax-kres$rmean, mdy=kres$rmean-kres$rmin,
                  pch=kres$pch, col=kres$clr, log='y', xlim=range(etas),
                  ylab="$\\|r\\|^2$",xlab="$\\eta$", las=1,
                  main=sprintf("$\\|r\\|^2$ after 500 iterations at $\\kappa=%.3f$",kappa))
    legend(x=positions[which(kappas==kappa)], legend=sprintf("conf.%04d",confs), 
           col=clrs, pch=pchs[which(kappas == kappa)], bty='n')
  }
  tikz.finalize(tikzfiles)
}


