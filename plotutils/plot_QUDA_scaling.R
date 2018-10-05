plot_QUDA_scaling <- function(datfile, basename){
  require("RColorBrewer")
  perf.dat <- read.table(file=datfile, header=TRUE)
  
  Ls <- unique(perf.dat$L)
  recons <- unqiue(perf.dat$recon)
  sloppys <- unique(perf.dat$sloppy)

  clr <- brewer.pal(n=length(recon), name="Dark2") 

  for( L in Ls ){
    Ldat <- perf.dat[ which(perf.dat$L == L), ]
    tikzfiles <- tikz.init(basename=sprintf("%s.L%d",basename,L),width=5,height=5)
    for( meas in c("gflop", "iters", "tts") ){
      for( bind in c(0,1) ){
        dat <- Ldat[ which(Ldat$bind == bind), ]
        
        clrs <- match( dat$recon, recons )
        pchs <- match( dat$sloppy, sloppys )-1
        pchoffset <- 15*dat$p2p
      
        plot(x=dat$np, y=dat[,meas], pch=pchs+pchoffset, col=clrs,
             xlab="no. GPUs", ylab=meas )  
      }
    }
    tikz.finalize(tikzfiles)
  }
}
