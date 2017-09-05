plot_peram_iters <- function(filename="iters.dat"){
  iters <- read.table(filename, header=TRUE); 
  
  conf <- unique(iters$conf)
  rnd <- unique(iters$rnd)
  
  require("RColorBrewer")
  
  pal <- brewer.pal(n=length(rnd), name="Dark2") 
  
  tikzfiles <- tikz.init(basename="iters");
  plot(y=iters$iters,
       x=iters$conf,
       type='n',
       xlab="$n_\\mathrm{conf}$",
       ylab="$N_\\mathrm{iters}$")
  for( cnf in conf ){
    for( rn in rnd ){
      i <- cnf + rn/length(rnd)
      liters <- iters[which(iters$conf == cnf & iters$rnd == rn), ]
      pt <- mean(liters$iters)
      pt.dy <- max(liters$iters)-pt
      pt.mdy <- pt-min(liters$iters)
      plotwitherror(x=i, y=pt, dy=pt.dy, mdy=pt.mdy, pch='.', cex=3, rep=TRUE, col=pal[rn+1])
    }
  }
  legend(x="topright", col=pal, legend=sprintf("$n_r=%d$",rnd), pt.cex=5, pch='.')
  tikz.finalize(tikzfiles)
}
