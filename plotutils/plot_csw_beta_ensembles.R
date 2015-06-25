plot_csw_beta_ensembles <- function(datafile, input,cswrange=c(1.00,2.5),betarange=c(1.6,1.85),Nt=500) {

  require("RColorBrewer")

  ensembles <- read.table(file=datafile,header=TRUE,colClasses="numeric")
  ensembles <- ensembles[order(ensembles$csw,ensembles$beta),]

  # one symbol per mass value
  mu <- unique(ensembles$mu)
  syms <- c()
  sym.cur <- 16
  for(i in mu){
    syms <- c(syms,rep(sym.cur, times=length(which(ensembles$mu==i))) )
    sym.cur <- sym.cur + 2
  }
  
  # one colour per mass value with some transparency
  clrs <- c()
  clr.idx <- 1
  clr.pal <- brewer.pal(n=9,name="Set1")
  for(i in mu){
    clr <- col2rgb(clr.pal[clr.idx])
    clr <- rgb(red=clr[1,1],green=clr[2,1],blue=clr[3,1],alpha=166,maxColorValue=255)
    clrs <- c(clrs,rep(clr, times=length(which(ensembles$mu==i))) )
    clr.idx <- clr.idx+1
    # we need to wrap around if there are too many
    if(clr.idx==(length(clr.pal)+1) ) clr.idx <- 1
  }

  tikzfiles <- tikz.init("csw_beta_ensembles",width=3.5,height=3.5,lwdUnit=0.7)
   
  # prepare plot area
  plot(x=ensembles$beta,y=ensembles$csw,
       xlim=betarange, ylim=cswrange,
       ylab="$ c_\\mathrm{sw} $",
       xlab="$ \\beta $", 
       pch=syms, col=clrs, cex=ensembles$Nt/Nt,
       xaxt='n',
       yaxt='n' )
  points(x=ensembles$beta,y=ensembles$csw,pch='.',col='black')
  axis(side=1,at=seq(betarange[1],betarange[2],by=0.025))
  axis(side=2,at=seq(cswrange[1],cswrange[2],by=0.25),las=1)

  lg <- c(sprintf("$ a\\mu = %s $", mu ),sprintf("%d traj.", Nt))
  legend(x="topright",legend=lg,pch=c(unique(syms),1),col=c(unique(clrs),"black"), bty='n',pt.cex=c(rep(1.5,length(unique(ensembles$mu))),1.0))
 
  tikz.finalize(tikzfiles)

}
  
