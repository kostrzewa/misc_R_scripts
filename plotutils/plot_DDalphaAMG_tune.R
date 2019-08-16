require(hadron)

plot_DDalphaAMG_tune <- function(datfile, basename, citer.ylim=c(2000,15000) ){ 
  require("RColorBrewer")

  perf.dat <- read.table(file=datfile, header=TRUE)

  perf.dat <- perf.dat[order(perf.dat$solvetime),]
  print(perf.dat)

  lvls <- unique(perf.dat$lvl)
  iter <- unique(perf.dat$setupiter)
  mucoarse <- unique(perf.dat$mucoarse)
  nvecs <- unique(perf.dat$nvec)

  max.perf <- NULL

  tikzfiles <- tikz.init(basename=basename,width=5,height=5)
  for( lvl in lvls ){
    dat <- perf.dat[ which(perf.dat$lvl == lvl), ]
    nvecs <- unique(dat$nvec)
    clr <- brewer.pal(name="Set2", n=length(nvecs))

    plot(x=dat$mucoarse, y=dat$setuptime, pch=dat$setupiter, col=clr[ match(x=dat$nvec, table=nvecs) ], 
         main=sprintf("$N_{\\mathrm{lvl}}=%d$",lvl), xlab="$\\mu_c$ factor", ylab="Setup time [s]",las=1,
         xaxt='n', log='y')
    axis(side=1,at=unique(dat$mucoarse))
    legend(x="topright",legend=c( sprintf("$N_{\\mathrm{setup}}^{\\mathrm{iter}}=%d$",iter),
                                  sprintf("$N_{\\mathrm{vec}}=%d$",nvecs) ),
           pch=c(iter,rep(15,length(nvecs))), col=c(rep("black",length(iter)),clr), bty='n')
    
    plot(x=dat$mucoarse, y=dat$solvetime, pch=dat$setupiter, col=clr[ match(x=dat$nvec, table=nvecs) ],  
         main=sprintf("$N_{\\mathrm{lvl}}=%d$",lvl), xlab="$\\mu_c$ factor", ylab="Solve time [s]",las=1,
         xaxt='n', log='y')
    axis(side=1,at=unique(dat$mucoarse))
    legend(x="topright",legend=c( sprintf("$N_{\\mathrm{setup}}^{\\mathrm{iter}}=%d$",iter),
                                  sprintf("$N_{\\mathrm{vec}}=%d$",nvecs) ),
           pch=c(iter,rep(15,length(nvecs))), col=c(rep("black",length(iter)),clr), bty='n')

    plot(x=dat$mucoarse, y=dat$solvefiter, pch=dat$setupiter, col=clr[ match(x=dat$nvec, table=nvecs) ],  
         main=sprintf("$N_{\\mathrm{lvl}}=%d$",lvl), xlab="$\\mu_c$ factor", ylab="Fine grid iterations",las=1,
         xaxt='n')
    axis(side=1,at=unique(dat$mucoarse))
    legend(x="topleft",legend=c( sprintf("$N_{\\mathrm{setup}}^{\\mathrm{iter}}=%d$",iter),
                                  sprintf("$N_{\\mathrm{vec}}=%d$",nvecs) ),
           pch=c(iter,rep(15,length(nvecs))), col=c(rep("black",length(iter)),clr), bty='n')
    
    plot(x=dat$mucoarse, y=dat$solveciter, pch=dat$setupiter, col=clr[ match(x=dat$nvec, table=nvecs) ],  
         main=sprintf("$N_{\\mathrm{lvl}}=%d$",lvl), xlab="$\\mu_c$ factor", ylab="Total coarse grid iterations",las=1,
         xaxt='n', ylim=citer.ylim)
    axis(side=1,at=unique(dat$mucoarse))
    legend(x="topright",legend=c( sprintf("$N_{\\mathrm{setup}}^{\\mathrm{iter}}=%d$",iter),
                                  sprintf("$N_{\\mathrm{vec}}=%d$",nvecs) ),
           pch=c(iter,rep(15,length(nvecs))), col=c(rep("black",length(iter)),clr), bty='n')

    plot(x=dat$mucoarse, y=dat$solvecgrid, pch=dat$setupiter, col=clr[ match(x=dat$nvec, table=nvecs) ],  
         main=sprintf("$N_{\\mathrm{lvl}}=%d$",lvl), xlab="$\\mu_c$ factor", ylab="Coarse grid \\%-age",las=1,
         xaxt='n')
    axis(side=1,at=unique(dat$mucoarse))
    legend(x="bottomleft",legend=c( sprintf("$N_{\\mathrm{setup}}^{\\mathrm{iter}}=%d$",iter),
                                  sprintf("$N_{\\mathrm{vec}}=%d$",nvecs) ),
           pch=c(iter,rep(15,length(nvecs))), col=c(rep("black",length(iter)),clr), bty='n')


  }
  tikz.finalize(tikzfiles)
}

