scaling_inverse_tts_solvercompare <- function(datfile,basename="scaling.invert",solver='cg',dbg=FALSE,show.max=FALSE) {
  require("RColorBrewer")
  allclrs <- brewer.pal(n=8, name="Dark2")
  allpchs <- c(17,8)

  dat <- read.table(datfile,stringsAsFactors=FALSE,header=TRUE)
  dat$one_ov_tts <- 1.0/dat$tts
  if(dbg) print(dat)

  action_types <- unique(dat$action)
  solvers <- unique(dat$solver)
  pars <- unique(dat$par)
  Ls <- unique(dat$L)
  N_nds <- unique(dat$nds)

  tikzfiles <- tikz.init(basename=basename,width=4.5,height=4.5,lwdUnit=1.0) 
  for( action in action_types ){
    ldat <- dat[ which(dat$action == action), ]
    plot(t='n',
         x=ldat$nds,
         y=ldat$one_ov_tts,
         xlab="$N_\\mathrm{nds}$",
         ylab="$1/t_s$ $\\mathrm{[sec]}^{-1}$")
    for( par in pars ){
      for( solver in solvers ){
        sidx <- which(ldat$solver == solver & ldat$par == par)
        points(x=ldat$nds[sidx],
               y=ldat$one_ov_tts[sidx],
               pch=allpchs[ which(solvers == solver) ],
               col=allclrs[ which(pars == par) ])
      }
      legend(x="topleft",
             bty='n',
             pch=rep(15,2),
             col=allclrs[1:2],
             legend=pars)

      legend(x="bottomright",
             bty='n',
             pch=rev(allpchs[1:2]),
             col=rep("black",2),
             legend=rev(solvers))
    }
  }
  tikz.finalize(tikzfiles)
}
  
