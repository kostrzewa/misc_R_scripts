#scaling_multi_barplot_tts <- function(datfile,basename="scaling.invert",solver='cg',action='tm',dbg=FALSE,show.max=FALSE) {
scaling_multi_barplot_tts <- function(datfile,basename="scaling.invert",solver='cg',dbg=FALSE,show.max=FALSE) {
  require("RColorBrewer")
  
  dat <- read.table(datfile,stringsAsFactors=FALSE,header=TRUE)

  actions <- sort(unique(dat$action))
  solvers <- unique(dat$solver)
  pars <- unique(dat$par)
  Ls <- unique(dat$L)
  N_nds <- unique(dat$nds)
  
  # for all parallelisation combinations for each lattice extent, parallelisation type and number of nodes, extract the maximum
  # performance with communication
  maxdat <- NULL
  mindat <- NULL
  for(L in Ls){
    for(action in actions){
      for(i_nds in 1:length(N_nds)){
        for(i_pars in 1:length(pars)){
          clsidx <- which( dat$nds==N_nds[i_nds] & dat$par==pars[i_pars] & dat$L==L & dat$solver==solver & dat$action==action )
          clsdat <- dat[clsidx,]
          maxdat <- rbind(maxdat,data.frame( clsdat[ which( clsdat$tts==max(clsdat$tts) )[1],] ) ) 
          mindat <- rbind(mindat,data.frame( clsdat[ which( clsdat$tts==min(clsdat$tts) )[1],] ) ) 
        }
      }
    }
  }

  mindat <- mindat[ order(mindat$nds,mindat$L,mindat$action,mindat$par), ]
  maxdat <- maxdat[ order(maxdat$nds,maxdat$L,mindat$action,maxdat$par), ]

  write.table(mindat,file=sprintf("%s.%s.dat",basename,solver),row.names=FALSE)

  if(dbg){
    cat("\nMaximum tts in each class\n")
    print(maxdat); cat("\n")
    cat("\nMinimum tts in each class\n")
    print(mindat); cat("\n")
  }

  l_Ls <- length(Ls)
  l_actions <- length(actions)
  l_solvers <- length(solvers)
  l_N_nds <- length(N_nds)
  l_pars <- length(pars)


  tts <- matrix(nrow=l_Ls*l_pars*l_actions,ncol=l_N_nds)
  max.tts <- matrix(nrow=l_Ls*l_pars*l_actions,ncol=l_N_nds)
  for( i in 1:l_N_nds ){
    if(dbg){
      cat(sprintf("On %d nodes, have data\n",N_nds[i] ))
      print( mindat[ mindat$nds==N_nds[i], ] )
      cat("\n")
    }
    tts[, i ] <- mindat[ mindat$nds==N_nds[i] ,]$tts
    max.tts[, i ] <- maxdat[ maxdat$nds==N_nds[i] ,]$tts
  }
  
  pal <- c(brewer.pal(n=l_pars,name="Blues"),brewer.pal(n=l_pars,"Reds"))
  tikzFiles <- tikz.init(basename=sprintf("%s.%s",basename,solver),width=4.0,height=4.0)
  par(mgp=c(2.5,0.3,0.5))
  par(mar=c(5,4,4,5)+.1)
  mids <- barplot(tts,beside=TRUE,las=1,ylab="",xlab="",
                  col=pal,
                  ylim=c(0,ceiling(1.2*max(tts))),xaxs='i',yaxs='i',
                  tck=0.015)
  
  if(show.max)
    points(x=mids,y=maxdat$tts,pch='-')

  #xtckidx <- seq(from=2, to=length(mids), by=l_Ls*l_pars)
  xtckidx <- seq(from=4, to=length(mids), by=l_Ls*l_pars*l_actions)
  xtckpos <- mids[xtckidx]-(mids[xtckidx]-mids[xtckidx-1])/2
  if(dbg){
    cat("x tickmark indices and positions\n")
    print(mids[xtckidx])
    print(xtckidx)
    print(xtckpos)
  }

  mtext(side=2,line=3,"time to solution/s")
  mtext(side=1,line=0.5,N_nds,at=xtckpos,adj=0.5)
  mtext(side=1,line=1.5,text="Jureca Nodes")

  legend('topright',cex=0.75,
         legend=sprintf("$L=%d$, %s, %s",rep(48,6), c(rep(actions[1],3),rep(actions[2],3)), rep(c("hybrid","hybrid (overlap comms)","MPI"),2) ),
         fill=pal,
         bty='n')

  xlims <- par("usr")[1:2]
  
  tikz.finalize(tikzFiles)
}
  
