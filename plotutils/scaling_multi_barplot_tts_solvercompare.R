#scaling_multi_barplot_tts <- function(datfile,basename="scaling.invert",solver='cg',action='tm',dbg=FALSE,show.max=FALSE) {
scaling_multi_barplot_tts <- function(datfile,basename="scaling.invert",solver='cg',dbg=FALSE,show.max=FALSE) {
  require("RColorBrewer")
  
  dat <- read.table(datfile,stringsAsFactors=FALSE,header=TRUE)
  dat <- dat[ dat$nds > 18, ]
  if(dbg) print(dat)

  action_types <- sort(unique(dat$action))
  solvers <- unique(dat$solver)
  pars <- unique(dat$par)
  Ls <- unique(dat$L)
  N_nds <- unique(dat$nds)
 
  for( actions in action_types ){ 
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

    write.table(mindat,file=sprintf("%s.%s.dat",basename,solver),row.names=FALSE,quote=FALSE)

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


    tts <- matrix(nrow=l_Ls*l_pars*l_solvers,ncol=l_N_nds)
    max.tts <- matrix(nrow=l_Ls*l_pars*l_solvers,ncol=l_N_nds)
    for( i in 1:l_N_nds ){
      if(dbg){
        cat(sprintf("On %d nodes, have data\n",N_nds[i] ))
        print( mindat[ mindat$nds==N_nds[i], ] )
        cat("\n")
      }
      tts[ , i ] <- dat[ dat$nds==N_nds[i] & dat$action==actions, ]$tts
      max.tts[ , i ] <- maxdat[ maxdat$nds==N_nds[i] ,]$tts
    }
   
    pal <- brewer.pal(name="Paired",n=8)
    pal <- pal[1:l_solvers]
    tikzFiles <- tikz.init(basename=sprintf("%s.%s",basename,actions),width=4.0,height=4.0)
    ## SPEEDUP
    par(mgp=c(2.5,0.3,0.5))
    par(mar=c(5,4,4,5)+.1)
    mids <- barplot(tts[1,1]/tts,beside=TRUE,las=1,ylab="",xlab="",
                    col=pal,
                    xaxs='i',yaxs='i',
                    yaxp=c(0,ceiling(signif(max(tts[1,1]/tts),2)),ceiling(signif(max(tts[1,1]/tts),2))),
                    tck=0.015)
    
    if(show.max)
      points(x=mids,y=maxdat$tts,pch='-')

    xtckidx <- seq(from=2, to=length(mids), by=l_solvers)
    xtckpos <- mids[xtckidx]-0.45
    if(dbg){
      cat("x tickmark indices and positions\n")
      print(mids[xtckidx])
      print(xtckidx)
      print(xtckpos)
    }

    mtext(side=2,line=3,"Speed-up")
    mtext(side=1,line=0.5,N_nds,at=xtckpos,adj=0.5)
    mtext(side=1,line=1.5,text="Jureca Nodes")

    legend('topleft',cex=0.95,
           legend=sprintf("%s, %s, $L=%d$", c('tmLQCD','QPhiX'), c(rep(actions[1],l_solvers)), rep(48,l_solvers)),
           fill=pal,
           bty='n')
    #############################

    ## Absolute TTS
    par(mgp=c(2.5,0.3,0.5))
    par(mar=c(5,4,4,5)+.1)
    mids <- barplot(tts,beside=TRUE,las=1,ylab="",xlab="",
                    col=pal,
                    ylim=c(0,ceiling(1.2*max(tts))),
                    xaxs='i',yaxs='i',
                    yaxp=c(0,signif(max(tts)/500,1)*500,10),
                    tck=0.015)
    
    if(show.max)
      points(x=mids,y=maxdat$tts,pch='-')

    xtckidx <- seq(from=2, to=length(mids), by=l_solvers)
    xtckpos <- mids[xtckidx]-0.45
    if(dbg){
      cat("x tickmark indices and positions\n")
      print(mids[xtckidx])
      print(xtckidx)
      print(xtckpos)
    }

    mtext(side=2,line=3,"Time to solution (s)")
    mtext(side=1,line=0.5,N_nds,at=xtckpos,adj=0.5)
    mtext(side=1,line=1.5,text="Jureca Nodes")

    legend('topright',cex=0.95,
           legend=sprintf("%s, %s, $L=%d$", c('tmLQCD','QPhiX'),c(rep(actions[1],l_solvers)), rep(48,l_solvers)),
           fill=pal,
           bty='n')
    ###############################
    
    tikz.finalize(tikzFiles)

    print( tts[1,1]/tts )
  } # action types
}
  
