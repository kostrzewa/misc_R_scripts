scaling_multi_barplot <- function(datfile,basename="scaling",dbg=FALSE,show.mins=FALSE,show.coms=FALSE,yscale=1e6,flopunit="Tflop/s") {
  require("RColorBrewer")
  
  dat <- read.table(datfile,stringsAsFactors=FALSE,header=TRUE)

  pars <- unique(dat$par)
  Ls <- unique(dat$L)
  N_nds <- unique(dat$nds)

  comms <- c(0,1)
  cmstext <- c("nocomms","comms")
 
  for( c_idx in 1:length(comms) ){
    # for all parallelisation combinations for each lattice extent, parallelisation type and number of nodes, extract the maximum
    # performance with communication
    maxdat <- NULL
    cmsdat <- NULL
    mindat <- NULL
    for(L in Ls){
      for(i_nds in 1:length(N_nds)){
        for(i_pars in 1:length(pars)){
          clsidx <- which( dat$nds==N_nds[i_nds] & dat$par==pars[i_pars] & dat$comm==comms[c_idx] & dat$L==L )
          clsdat <- dat[clsidx,]
          maxdat <- rbind(maxdat,data.frame( clsdat[ which( clsdat$mflops==max(clsdat$mflops) )[1],] ) ) 
          mindat <- rbind(mindat,data.frame( clsdat[ which( clsdat$mflops==min(clsdat$mflops) )[1],] ) )
          
          com.clsidx <- which( dat$nds==N_nds[i_nds] & dat$par==pars[i_pars] & dat$comm==xor(comms[c_idx],1) & dat$L==L )
          com.clsdat <- dat[com.clsidx,]
          cmsdat <- rbind(cmsdat,data.frame( com.clsdat[ which( com.clsdat$mflops==max(com.clsdat$mflops) )[1],] ) ) 

        }
      }
    }
    #save(maxdat,file=sprintf("max.%s.Rdata",cmstext[c_idx]))
    #save(mindat,file=sprintf("min.%s.Rdata",cmstext[c_idx]))
    
    maxdat <- maxdat[ order(maxdat$nds,maxdat$L,maxdat$par), ]
    mindat <- mindat[ order(mindat$nds,mindat$L,mindat$par), ]
    cmsdat <- cmsdat[ order(cmsdat$nds,cmsdat$L,cmsdat$par), ]

    write.table(maxdat,row.names=FALSE,file=sprintf("%s.%s.dat",basename,cmstext[c_idx]))

    maxdat$mflops <- maxdat$mflops/yscale
    mindat$mflops <- mindat$mflops/yscale
    cmsdat$mflops <- cmsdat$mflops/yscale

    if(dbg){
      cat("\nMaximum performance in each class\n")
      print(maxdat); cat("\n")
      cat("\nMinimum performance in each class\n")
      print(mindat); cat("\n")
    }

    l_Ls <- length(Ls)
    l_N_nds <- length(N_nds)
    l_pars <- length(pars)


    mflops <- matrix(nrow=l_Ls*l_pars,ncol=l_N_nds)
    full.mflops <- matrix(nrow=nrow(dat)/l_N_nds,ncol=l_N_nds)
    for( i in 1:l_N_nds ){
      if(dbg){
        cat(sprintf("On %d nodes, have data\n",N_nds[i] ))
        print( maxdat[ maxdat$nds==N_nds[i], ] )
        cat("\n")
      }
      mflops[, i ] <- maxdat[ maxdat$nds==N_nds[i] ,]$mflops
      full.mflops[, i ] <- dat[ dat$nds==N_nds[i] ,]$mflops
    }
    
    pal <- c(brewer.pal(n=l_pars,name="Blues"),brewer.pal(n=l_pars,name="Reds"))


    tikzFiles <- tikz.init(basename=sprintf("%s.%s",basename,cmstext[c_idx]),width=4.0,height=4.0)
    par(mgp=c(2.5,0.3,0.5))
    par(mar=c(5,4,4,5)+.1)
    mids <- barplot(mflops,beside=TRUE,las=1,ylab="",xlab="",
                    col=pal,
                    ylim=c(0,1.3*max(maxdat$mflops)),xaxs='i',yaxs='i',
                    tck=0.015)
    
    if(show.mins)
      points(x=mids,y=mindat$mflops,pch='-')
    
    if(show.coms)
      points(x=mids,y=cmsdat$mflops,pch='-')

    xtckidx <- seq(from=4, to=length(mids), by=l_Ls*l_pars)
    xtckpos <- mids[xtckidx]-(mids[xtckidx]-mids[xtckidx-1])/2
    if(dbg){
      cat("x tickmark indices and positions\n")
      print(xtckidx)
      print(xtckpos)
    }

    mtext(side=2,line=2,flopunit)
    mtext(side=1,line=0.5,N_nds,at=xtckpos)
    mtext(side=1,line=1.5,text="Jureca Nodes")

    legend('topleft',cex=0.77,
           legend=sprintf("$L=%d$, %s",c(rep(32,3),rep(48,3)), rep(c("hybrid","hybrid (overlap comms)","MPI"),2)),
           fill=pal,
           bty='n')

    # fit straight lines to each parallelisation and select the one with the maximum slope
    lns <- list()
    for( i_L in 1:l_Ls ){
      lslop <- NULL
      t_lns <- list()
      for( i_pars in 1:l_pars ){
        clsidx <- which( maxdat$par==pars[i_pars] & maxdat$L==Ls[i_L] )
        t_lns[[length(t_lns)+1]] <- lm(mflops~x, data=cbind(maxdat[clsidx,], x=mids[clsidx]), subset=which(maxdat[clsidx,]$mflops > 100/yscale) )
        lslop <- c(lslop,summary(t_lns[[length(t_lns)]])$coefficients[2])
      }
      print(lslop)
      maxidx <- which( lslop==max(lslop) )
      lns[[length(lns)+1]] <- t_lns[[ maxidx ]]
    }

    abline(lns[[1]],col="blue")
    abline(lns[[2]],col="red")

    tikz.finalize(tikzFiles)

##########################################################

    tikzFiles <- tikz.init(basename=sprintf("full.%s.%s",basename,cmstext[c_idx]),width=8.0,height=4.0)
    mids <- barplot(full.mflops/yscale,beside=TRUE,las=1,ylab="",xlab="",
                    #col=pal,
                    ylim=c(0,1.3*max(dat$mflops)/yscale),xaxs='i',yaxs='i',
                    tck=0.015)
    
    xtckidx <- seq(from=nrow(dat)/l_N_nds/2, to=length(mids), by=nrow(dat)/l_N_nds)
    xtckpos <- mids[xtckidx]-(mids[xtckidx]-mids[xtckidx-1])/2
    if(dbg){
      cat("x tickmark indices and positions\n")
      print(xtckidx)
      print(xtckpos)
    }

    mtext(side=2,line=2,text=flopunit)
    mtext(side=1,line=0.5,N_nds,at=xtckpos)
    mtext(side=1,line=1.5,text="Jureca Nodes")
    tikz.finalize(tikzFiles)
###########################################################
  }

}
  
