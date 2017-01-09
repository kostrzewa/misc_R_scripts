source("~/code/R/misc_R_scripts/plotutils/plot_multihist.R")

plot_force_histograms <- function(datfiles,lg,basename="forces_hist_comparison",width=5,height=3) {
  fdat <- list()
  mon <- list()
  fdat[[1]] <- read.table(header=TRUE,file=datfiles[1],stringsAsFactors=FALSE)
  
  mon <- unique(fdat[[1]]$mon)
  
  if(length(datfiles)>1){
    for( i in 2:length(datfiles) ){
      fdat[[length(fdat)+1]] <- read.table(header=TRUE,file=datfiles[i],stringsAsFactors=FALSE)
      tmon <- unique(fdat[[i]]$mon)
      if(any(mon!=tmon))
        stop("plot_force_histograms: monomials inconsistent!!")
    }
  }
  
  x <- 1
  max.df <- NULL
  aver.df <- NULL
  for( m in mon ){
    maxdat <- list()
    averdat <- list()
    dx <- 1
    for( f in fdat ){
      max <- f$max[f$mon==m]
      aver <- f$aver[f$mon==m]

      maxdat[[length(maxdat)+1]] <- max
      averdat[[length(averdat)+1]] <- aver
      
      max.qt <- quantile(max,probs=c(0.0,0.1573,0.5,0.8427,1.0))
      aver.qt <- quantile(aver,probs=c(0.0,0.1573,0.5,0.8427,1.0))
      max.df <- rbind(max.df,data.frame(x=x-0.4+dx*0.2,
                                            val=max.qt[3],
                                            dval1=max.qt[4]-max.qt[3],
                                            dval2=max.qt[5]-max.qt[4],
                                            mdval1=max.qt[3]-max.qt[2],
                                            mdval2=max.qt[2]-max.qt[1],
                                            mon=mon,clr.idx=dx,
                                            stringsAsFactors=FALSE))
      aver.df <- rbind(aver.df,data.frame(x=x-0.4+dx*0.2,
                                            val=aver.qt[3],
                                            dval1=aver.qt[4]-aver.qt[3],
                                            dval2=aver.qt[5]-aver.qt[4],
                                            mdval1=aver.qt[3]-aver.qt[2],
                                            mdval2=aver.qt[2]-aver.qt[1],
                                            mon=mon,clr.idx=dx,
                                            stringsAsFactors=FALSE))
      dx <- dx + 1
    }
    x <- x+1
    if(!missing(lg)){
      try(plot_multihist(dat=maxdat,basename=sprintf("max.%s.%s",basename,m),lg=lg,width=width,height=height,main="",xlab="$\\mathrm{max}(F^2)$"))
      try(plot_multihist(dat=averdat,basename=sprintf("avg.%s.%s",basename,m),lg=lg,width=width,height=height,main="",xlab="$\\mathrm{avg}(F^2)$",factor=10,xlim.probs=c(0,1)))
    } else {
      try(plot_multihist(dat=maxdat,basename=sprintf("max.%s.%s",basename,m),width=width,height=height,main="",xlab="$\\mathrm{max}(F^2)$"))
      try(plot_multihist(dat=averdat,basename=sprintf("avg.%s.%s",basename,m),width=width,height=height,main="",xlab="$\\mathrm{avg}(F^2)$"))
    }
  }
  require("RColorBrewer")
  bpal <- brewer.pal(n=length(fdat),name="Dark2")
  tikzfiles <- tikz.init(basename=sprintf("%s.comparison",basename),width=6,height=6,lwdUnit=0.7)
  plotwitherror(x=max.df$x,y=max.df$val, dy=cbind(max.df$dval1,max.df$dval2), mdy=cbind(max.df$mdval1,max.df$mdval2), 
                col=bpal[max.df$clr.idx],
                log='y',
                yaxt='n',
                ylab="",
                xlab="",xaxt='n',
                errsum.method="linear",
                ylim=c(10^(-9),8*10^2),pch=16)
  plotwitherror(rep=TRUE,
                x=aver.df$x,y=aver.df$val, dy=cbind(aver.df$dval1,aver.df$dval2), mdy=cbind(aver.df$mdval1,aver.df$mdval2), 
                col=bpal[aver.df$clr.idx],
                pch=1,errsum.method="linear")
  legend(legend=lg,fill=bpal,bty='n',x="bottomleft")
  legend(legend=c("$\\mathrm{max}(F^2_i)$","$\\mathrm{avg}(F^2_i)$"),pch=c(16,1),x="topright",bty='n',pt.cex=1.5)
  axis(side=2, labels=TRUE, at=10^((-9):4), tck=0.02,las=1)
  axis(side=2, labels=FALSE, at=outer(2:9,10^((-9):4)), tck=0.01)

  tikz.finalize(tikzfiles)
}


