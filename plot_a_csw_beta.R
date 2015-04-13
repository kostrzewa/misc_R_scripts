plot_a_csw_beta <- function(datafile,cswrange=c(1.0,2.5),betarange=c(1.6,1.8),arange=c(0.05,0.13),kapparange=c(0.125,0.145)) {
  a <- read.table(file=datafile,header=TRUE,colClasses="numeric")
  a <- a[order(a$csw,a$beta),]
  a <- cbind(a,loga=log(a$a))
  a <- cbind(a,dloga=a$da/a$a)
  a <- cbind(a,logkappac=log(a$kappac))
  save(a,file="a_csw_beta_kappac.Rdata")

  # model for log(a)
  w <- 1/a$dloga^2
  amodel <- lm(loga ~ csw + beta, data=a, weights=w )
  save(amodel,file="a_csw_beta.model.Rdata")
  a.cfs <- amodel$coefficients[2:3]
  
  # model for kappa_c
  kappacmodel <- lm(logkappac ~ csw + beta, data=a)
  save(kappacmodel,file="kappac_csw_beta.model.Rdata")
  k.cfs <- kappacmodel$coefficients[2:3]

  csw <- unique(a$csw)
  syms <- c()
  sym <- 0
  for(i in csw){
    syms <- c(syms,rep(sym, times=length(which(a$csw==i))) )
    sym <- sym + 1
  }

  newdata <- data.frame(csw=seq(cswrange[1],cswrange[2],length.out=100),
                        beta=seq(betarange[1],betarange[2],length.out=100))
  a.prediction <- predict(amodel,newdata=newdata,
                        interval="confidence", level=0.68)
  k.prediction <- predict(kappacmodel,newdata=newdata,
                          interval="confidence", level=0.68)


  library(tikzDevice)
  temp <- sprintf("%s.%s","a_kappac_csw_beta",c("tex","pdf","aux","log"))
  tikzfiles <- list(tex=temp[1],pdf=temp[2],aux=temp[3],log=temp[4])
  rm(temp)
  tikz(tikzfiles$tex, standAlone = TRUE, width=3.5, height=3.5)
  
  poly.x <- -(a.cfs[1]*newdata$csw+a.cfs[2]*newdata$beta)
  poly.x <- c(poly.x,rev(poly.x))
  poly.y <- c(a.prediction[,2],rev(a.prediction[,3]))
  poly.col <- rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4)
  xpts <- -(a.cfs[1]*a$csw+a.cfs[2]*a$beta)
  xlims <- -c( a.cfs[1]*cswrange[1]+a.cfs[2]*betarange[1], a.cfs[1]*cswrange[2]+a.cfs[2]*betarange[2] )
  plotwitherror(x=xpts,y=a$a,dy=a$da,log='y',
                xlim=xlims,ylim=arange,
                xlab="$ -(\\alpha c_\\mathrm{sw} + b \\beta) $",
                ylab="$ a (\\mathrm{fm}) $",
                pch=syms )
  polygon(x=poly.x,y=exp(poly.y),col=poly.col)
  csw.legend <- sprintf("$ c_\\mathrm{sw} = %s $", csw)
  legend(x=xlims[1],y=arange[1]+0.3*(arange[2]-arange[1]),legend=csw.legend,pch=unique(syms))
  
  poly.x <- -(k.cfs[1]*newdata$csw+k.cfs[2]*newdata$beta)
  poly.x <- c(poly.x,rev(poly.x))
  poly.y <- c(k.prediction[,2],rev(k.prediction[,3]))
  xpts <- -(k.cfs[1]*a$csw+k.cfs[2]*a$beta)
  xlims <- -c( k.cfs[1]*cswrange[1]+k.cfs[2]*betarange[1], k.cfs[1]*cswrange[2]+k.cfs[2]*betarange[2] )
  plot(x=xpts,y=a$kappac,log='y',
       xlim=xlims,ylim=kapparange,
       xlab="$ -(\\alpha c_\\mathrm{sw} + b \\beta) $",
       ylab="$ \\kappa_c $", 
       pch=syms )
  polygon(x=poly.x,y=exp(poly.y),col=poly.col)
  csw.legend <- sprintf("$ c_\\mathrm{sw} = %s $", csw)
  legend(x=xlims[1],y=kapparange[1]+0.33*(kapparange[2]-kapparange[1]),legend=csw.legend,pch=unique(syms))


  dev.off()                                                                                                                                                                                                  
  tools::texi2dvi(tikzfiles$tex,pdf=T)
  # use pdfcrop tool to remove plot borders
  command <- sprintf("pdfcrop %s %s",tikzfiles$pdf,tikzfiles$pdf)
  system(command)
  # remove temporary files 
  command <- sprintf("rm %s %s %s", tikzfiles$tex, tikzfiles$log, tikzfiles$aux)
  system(command)

}
  
