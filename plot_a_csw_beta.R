plot_a_csw_beta <- function(datafile,cswrange=c(1.0,2.5),betarange=c(1.6,1.8),arange=c(0.05,0.13)) {
  a <- read.table(file=datafile,header=TRUE,colClasses="numeric")
  a <- a[order(a$csw,a$beta),]
  a <- cbind(a,loga=log(a$a))
  a <- cbind(a,dloga=a$da/a$a)
  w <- 1/a$dloga^2
  amodel <- lm(loga ~ csw + beta, data=a, weights=w )
  cfs <- amodel$coefficients[2:3]
  
  csw <- unique(a$csw)
  syms <- c()
  sym <- 0
  for(i in csw){
    syms <- c(syms,rep(sym, times=length(which(a$csw==i))) )
    sym <- sym + 1
  }

  xpts <- -(cfs[1]*a$csw+cfs[2]*a$beta)
  xlims <- -c( cfs[1]*cswrange[1]+cfs[2]*betarange[1], cfs[1]*cswrange[2]+cfs[2]*betarange[2] )
  
  newdata <- data.frame(csw=seq(cswrange[1],cswrange[2],length.out=100),
                        beta=seq(betarange[1],betarange[2],length.out=100))
  prediction <- predict(amodel,newdata=newdata,
                        interval="confidence", level=0.68)

  poly.x <- -(cfs[1]*newdata$csw+cfs[2]*newdata$beta)
  poly.x <- c(poly.x,rev(poly.x))
  poly.y <- c(prediction[,2],rev(prediction[,3]))
  poly.col <- rgb(red=1.0,green=0.0,blue=0.0,alpha=0.4)
 
  library(tikzDevice)
  temp <- sprintf("%s.%s","a_csw_beta",c("tex","pdf","aux","log"))
  tikzfiles <- list(tex=temp[1],pdf=temp[2],aux=temp[3],log=temp[4])
  rm(temp)

  tikz(tikzfiles$tex, standAlone = TRUE, width=5, height=5)
  plotwitherror(x=xpts,y=a$a,dy=a$da,log='y',
                xlim=xlims,ylim=arange,
                xlab="$ \\alpha c_\\mathrm{sw} + b \\beta $",
                ylab="$ a (\\mathrm{fm}) $",
                pch=syms )
  polygon(x=poly.x,y=exp(poly.y),col=poly.col)


  csw.legend <- sprintf("$ c_\\mathrm{sw} = %s $", csw)
  legend(x=xlims[1],y=arange[1]+0.2*(arange[2]-arange[1]),legend=csw.legend,pch=unique(syms))
  dev.off()                                                                                                                                                                                                  
  tools::texi2dvi(tikzfiles$tex,pdf=T)
  # use pdfcrop tool to remove plot borders
  command <- sprintf("pdfcrop %s %s",tikzfiles$pdf,tikzfiles$pdf)
  system(command)
  # remove temporary files 
  command <- sprintf("rm %s %s %s", tikzfiles$tex, tikzfiles$log, tikzfiles$aux)
  system(command)

}
  
