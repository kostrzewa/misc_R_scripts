source("~/code/R/misc_R_scripts/lat_phys_ratios/compute_ratio.R")

analysis_gradient_flow <- function(path,read.data=TRUE,plot=FALSE,skip=0) {
  if(read.data) {
    raw.gradflow <- readgradflow(path=path,skip=skip)
    save(raw.gradflow,file="raw.gradflow.Rdata",compress=FALSE)
  }else{
    cat("Warning, reading data from raw.gradflow.Rdata, if the number of samples changed, set read.data=TRUE to reread all output files\n")
    load("raw.gradflow.Rdata")
  }

  t_vec <- unique(raw.gradflow$t)
  Ncol <- ncol(raw.gradflow)
  Nrow <- length(t_vec)
  # allocate some memory, for each observable, have space for the value, the error and the autocorrelation time
  # create a list with NULL rownams and sensible column names for this purpose
  subnames <- c("value","dvalue","ddvalue","tauint","dtauint")
  cnames <- colnames(raw.gradflow[,3:Ncol])
  outer.cnames <- t(outer(cnames,subnames,FUN=function(x,y) { paste(x,y,sep=".") }))
  grad.dimnames <- list()
  grad.dimnames[[1]] <- NULL
  grad.dimnames[[2]] <- c("t",as.vector(outer.cnames) )
  
  gradflow <- as.data.frame(matrix(data=NA,nrow=Nrow,ncol=(length(subnames)*(Ncol-2)+1),dimnames=grad.dimnames))
  for(i_row in 1:length(t_vec)){
    uwerr.gradflow <- apply(X=raw.gradflow[which(raw.gradflow$t==t_vec[i_row]),3:Ncol],MARGIN=2,FUN=uwerrprimary)
    summaryvec <- c(t_vec[i_row])
    for(i_col in 1:length(cnames)) {
      obs <- uwerr.gradflow[[cnames[i_col]]]
      summaryvec <- c(summaryvec, obs$value, obs$dvalue, obs$ddvalue, obs$tauint, obs$dtauint)
    }
    gradflow[i_row,] <- summaryvec
  }
  
  save(gradflow,file="gradflow.Rdata",compress=FALSE)
   
  # find w_0 and its lower and upper values
  w0sq <- c( approx(x=gradflow$Wsym.value+gradflow$Wsym.dvalue,y=gradflow$t,xout=0.3)$y, 
             approx( x=gradflow$Wsym.value, y=gradflow$t, xout=0.3 )$y, 
             approx( x=gradflow$Wsym.value-gradflow$Wsym.dvalue, y=gradflow$t, xout=0.3)$y )
  
  cat(sprintf("w0^2/a^2: %f (+%f -%f)\n",w0sq[2],w0sq[3]-w0sq[2],w0sq[2]-w0sq[1]))
   
  w0 <- sqrt(w0sq)
  cat(sprintf("w0/a: %f (+%f -%f)\n",w0[2],w0[3]-w0[2],w0[2]-w0[1]))
   

  w0fm <- list(val=0.1755,dval=0.0019,name="w_0(fm)")
  a <- matrix(nrow=3,ncol=2)
  for(i in 1:length(w0)){
    a_fm <- compute_ratio(dividend=w0fm,divisor=list(val=w0[i],dval=0,name="w_0"),name="a(fm)")
    a[i,] <- c(a_fm$val,a_fm$dval)
  }
  cat(sprintf("Using w0 = %f (%f)\n",w0fm$val,w0fm$dval))
  cat("a(fm):", a[2,1], " (", a[3,1]-a[2,1], ",", a[1,1]-a[2,1]," ) ( ", a[2,2], ")\n")
  cat("a(fm) ~ ", a[2,1], "(", sqrt( 0.5*(abs(a[3,1]-a[2,1])+abs(a[1,1]-a[2,1]))^2 + a[2,2]^2 ), ")\n");   
  
  if(plot) {
    # set up plot
    plot(x=gradflow$t, y=gradflow$Wsym.value,
         type='n',xlim=c(0,4),ylim=c(0.0,0.4),
         xlab="t/a^2",ylab="W(t)")
    # draw errorband
    poly.col <- rgb(red=1.0,green=0.0,blue=0.0,alpha=0.6)
    poly.x <- c(gradflow$t,rev(gradflow$t))
    poly.y <- c(gradflow$Wsym.value+gradflow$Wsym.dvalue,rev(gradflow$Wsym.value-gradflow$Wsym.dvalue))
    polygon(x=poly.x,y=poly.y,col=poly.col)
    abline(h=0.3)
    abline(v=w0sq)
  }
  
}

