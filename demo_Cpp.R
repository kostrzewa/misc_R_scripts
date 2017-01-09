# function fits to 

# to get a fit in physical units, simply replace with the physical hbarc 

#hbarc <- 197.326972
#dhbarc <-  0.000005

hbarc <- 1
dhbarc <- 0

a <- 1
da <- 0

#a <- 0.091
#da <- 0.006

# m_ps_phys <- 130.4 

m_ps_phys <- 0

demo_Cpp <- function(prefix,T,i_start,i_end,addon,t1,t2,single,plateau.t1=t1,plateau.t2=t2,debug=FALSE,sym=FALSE)
{
  Tover2 <- floor(T/2)

  Cpp.raw <- NULL
  Cpp.lengths <- NULL
  Cpp.missing <- NULL

  # read the data files
  for(i in i_start:i_end){
    filename <- sprintf("%s/Cpp.data.%06d%s",prefix,i,addon)
    readout <- try(read.table(filename)$V2,silent=TRUE)
    # skip filenames which can't be found by misusing the class attribute
    if( length(attr(readout,"class")) != 0) {
      Cpp.missing <- c(Cpp.missing,i)
      next
    } else {
      Cpp.data <- readout
    }
    Cpp.raw <- rbind(Cpp.raw,Cpp.data)
  }

  if(length(Cpp.missing)>0) {
    print(sprintf("There were %d missing files!",length(Cpp.missing)))
  }

  nomeas <- length(Cpp.raw[,1])
  Cpp.mean <- NULL
  Cpp.error <- NULL
  m_eff <- NULL
  # m_eff.error <- NULL
  # m_eff.estimate <- NULL # keep track of bootstrap estimate of m_eff
  m_eff.estimate.uwerr <- NULL
  m_eff.error.uwerr <- NULL # do gamma method error analysis for comparison
  for(t in 1:(Tover2+1)){
    meanval <- mean(Cpp.raw[,t])
    errorval <- sqrt(var(Cpp.raw[,t])/(nomeas))
    Cpp.mean <- c(Cpp.mean,meanval)
    Cpp.error <- c(Cpp.error,errorval)
  }
  
  # function to solve for effective mass in symmetric or non-symmetric manner
  # the 'value' input is for finding the root of f(E,t) - value(t) = 0
  corratio_fun <- function(E,value=0,t,N_t) {
    Thalf <- N_t/2 
    if(sym){
      return( ( cosh((Thalf-t-1)*E) - cosh((Thalf-t+1)*E) ) / ( cosh((Thalf-t)*E) ) - value )
    } else {
      return( (cosh((Thalf-t)*E) / cosh((Thalf-t-1)*E)) - value )
    }
  }
  
  # remember that array index [1] is equivalent to t=0!
  corratio_val <- function(t) {
    if(sym) {
      return( ( Cpp.mean[t+2] - Cpp.mean[t] ) / Cpp.mean[t+1] )
    } else {
      return( Cpp.mean[t+1]/Cpp.mean[t+2] )
    }
  }

  corratio_err <- function(t) {
    if(sym) {
      return( (1/Cpp.mean[t+1])*(Cpp.error[t+2] - Cpp.error[t] - Cpp.error[t+1]*corratio_val(t)))
    } else {
      return( (1/Cpp.mean[t+2])*(Cpp.error[t+1]-Cpp.error[t+2]/Cpp.mean[t+2]) )
    }
  } 

  uniroot_boot <- function(data,t_val) {
    if(sym) {
      val <- (mean(data[,t_val+2]) - mean(data[,t_val]))/ mean( data[,t_val+1] )
    } else {
      val <- (mean(data[,t_val+1])/mean(data[,t_val+2]))
    }
    solution <- uniroot(corratio_fun,value=val,t=t_val,N_t=T,lower=0,upper=3)
    return(solution$root)
  }

  # data should be in the form data[1]=C(t+1),data[2]=C(t-1),data[3]=C(t)
  uniroot_uwerr <- function(data,t_val) {
    val <- NULL
    if(sym) {
      val <- (data[1]-data[2])/data[3]
    } else {
      val <- (data[3]/data[1])
    }
    solution <- uniroot(corratio_fun,value=val,t=t_val,N_t=T,lower=0,upper=3)
    return(solution$root)
  }

  noop_tsboot <- function(data) {
    return(data)
  }

  noop_boot <- function(original,indices) {
    return(original[indices])
  }

  # a slightly more elaborate fit which fits against a functional form
  picorrel <- function(t,a,E,N_t) {
    return( a*( exp(-t*E) + exp(-(N_t-t)*E) ) )
  }

  # from the range used for the effective mass estimate, derive a range in which
  # the correlator fit will be attempted directly
  picorfit.t1 <- plateau.t1
  picorfit.t2 <- plateau.t2
  # generate a data frame for this range
  Cppdata <- data.frame(t=seq(picorfit.t1,picorfit.t2,1),
    Cpp=Cpp.mean[(picorfit.t1+1):(picorfit.t2+1)],
    weight=(1/Cpp.error^2)[(picorfit.t1+1):(picorfit.t2+1)])
  
  # plot raw correlator data points 
  cppmax <- max(Cpp.mean)
  cppmin <- min(Cpp.mean)
  # sometimes the correlator might become negative, here we catch that case
  # for the log plot
  if( cppmin < 0 ){
    cppmin <- 0.00000001
  }

  pp_filename <- sprintf("pp_correl_%s.pdf",prefix);

  pdf(pp_filename,title=prefix,width=5,height=5)
  op <- par(family="Palatino")
  plotwitherror(x=seq(0,Tover2,1),
    ylim=c(cppmin,8),
    log="y"
    ,y=Cpp.mean,
    dy=20*Cpp.error,
    xlab=NA,ylab=expression(log(C[pp](tau))),
    xaxt="n",yaxt="n",
    pch=NA,col="forestgreen",
    lwd=3,mgp=c(0,0,0))
    
#  abline(v=picorfit.t1,col='gray')
#  abline(v=picorfit.t2,col='gray')
  
  picorfit <- nls(Cpp~picorrel(t,a,E,T),data=Cppdata,
    weights=Cppdata$weight,
    start=list(a=1,E=0.06),
    model=TRUE)

  print(summary(picorfit))

  picorfit$a <- coef(summary(picorfit))["a","Estimate"]
  picorfit$da <- coef(summary(picorfit))["a","Std. Error"]
  picorfit$dE <- coef(summary(picorfit))["E","Std. Error"]
  picorfit$E <- coef(summary(picorfit))["E","Estimate"]
  
  cppmax <- max(Cpp.mean)
  cppmin <- min(Cpp.mean)
  if( cppmin < 0 ){
    cppmin <- 0.00001
  }

  plot(xlim=c(0,Tover2),function(x) picorrel(t=x,a=picorfit$a,E=picorfit$E,N_t=T),add=TRUE,col="violetred",lwd=2)
#  plot(xlim=c(0,Tover2),function(x) picorrel(t=x,a=picorfit$a-picorfit$da,E=picorfit$E-picorfit$dE,N_t=T),add=TRUE,col="pink")
#  plot(xlim=c(0,Tover2),function(x) picorrel(t=x,a=picorfit$a+picorfit$da,E=picorfit$E+picorfit$dE,N_t=T),add=TRUE,col="pink")
  if(!missing(single)){
    points(xlim=c(0,Tover2),Cpp.raw[single,],col="black",pch="*",cex=1.2)
  }

  legendlabels <- c("on single gauge conf.","ensemble average",expression(A(exp(-m*tau)+exp(-m(T-tau)))) )
  legend(x=9,y=8,pch=c("*","I",NA),lty=c(0,0,1),col=c("black","forestgreen","violetred"),legend=legendlabels,bty="n",lwd=c(0,0,2))

  label <- sprintf(" = %.5f (%.5f)",picorfit$E,picorfit$dE)
  #text(x=T*2/3,y=cppmax*2/3,substitute(paste(am[eff],label)))

  dev.off() # pp_filename
}
