eval_Cppdata <- function(prefix,T,i_start,i_end,addon,t1,t2,plateau.t1=t1,plateau.t2=t2,debug=FALSE)
{
  # produce one plot with two columns
  par(mfrow=c(1,2))
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
    # in the case where only half the time extent has been measured (online f.ex.)
    # mirror the correlator  
    if( length(Cpp.data) < T ) {
      for(t in 0:(Tover2-1)) {
        Cpp.data <- c(Cpp.data,Cpp.data[Tover2-t])
      }
    }
    Cpp.raw <- rbind(Cpp.raw,Cpp.data)
  }

  if(length(Cpp.missing)>0) {
    print("There were missing files!")
  }

  Cpp.mean <- NULL
  Cpp.error <- NULL
  m_eff <- NULL
  m_eff.error <- NULL
  m_eff.estimate <- NULL # keep track of bootstrap estimate of m_eff
  m_eff.estimate.uwerr <- NULL
  m_eff.error.uwerr <- NULL # do gamma method error analysis for comparison
  for(t in 1:T){
    meanval <- mean(Cpp.raw[,t])
    errorval <- sqrt(var(Cpp.raw[,t])/(i_end-i_start))
    Cpp.mean <- c(Cpp.mean,meanval)
    Cpp.error <- c(Cpp.error,errorval)
  }
  
  # function to solve for effective mass in symmetric manner
  # the 'value' input is for finding the root of f(E,t) - value(t) = 0
  corratio_fun <- function(E,value=0,t,N_t,sym=TRUE) {
    Thalf <- N_t/2 
    if(sym){
      return( ( cosh((Thalf-t-1)*E) - cosh((Thalf-t+1)*E) ) / ( cosh((Thalf-t)*E) ) - value )
    } else {
      return( (cosh((Thalf-t)*E) / cosh((Thalf-t-1)*E)) - value )
    }
  }
  
  # remember that array index [1] is equivalent to t=0!
  corratio_val <- function(t,sym=TRUE) {
    if(sym) {
      return( ( Cpp.mean[t+2] - Cpp.mean[t] ) / Cpp.mean[t+1] )
    } else {
      return( Cpp.mean[t+1]/Cpp.mean[t+2] )
    }
  }

  corratio_err <- function(t,sym=TRUE) {
    if(sym) {
      return( (1/Cpp.mean[t+1])*(Cpp.error[t+2] - Cpp.error[t] - Cpp.error[t+1]*corratio_val(t)))
    } else {
      return( (1/Cpp.mean[t+2])*(Cpp.error[t+1]-Cpp.error[t+2]/Cpp.mean[t+2]) )
    }
  } 

  uniroot_boot <- function(data,t_val,sym=TRUE) {
    if(sym) {
      val <- (mean(data[,t_val+2]) - mean(data[,t_val]))/ mean( data[,t_val+1] )
    } else {
      val <- (mean(data[,t_val+1])/mean(data[,t_val+2]))
    }
    solution <- uniroot(corratio_fun,sym=sym,value=val,t=t_val,N_t=T,lower=0,upper=3)
    return(solution$root)
  }

  uniroot_uwerr <- function(data,t_val,sym=TRUE) {
    val <- (data[1]-data[2])/data[3]
    solution <- uniroot(corratio_fun,sym=sym,value=val,t=t_val,N_t=T,lower=0,upper=3)
    return(solution$root)
  }

  noop_boot <- function(data) {
    return(data)
  }

  if(debug) {
    for(t in t1:t2) {
      print(sprintf("plotting 'f(E,t) - value(t)' for t = %d",t))
      plot(function(x) corratio_fun(x,value=corratio_val(t,sym=TRUE),t=t,N_t=T,sym=TRUE),xlim=c(-1,3),ylab="C(E,t)-value(t)",xlab="E")
      abline(h=0,col='gray')
      readline("press any key for next plot")
    }
  }

  # effective mass through solution of symmetric correlator ratio
  for(t in t1:t2) {
    print(t)
    
    # -- debug begin
    if(debug){
      #hist(Cpp.raw[,t+1],main=paste("Histogram of Cpp at",t),breaks=50)
      plot(density(Cpp.raw[,t+1]),main=paste("density plot of Cpp at",t))
      abline(v=mean(Cpp.raw[,t+1]))
      abline(v=median(Cpp.raw[,t+1]),col='red')
      #plot(density(tsboot(Cpp.raw[,t+1],noop_boot,sim="fixed",l=floor(NROW(Cpp.raw)/10),R=500)$t))
      readline()
    }
    # -- debug end
    
    solution <- uniroot(corratio_fun,value=corratio_val(t,sym=TRUE),sym=TRUE,t=t,N_t=T,lower=0,upper=6)
    m_eff <- c(m_eff,solution$root)
    # do bootstrap analysis for effective mass value for each timeslice
    boots <- tsboot(Cpp.raw,statistic=uniroot_boot,sim="fixed",l=floor(NROW(Cpp.raw)/10),R=1200,t_val=t,sym=TRUE)
    
    # -- debug begin
    #plot(density(boots$t),main="density of m_eff as measured on bootstrap")
    #readline()
    # -- debug end
    
    m_eff.error <- c(m_eff.error,sd(boots$t))
    m_eff.estimate <- c(m_eff.estimate,mean(boots$t))
    tempdat <- t(cbind(Cpp.raw[,t+2],Cpp.raw[,t],Cpp.raw[,t+1]))
    #print(tempdat)
    temp <- uwerr(uniroot_uwerr,data=tempdat,pl=FALSE,t_val=t)
    m_eff.estimate.uwerr <- c(m_eff.estimate.uwerr,temp$value)
    m_eff.error.uwerr <- c(m_eff.error.uwerr,temp$dvalue)
    #m_eff <- c(m_eff,log(Cpp.mean[t]/Cpp.mean[t+1]))
    #m_eff.error <- c(m_eff.error, abs(Cpp.error[t]/Cpp.mean[t])+abs(Cpp.error[t+1]/Cpp.mean[t+1])) 
  }

  # a slightly more elaborate fit which fits against a functional form
  picorrel <- function(t,a,E,N_t) {
    return( a*( exp(-t*E) + exp(-(N_t-t)*E) ) )
  }

  # from the range used for the effective mass estimate, derive a range in which
  # the correlator fit will be attempted directly
  picorfit.t1 <- t1
  picorfit.t2 <- t2
  # generate a data frame for this range
  Cppdata <- data.frame(t=seq(picorfit.t1,picorfit.t2,1),
    Cpp=Cpp.mean[(picorfit.t1+1):(picorfit.t2+1)],
    weight=(1/Cpp.error^2)[(picorfit.t1+1):(picorfit.t2+1)])
  
  # debug statements
  #print(Cppdata)
  #print(Cpp.mean)
  
  # plot raw correlator data points 
  cppmax <- max(Cpp.mean)
  cppmin <- min(Cpp.mean)
  # sometimes the correlator might become negative, here we catch that case
  # for the log plot
  if( cppmin < 0 ){
    cppmin <- 0.00000001
  }

  plotwitherror(x=seq(0,T-1,1),
    ylim=c(cppmin,1.2*cppmax),
    log="y"
    ,y=Cpp.mean,
    dy=Cpp.error,
    xlab="t",ylab=expression(C[pp](t)),
    xaxp=c(0,T-1,T-1))
    
  abline(v=picorfit.t1,col='gray')
  abline(v=picorfit.t2,col='gray')
  
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

  plot(xlim=c(0,T),function(x) picorrel(t=x,a=picorfit$a,E=picorfit$E,N_t=T),add=TRUE,col="red")
  plot(xlim=c(0,T),function(x) picorrel(t=x,a=picorfit$a-picorfit$da,E=picorfit$E-picorfit$dE,N_t=T),add=TRUE,col="pink")
  plot(xlim=c(0,T),function(x) picorrel(t=x,a=picorfit$a+picorfit$da,E=picorfit$E+picorfit$dE,N_t=T),add=TRUE,col="pink")

  label <- sprintf(" = %.5f (%.5f)",picorfit$E,picorfit$dE)
  text(x=T*2/3,y=cppmax*2/3,substitute(paste(am[eff],label)))

  # effective mass plot
  #meffmax <- max(m_eff)+max(m_eff.error)
  #meffmin <- min(m_eff)-max(m_eff.error)

  plotwitherror(x=seq(t1,t2,1),
    #ylim=c(meffmin,meffmax),
    y=m_eff.estimate.uwerr,dy=m_eff.error.uwerr,
    ylab=expression(am[eff](t)),xlab="t",
    xaxp=c(0,T-2,T-2),col='gray')

  plotwitherror(rep=TRUE,x=seq(t1,t2,1),
    y=m_eff,dy=m_eff.error)
      
  points(x=seq(t1,t2,1),y=m_eff.estimate,col='gray')

  m_ps <- mean(m_eff[(plateau.t1+1-t1):(plateau.t2+1-t1)])
  dm_ps <- sqrt(var(m_eff[(plateau.t1+1-t1):(plateau.t2+1-t1)])/(plateau.t2-plateau.t1)+mean(m_eff.error[(plateau.t1+1-t1):(plateau.t2+1-t1)]/sqrt(plateau.t2-plateau.t1))^2)
 
  print("pion mass from correlator fit")
  print(c(m_ps=picorfit$E,dm_ps=picorfit$dE))

  print("pion mass from explicit correlator ratio matching")
  print(c(m_ps=m_ps,dm_ps=dm_ps))

  abline(v=plateau.t1,col="gray")
  abline(v=plateau.t2,col="gray")

  abline(h=m_ps)
  abline(h=m_ps+dm_ps,col="gray")
  abline(h=m_ps-dm_ps,col="gray")

  abline(h=picorfit$E,col="red")
  abline(h=picorfit$E+picorfit$dE,col="pink")
  abline(h=picorfit$E-picorfit$dE,col="pink")
}
