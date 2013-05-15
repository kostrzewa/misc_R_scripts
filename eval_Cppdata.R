eval_Cppdata <- function(prefix,T,i_start,i_end,addon,t1,t2,debug=FALSE)
{
  Tover2 <- floor(T/2)

  Cpp.raw <- NULL
  Cpp.lengths <- NULL

  # read the data files
  for(i in i_start:i_end){
    filename <- sprintf("%s/Cpp.data.%06d%s",prefix,i,addon)
    Cpp.data <- read.table(filename)$V2
    # in the case where only half the time extent has been measured (online f.ex.)
    # mirror the correlator  
    if( length(Cpp.data) < T ) {
      for(t in 0:(Tover2-1)) {
        Cpp.data <- c(Cpp.data,Cpp.data[Tover2-t])
      }
    }
    Cpp.raw <- rbind(Cpp.raw,Cpp.data)
  }

  Cpp.mean <- NULL
  Cpp.error <- NULL
  m_eff <- NULL
  m_eff.error <- NULL
  for(t in 1:T){
    meanval <- mean(Cpp.raw[,t])
    errorval <- sqrt(var(Cpp.raw[,t])/(i_end-i_start))
    Cpp.mean <- c(Cpp.mean,meanval)
    Cpp.error <- c(Cpp.error,errorval)
  }
  
  # function to solve for effective mass in symmetric manner
  # the 'value' input is for finding the root of f(E,t) - value(t) = 0
  corratio_fun <- function(E,value,t,N_t,sym=TRUE) {
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

  uniroot_boot <- function(data,indices,t,sym=TRUE) {
    val <- (mean(data[,t+2][indices]) - mean(data[,t][indices]))/ mean( data[,t][indices] )
    solution <- uniroot(corratio_fun,value=val,t=t,N_t=T,lower=0,upper=3)
    return(solution$root)
  }

  if(debug) {
    for(t in t1:t2) {
      print(sprintf("plotting 'f(E,t) - value(t)' for t = %d",t))
      plot(function(x) corratio_fun(x,value=corratio_val(t,sym=TRUE),t=t,N_t=T,sym=TRUE),xlim=c(-1,3),ylab="C(E,t)-value(t)",xlab="E")
      readline("press any key for next plot")
    }
  }

  # effective mass through solution of symmetric correlator ratio
  for(t in t1:t2) {
    solution <- uniroot(corratio_fun,value=corratio_val(t,sym=TRUE),sym=TRUE,t=t,N_t=T,lower=0,upper=6)
    m_eff <- c(m_eff,solution$root)
    # do bootstrap analysis for effective mass value for each timeslice
    boots <- boot(Cpp.raw,uniroot_boot,R=200,t=t)
    m_eff.error <- c(m_eff.error,sd(boots$t))
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

  par(mfrow=c(1,2))
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
  text(x=T*2/3,y=cppmax*2/3,substitute(paste(m[eff],label)))

  # effective mass plot
  meffmax <- max(m_eff)+max(m_eff.error)
  meffmin <- min(m_eff)-max(m_eff.error)

  plotwitherror(x=seq(t1,t2,1),
    ylim=c(meffmin,meffmax),
    y=m_eff,dy=m_eff.error,
    ylab=expression(am[eff](t)),xlab="t",
    xaxp=c(0,T-2,T-2))

  abline(v=t1,col='gray')
  abline(v=t2,col='gray')

  m_ps <- mean(m_eff)
  dm_ps <- sqrt((var(m_eff)/(length(m_eff)))+mean(m_eff.error/sqrt(length(m_eff.error)))^2)

  print("pion mass from explicit correlator ratio matching")
  print(c(m_ps=m_ps,dm_ps=dm_ps))

  print("pion mass from correlator fit")
  print(c(m_ps=picorfit$E,dm_ps=picorfit$dE))

  abline(h=m_ps)
  abline(h=m_ps+dm_ps,col="gray")
  abline(h=m_ps-dm_ps,col="gray")

  abline(h=picorfit$E,col="red")
  abline(h=picorfit$E+picorfit$dE,col="pink")
  abline(h=picorfit$E-picorfit$dE,col="pink")
}
