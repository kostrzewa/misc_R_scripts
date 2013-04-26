eval_Cppdata <- function(T,i_start,i_end,addon,t1,t2)
{
  Cpp.raw <- NULL
  for(i in i_start:i_end){
    filename <- sprintf("Cpp.data.%06d%s",i,addon)
    Cpp.data <- read.table(filename)
    Cpp.raw <- rbind(Cpp.raw,Cpp.data$V2)  
  }

  Cpp.mean <- NULL
  Cpp.error <- NULL
  m_eff <- NULL
  m_eff.error <- NULL
  for(t in 1:T){
    meanval <- mean(Cpp.raw[,t])
    errorval <- sqrt(var(Cpp.raw[,t]))/sqrt(i_end-i_start)
    Cpp.mean <- c(Cpp.mean,meanval)
    Cpp.error <- c(Cpp.error,errorval)
  }

  # simple exponential fit, taking m_eff to be the logarithm
  # of the forward derivative of C(t)
  for(t in 1:(T-1)) {
    m_eff <- c(m_eff,log(Cpp.mean[t]/Cpp.mean[t+1]))
    m_eff.error <- c(m_eff.error, abs(Cpp.error[t]/Cpp.mean[t])+abs(Cpp.error[t+1]/Cpp.mean[t+1])) 
  }

  # a slightly more elaborate fit which fits against a functional form
  picorrel <- function(t,a,E,N_t) {
    return( a*( exp(-t*E) + exp(-(N_t-t)*E) ) )
  }

  Cppdata <- data.frame(t=seq(6,14,1),Cpp=Cpp.mean[7:15])
  picorfit <- nls(Cpp~picorrel(t,a,E,20),data=Cppdata,weights=(1/Cpp.error^2)[7:15],start=list(a=1000,E=1))

  print(summary(picorfit))

  picorfit$a <- coef(summary(picorfit))["a","Estimate"]
  picorfit$da <- coef(summary(picorfit))["a","Std. Error"]
  picorfit$dE <- coef(summary(picorfit))["E","Std. Error"]
  picorfit$E <- coef(summary(picorfit))["E","Estimate"]

  #print(cbind(Cpp.mean,Cpp.error))

  cppmax <- max(Cpp.mean)
  cppmin <- min(Cpp.mean)
  if( cppmin < 0 ){
    cppmin <- 0.00001
  }

  par(mfrow=c(1,2))
  plotwitherror(x=seq(0,T-1,1),
    ylim=c(cppmin,1.2*cppmax),
    log="y",y=Cpp.mean,
    dy=Cpp.error,
    xlab="t",ylab=expression(C[pp](t)),
    xaxp=c(0,T-1,T-1))

  plot(xlim=c(0,20),function(x) picorrel(t=x,a=picorfit$a,E=picorfit$E,N_t=20),add=TRUE,col="red")
  plot(xlim=c(0,20),function(x) picorrel(t=x,a=picorfit$a-picorfit$da,E=picorfit$E-picorfit$dE,N_t=20),add=TRUE,col="pink")
  plot(xlim=c(0,20),function(x) picorrel(t=x,a=picorfit$a+picorfit$da,E=picorfit$E+picorfit$dE,N_t=20),add=TRUE,col="pink")
  #curve(y=picorrel(x,picorfit$a,picorfit$E,20),xlim=c(0,19))

  meffmax <- max(m_eff)
  meffmin <- min(m_eff)

  plotwitherror(x=seq(0,T-2,1),
    ylim=c(meffmin,meffmax),
    y=m_eff,dy=m_eff.error,
    ylab=expression(m[eff](t)),xlab="t",
    xaxp=c(0,T-2,T-2))

  m_ps <- mean(m_eff[(t1+1):(t2+1)])
  dm_ps <- sqrt((var(m_eff[(t1+1):(t2+1)])/(t2+t1))+mean(m_eff.error[(t1+1):(t2+1)]/sqrt(t2-t1))^2)

  print("pion mass from correlator derivative")
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
