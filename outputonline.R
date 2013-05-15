outputonline <- function(type,beta,L,T,kappa,mu,t1,t2,skip,csw=0)
{
  rundir <- NULL
  if(csw==0) {
    rundir <- sprintf("%s_b%g-L%dT%d-k%.7g-mu%g",type,beta,L,T,kappa,mu)
  } else {
    rundir <- sprintf("%s_b%g-L%dT%d-csw%g-k%g-mu%g",type,beta,L,T,csw,kappa,mu)
  }
  filename <- sprintf("%s/piononline.dat",rundir)
  plaqfile <- sprintf("%s/output.data",rundir)
  pioncor <- readcmicor(filename)
  onlineout <- onlinemeas(pioncor,t1=t1,t2=t2,kappa=kappa,mu=mu,skip=skip,method="uwerr")
  par(mfrow=c(1,2))
#  print("plotting dpaopp")
  plot(x=seq(skip+1,(skip+length(onlineout$MChist.dpaopp)),1),onlineout$MChist.dpaopp,t='l',ylab=expression(am[PCAC]),xlab=expression(t[HMC]))
  abline(h=onlineout$uwerrresultmpcac$value,col="black")
  abline(h=mean(onlineout$MChist.dpaopp),col="red")
  lengthdpaopp <- length(onlineout$MChist.dpaopp)
  mindpaopp <- min(onlineout$MChist.dpaopp)
  maxdpaopp <- max(onlineout$MChist.dpaopp)
  legend( 
    x=skip+0.5*lengthdpaopp,
    y=mindpaopp-sign(mindpaopp)*(maxdpaopp-mindpaopp)*0.17,
    lwd=c(1,1),col=c("black","red"),
    legend=c(
    expression(paste(m[PCAC]," from fit")),
    expression(paste(m[PCAC]," from derivative"))
    ),bg="white"
  )

#  print("plotting plaquette")
  # something in the skip computation is odd, let's just solve it like this
  if(skip==0){
    shift <- 1
  } else {
    shift <- 0
  }
  plot(x=seq(skip+shift,length(read.table(plaqfile)$V2),1),read.table(plaqfile)$V2[skip:length(read.table(plaqfile)$V2)],t='l',ylab=expression("<P>"),xlab=expression(t[HMC])) 
  print("m_pcac from fit")
  print(c(mcpac=onlineout$uwerrresultmpcac$value,dmpcac=onlineout$uwerrresultmpcac$dvalue))
  print("m_pcac from derivative")
  print(c(mpcac=mean(onlineout$MChist.dpaopp),dmpcac=sqrt(var(onlineout$MChist.dpaopp)/length(onlineout$MChist.dpaopp))))

  print("best estimate")
  print(c(mpcac=mean(c(onlineout$uwerrresultmpcac$value,mean(onlineout$MChist.dpaopp))),
        dmpcac=sqrt(onlineout$uwerrresultmpcac$dvalue^2+sqrt(var(onlineout$MChist.dpaopp)/length(onlineout$MChist.dpaopp))^2)))
}
