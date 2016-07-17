study_relerr <- function(basename,nosamples=12,T=48,boot.R=1000,boot.l=1,seed=123456,sym=TRUE,label="C_\\pi",log="",plot.cf=FALSE){
  files <- getorderedfilelist(path=".",basename=basename, last.digits=4)
  cfsamples <- readbinarysamples(files=files, T=T, nosamples=nosamples, sym=sym )
  for( i in 1:length(cfsamples) ){
    cfsamples[[i]] <- bootstrap.cf(cf=cfsamples[[i]],boot.R=boot.R,boot.l=boot.l,seed=seed)
  }

  plotname <- strsplit(x=basename,split='/')[[1]]
  plotname <- plotname[length(plotname)]

  tikzfiles <- tikz.init(basename=sprintf("relerrs.%s",plotname),width=5,height=5)
  if(plot.cf){
    for( i in 1:length(cfsamples) ){
      plot(cfsamples[[i]],xlab="$t/a$",log=log,ylab=sprintf("$%s(t/a,N_\\eta=%d)$",label,i),xaxp=c(0,T/2+1,T/2+1))
    }
  }
  
  t20relerr <- double()
  cferr <- NULL
  for( i in 1:length(cfsamples) ){ 
    cferr <- apply(X=cfsamples[[i]]$cf.tsboot$t,MARGIN=2,FUN=sd)
    cferr.uw.list <- apply(X=cfsamples[[i]]$cf,FUN=function(x){ uwerrprimary(data=x,pl=FALSE,S=3) }, MARGIN=2)
    cferr.uw <- NULL
    for(uw in cferr.uw.list){
      cferr.uw <- rbind(cferr.uw,data.frame(value=uw$value,dvalue=uw$dvalue,ddvalue=uw$ddvalue))
    }
    plot(ylim=c(0,0.18),y=cferr/cfsamples[[i]]$cf0,x=0:(T/2),las=1,xlab="$t/a$",ylab=sprintf("$\\sigma_{%s}/%s\\,(t/a,N_\\eta=%d)$",label,label,i),pch=15,xaxp=c(0,T/2+1,T/2+1))
    plotwitherror(x=0:(T/2),y=cferr.uw$dvalue/cferr.uw$value,dy=cferr.uw$ddvalue/cferr.uw$value,rep=TRUE,col="red",pch=16,cex=0.5)
    legend("topright",col=c("black","red"),pch=c(15,16),bty='n',legend=c("bootstrap","Gamma method"))
    t20relerr <- c(t20relerr,cferr[21]/cfsamples[[i]]$cf0[21])
  }

  etafit <- lm( "t20relerr ~ oneoversqrtNeta + 0", data=data.frame(t20relerr=t20relerr,oneoversqrtNeta=1/sqrt(1:12)) )

  print(etafit)
  
  Neta <- seq(from=1,to=15,by=0.2)
  newdata <- data.frame(oneoversqrtNeta=1/sqrt(Neta))

  pred.y <- predict(etafit,newdata=newdata)


  plot(y=t20relerr,x=1:nosamples,xlab="$N_\\eta$",ylab=sprintf("$\\sigma_{%s}/%s\\,(t/a=20,N_\\eta)$",label,label),las=1,xaxp=c(1,nosamples,nosamples-1),pch=16)
  lines(x=Neta,y=pred.y)
  legend("topright",lty=1,pch=NA,col="black",legend=sprintf("$%f/\\sqrt{N_\\eta}$",etafit$coefficients[1]),bty='n',lwd=2)

  tikz.finalize(tikzfiles)
}

