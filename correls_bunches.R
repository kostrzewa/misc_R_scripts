correls_bunches <- function(dir,skip,debug=F) {
  filename <- sprintf("%s/piononline.dat",dir)
  omeas <- read.table(filename)
  Thalf <- max(omeas$V3)
  increment <- omeas[(omeas$V1==1 & omeas$V3==0),6][2]-omeas[(omeas$V1==1 & omeas$V3==0),6][1]
  max_sample <- max(omeas$V6)-1
  min_sample <- min(omeas$V6)+skip*increment
  pp_norm <- max(omeas[(omeas$V6>min_sample & omeas$V6<max_sample & omeas$V1==1),4])
  pp_min <- min(omeas[(omeas$V6>min_sample & omeas$V6<max_sample & omeas$V1==1),4])/pp_norm

  pa_norm <- max(omeas[(omeas$V6>min_sample & omeas$V6<max_sample & omeas$V1==2),4])

  if(debug) {
    print( sprintf("Measurements taken every %d trajectories",increment) )
    print( sprintf("PP normalization: %g",pp_norm) )
    print( sprintf("PP minimum: %g",pp_min) )
    print( sprintf("PA normalization: %g", pa_norm) )
  }
  
  
  step <- 10*increment

  # PP correlator in bunches of 10
  pdf(sprintf("%s_Pion_correls.pdf",dir),onefile=T,width=10,height=5,title=dir)
  par(mfrow=c(1,2),family="Palatino")
  for( i in seq(from=min_sample,to=max_sample,by=step) ) {
    if(debug) { print(i) }
    plot(x=omeas[(omeas$V6>i & omeas$V6<(i+step) & omeas$V1==1),3],
         y=omeas[(omeas$V6>i & omeas$V6<(i+step) & omeas$V1==1),4]/pp_norm,
         log='y',ylim=c(pp_min,1.2),
         main=paste("PP ",i),xlab="t",ylab=expression(C[PP])) 
    plot(x=omeas[(omeas$V6>i & omeas$V6<(i+step) & omeas$V1==2 & omeas$V3 > 1),3],
         y=omeas[(omeas$V6>i & omeas$V6<(i+step) & omeas$V1==2 & omeas$V3 > 1),4]/pa_norm,
         main=paste("PA ",i),xlab="t",ylab=expression(C[PA]))
         #,ylim=c(-0.02,0.04))
#         log='y') 
  }
  dev.off()

  # PA correlator in bunches of 10
  #pdf("PA_correls_bunches_of_10.pdf")
  #for( i in seq(from=min_sample,to=max_sample,by=10) ) {
  #  if(debug) { print(i) }
  #  plot(x=omeas[(omeas$V6>i & omeas$V6<(i+10) & omeas$V1==2),3],
  #  y=omeas[(omeas$V6>i & omeas$V6<(i+10) & omeas$V1==2),4],
  #  log='y',main=paste("PA ",i),xlab="t",ylab=expression(C[PA])) 
  #}
  #dev.off()
}
  