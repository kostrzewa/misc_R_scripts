
# compare_timeseries is a function which attempts to automatically read out
# a number of HMC histories from "output.data" or "monomial-xx.data" files
# in different subdirectories of a directory path

# It would be called via the "compare_driver" function, where this function takes
# only a path argument (and an optional "debug")

# In order to use it, one has to adjust the legend labels in the correct order
# for which it helps to run with "debug=T" once.

# "observables", "files", "columns" and "ylims" are vectors and data.frames respectively
# which hold information about the timeseries to be extracted and plotted
# the number of rows in the data frame and the length of the vectors must be equal

# observables: names of the timeseries to be extracted. This is also used as a filename
#              for the resulting pdf
# files      : files which are to be read to extract a given observable
# columns    : columns in the given file (c1: x-axis, c2: y-axis) which are to be extracted for the given observable
# ylims      : lower and upper y plot boundaries for the resulting plot

compare_driver <- function(path,debug=F) {
  leg.labels <- c("poly O(192) cold", "poly O(100) continue", "poly O(192) hot", "rational O(12) continue", "poly O(192) mu=0.006, cold")
  observables <- c("<P>","dH","exp(-dH)","minEV","maxEV")
  files <- c(rep("output.data",3),rep("monomial-04.data",2))
  columns <- data.frame(c1=c(1,1,1,1,1),c2=c(2,3,4,2,3))
  ylims <- data.frame(y1=c(0.532,-3,0,1.6e-5,0.62),y2=c(0.536,6,10,4e-5,0.9))

  for( i in 1:length(observables) ) {
    pdf(file=paste(observables[i],".pdf",sep=""),height=6,width=6.5,title=observables[i])
    par(family="Palatino")
    compare_timeseries(path=path,filename=files[i],name=observables[i],
                       columns=c(columns$c1[i],columns$c2[i]),ylim=c(ylims$y1[i],ylims$y2[i]),
                       debug=debug, leg.labels=leg.labels)
    dev.off()
  }
}

# compare_timeseries will enter "path" and list all subdirectories contained within that. 
# It will subsequently attempt to descend into all these subdirectories, reading the file
# given by the "filename" argument
# the "columns" argument refers to the columns in the file given by "filename" where 
# columns[1] will be the x-axis and columns[2] the y-axis in the resulting plots

compare_timeseries <- function(path,filename,ylim,leg.labels,columns=c(1,2),name="<P>",debug=F) {
  if(debug) {
    print(path)
    print(filename)
    print(cols)
  }
  dirs <- list.dirs(path,recursive=F)
  files <- sprintf("%s/%s",dirs,filename)
  if(debug){
    print(files)
  }
  library("RColorBrewer")
  colors <- brewer.pal(n=length(files),name="Set1")
  ts <- NULL
  xmax <- 0
  for( i in 1:length(files) ) {
    if(debug){
      cat("Loading file", i, "\n")
    }
    no_columns <- max(count.fields(files[i]))
    ts[[i]] <- read.table(files[[i]],fill=TRUE,col.names=1:no_columns)
    if( max(ts[[i]][,1]) > xmax ) {
      xmax <- max(ts[[i]][,1])
    }
  }
  
  plot(NA,xlim=c(0,xmax),ylim=ylim,xlab=expression(t[HMC]),ylab=name,main=name)
  for( i in 1:length(files) ) {
    lines(x=ts[[i]][,columns[1]],y=ts[[i]][,columns[2]],col=colors[i])
  }
 
  legend(x=xmax/3,y=ylim[2],col=colors,legend=leg.labels,lty=rep(1,length(files)),bg="white")
       
}
