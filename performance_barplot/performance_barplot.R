performance_barplot <- function(datafile,debug=F,offset=20000,xmax=65,title="") {
  require(RColorBrewer)
  performance <- read.table(header=T,file=datafile)
  performance <- performance[ order(performance$comm),  ]

  if(debug) {
    print(performance)
  }
  
  vlength <- length(performance$name)+1
  names <- NULL
  for( i in 1:length(performance$name) ) {
    names[i] <- sprintf("(%d) %s",vlength-i,performance$name[i])
    if( performance$comm[i] == 0 ) {
      performance$nocomm[i] <- performance$nocomm[i]
      performance$comm[i] <- offset
    }
    if(debug) {
      print(names[i])
    }
  }

  cols <- brewer.pal(n=8,name="Paired")
  #cols <- heat.colors(n=8)
  par(mar=c(4,15.1,3,2),family="Palatino",mgp=c(1.0,0.3,0))
  values <- matrix(c(performance$comm-offset,performance$nocomm-performance$comm),ncol=2)
  if(debug) {
    print(values)
  }
  barplot(xlim=c(offset/1000,xmax),offset=offset/1000,t(values)/1000,names=names,col=rev(cols[1:2]),horiz=T,las=1,xlab="GFlop/s per node",border=NA,main=title)
}
