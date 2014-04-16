source("~/code/R/misc_R_scripts/plot_timeseries.R")
plot_forces <- function(path) {
  files <- list.files(path,pattern="*.data")
  for( i in seq(1,length(files))) {
    print(files[i])
    name <- strsplit(files[i],split=c("."),fixed=T)
    forcedat <- read.table(files[i])
    offset <- 0
    for( measure in c("avg","max") ) {
      titletext <- sprintf("%s.%s.force",name[[1]][1],measure)
      pdf.filename <- sprintf("%s.pdf",titletext)
      plot_timeseries(dat=forcedat[100:max(forcedat[,1]),(3+offset)],trange=c(100,max(forcedat[,1])),ylab=paste(measure,"F^2"),pdf.filename=pdf.filename,
      name=paste(measure,"F^2"), plotsize=5,filelabel=titletext,titletext=titletext)
      offset <- offset + 1
    }
  }
}
