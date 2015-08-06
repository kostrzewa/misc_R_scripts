source("~/code/R/misc_R_scripts/plot_timeseries.R")
plot_forces <- function(path,skip=0,debug=FALSE) {
  files <- list.files(path,pattern="*.force.data",full.names=T)
  for( i in seq(1,length(files))) {
   
    # extract some strings from the filename 
    print(files[i])
    filename <- strsplit(files[i],split=c("/"),fixed=T)
    name <- filename[[length(filename)]][1]
    monomial <- filename[[length(filename)]][ length( filename[[length(filename)]] ) ]
    monomial <- strsplit(monomial,split=c("."), fixed=T)
    monomial <- monomial[[1]][1]
    if(debug) {
      cat("Extracting name and monomial name from filename\n")
      cat("Filename:", files[i], "\n")
      cat("Name:", name, "\n")
      cat("Monomial:", monomial, "\n")
    }

    forcedat <- read.table(files[i])
    offset <- 0
    for( measure in c("avg","max") ) {
      titletext <- sprintf("%s.%s.force",monomial,measure)
      pdf.filename <- sprintf("%s.%s.pdf",titletext,name)
      if(debug) {
        cat("PDF filename:", pdf.filename, "\n")
      }
      plot_timeseries(dat=forcedat[skip:length(forcedat[,1]),(2+offset)],trange=c((skip),length(forcedat[,1])),ylab=paste(measure,"F^2"),pdf.filename=pdf.filename,
      name=paste(measure,"F^2"), plotsize=5,filelabel=titletext,titletext=titletext)
      offset <- offset + 1
    }
  }

  if( Sys.which("pdfcat") != "" ) {
    command <- sprintf("pdfcat forces.%s.pdf *.force.%s.pdf",name,name)
    print(paste("calling",command))
    system(command=command)
  } else {
    print("pdfcat not found, not concatenating plots!")
  }
  
}
