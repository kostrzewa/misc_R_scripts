plot_omeas_resultsum <- function(debug=FALSE){
  resultsfile <- "omeas.summary.RData"
  resultsum <- list()
  if(file.exists(resultsfile)){
    load(resultsfile)
  } else {
    stop("plot_omeas_resultsum: 'omeas.summary.RData' could not be found!")
  }

  keys <- names(resultsum)
  observables <- names(resultsum[[1]]$obs)

  # change ordering of data in the list, such that the observable is "the slowest running index"
  res <- list()
  for( obs in observables ){
    res[[obs]] <- data.frame()
    for( key in keys ){
      res[[obs]] <- rbind( res[[obs]], data.frame(key=key, data.frame(t(resultsum[[key]]$obs[,obs])) ) )
    }
  }
  if(debug) print(res); print(res[["mpi"]]$val)

  tikzfiles <- tikz.init(basename="omeas.summary", width=4.5, height=4.5, sanitize=TRUE)
  op <- par(mar=c(6,4,1,1))
  for( obs in observables ){
    N <- length( res[[obs]]$val )
    plotwitherror(x=1:N, y=res[[obs]]$val, dy=res[[obs]]$dval, ylab=obs, xlab="",
                  xaxt='n')
    axis(side=1, at=1:N, labels=NA)
    text(x=1:N, labels=keys, y=par()$usr[3]-0.20*(par()$usr[4]-par()$usr[4]), srt=25, adj=1.0, xpd=TRUE )
  }
  tikz.finalize(tikzfiles)
  par(op)
}
