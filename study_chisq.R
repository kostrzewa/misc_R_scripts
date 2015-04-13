study_chisq <- function(file) {
  chisq <- read.table(header=T,file=file)
  pdf("study_chisq.pdf")
  par( mfrow = c( 2, 2 ) )
  plotwitherror(main="matrixfit mass", xlab="boot.l", ylab="aM",
                x=chisq$boot.l, y=chisq$M, dy=chisq$dM, col=chisq$col, pch=chisq$pch+chisq$pchoffset)
  plotwitherror(main="matrixfit reduced chisq", xlab="boot.l", ylab="chisq/df",
                x=chisq$boot.l, y=chisq$M.chisq, col=chisq$col, pch=chisq$pch+chisq$pchoffset)
  plotwitherror(main="matrixfit Q", xlab="boot.l", ylab="Q",
                x=chisq$boot.l, y=chisq$M.Q, col=chisq$col, pch=chisq$pch+chisq$pchoffset,ylim=c(0,1))
  plot.new()

  par( mfrow = c( 2, 2 ) )
  plotwitherror(main="effectivemass", xlab="boot.l", ylab="aM",
                x=chisq$boot.l, y=chisq$M, dy=chisq$dM, col=chisq$col, pch=chisq$pch+chisq$pchoffset)
  plotwitherror(main="effectivemass reduced chisq", xlab="boot.l", ylab="chisq/df",
                x=chisq$boot.l, y=chisq$M.chisq, col=chisq$col, pch=chisq$pch+chisq$pchoffset)
  plotwitherror(main="effectivemass Q", xlab="boot.l", ylab="Q",
                x=chisq$boot.l, y=chisq$Meff.Q, col=chisq$col, pch=chisq$pch+chisq$pchoffset, ylim=c(0,1))
  plot.new()

  dev.off()

}
