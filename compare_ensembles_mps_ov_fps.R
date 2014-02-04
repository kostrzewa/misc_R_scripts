compute_ratio <- function(dividend,divisor,debug=FALSE) {                                                                             
  ratio <- list( value=dividend$value / divisor$value, 
                 error=sqrt( (dividend$error/divisor$value)^2 + (divisor$error*dividend$value/divisor$value^2)^2 ))
  if(debug) {
    print(ratio)
  }
  return(ratio)
}

compare_ensembles_mps_ov_fps <- function(datafile,debug=FALSE) {
  dat <- read.table(datafile,header=T)
  dividends <- list(value=dat$m_ps,error=dat$dm_ps)
  divisors <- list(value=dat$f_ps,error=dat$df_ps)
  ratios <- compute_ratio(dividend=dividends,divisor=divisors,debug=T)
  
  require(tikzDevice) 
  tikz('m_ps_over_f_ps.tex', standAlone = TRUE, width=4, height=4)

  plotwitherror(x=seq(1,length(ratios$value)),
    y=ratios$value,
    dy=ratios$error,xaxt='n',
    main="$\\frac{m_{PS}}{f_{PS}}$ for $N_f=2+2$ analyses at $\\beta=1.85$ and $N_f=2$ at $\\beta=4.35$",
    ylab="",
    xlab="",las=1,cex.main=0.85)

  axis(1,at=seq(1,length(ratios$value)),labels=dat$run)

  dev.off()
  tools::texi2dvi("m_ps_over_f_ps.tex",pdf=T)
#  system("pdfcrop m_ps_over_f_ps.pdf m_ps_over_f_ps.pdf")
}
