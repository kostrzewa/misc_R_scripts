source("~/code/R/misc_R_scripts/lat_phys_ratios/compute_ratio.R")

plot_overview <- function(file,pdf.filename="csw_overview.pdf") {
  cswdat <- read.table(file,header=T)

  mps <- list(val=cswdat$mps, dval=cswdat$dmps, name="Mps")
  fps <- list(val=cswdat$fps, dval=cswdat$dfps, name="fps")

  mps_ov_fps <- compute_ratio(dividend=mps, divisor=fps, name=expression(M[ps]/f[ps]))

  pdf(file=pdf.filename,width=6,height=6,onefile=T)
  plotwitherror(x=cswdat$csw, y=mps_ov_fps$val, dy=mps_ov_fps$dval, pch=cswdat$sym,
                ylab=mps_ov_fps$name, xlab=expression(C[sw]), main="nf=2, beta=2.1, amu=0.006")
  plotwitherror(x=cswdat$csw, y=mps$val, dy=mps$dval, pch=cswdat$sym,
                ylab=expression(M[ps]), xlab=expression(C[sw]), main="nf=2, beta=2.1, amu=0.006")
  plotwitherror(x=cswdat$csw, y=fps$val, dy=fps$dval, pch=cswdat$sym,
                ylab=expression(f[ps]), xlab=expression(C[sw]), main="nf=2, beta=2.1, amu=0.006")
  dev.off()
}
