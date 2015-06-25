load("extrapolations.iwa_b2.1_L48T96_k0.13729_mul0.0009.Rdata")
load("extrapolations.sys.iwa_b2.1_L48T96_k0.13729_mul0.0009.Rdata")
physvales <- read.csv("phys_values.csv",stringsAsFactors=TRUE)

source("~/code/R/misc_R_scripts/plotutils/extract_data_summary_listplot.R");
lplotdat <- extract_data_summary_listplot(ext=extrapolations.iwa_b2.1_L48T96_k0.13729_mul0.0009,
                                          ext.sys=extrapolations.sys.iwa_b2.1_L48T96_k0.13729_mul0.0009, 
                                          physvalues=physvalues,
                                          a=list(val=0.0914,dval=0.0022))

source("~/code/R/misc_R_scripts/ratios_and_interpolations/reshuffle_ext_table.R")
lplotdat <- reshuffle_ext_table(lplotdat)

source("~/code/R/misc_R_scripts/plotutils/plot_list_vertical.R");

for(errsum.method in c("linear","linear.quadrature","quadrature")){
   
  plot_list_vertical(basename=sprintf("Qlat_ov_Qphys.%s",errsum.method),
                     x=lplotdat$x,dx=lplotdat$dx,mdx=lplotdat$mdx,
                     labels=lplotdat$labels,labelpos=0.86,
                     errsum.method=errsum.method,
                     item.height=0.4,mar=c(5,5,5,5),pch='.',cex=3.0, 
                     xlim=c(0.89,1.15),ylim=c(0.5,length(lplotdat$x)+0.5 ),
                     xlab="$Q_\\mathrm{lat}/Q_\\mathrm{phys}$")
}
