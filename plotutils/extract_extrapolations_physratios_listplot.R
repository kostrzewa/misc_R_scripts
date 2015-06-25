# given extrapolations and extrapolations.sys from the ratios_and_interpolations
# analysis, extract the various values from the last row for each observable
# as this corresponds to the situation of interest

source("~/code/R/misc_R_scripts/prop_error/binary_operations.R")

extract_extrapolations_physratios_listplot <- function(ext, ext.sys, physratios) {
  if( length(ext) != length(ext.sys) ){
    stop("extract_interpolations_listplot: ext and ext.sys should be of the same length!\n")
  }
  
  names <- c(#"f_K", "f_D", "f_Ds", "m_Ds",
             "m_K_ov_f_K","m_D_ov_f_D","m_Ds_ov_f_Ds",
             "f_K_ov_f_pi","f_D_ov_f_pi","f_Ds_ov_f_pi",
             "f_D_ov_f_K","f_Ds_ov_f_K",
             "f_Ds_ov_f_D")
  
  x <- NULL
  dx <- NULL
  mdx <- NULL
  labels <- NULL
  for(i in 1:length(ext)){
    if(!any(names==ext[[i]]$name[1])) next
    phys <- physratios[ext[[i]]$name[1]==physratios$name,]
    phys.noerr <- phys
    phys.noerr$dval <- 0 
    if(nrow(phys)==0) next
    nr <- nrow(ext[[i]])
    labels <- c(labels,ext[[i]]$texlabel[nr])
    ratio <- compute_ratio(ext[[i]][nr,],phys,name=phys$name)
    ratio.noerr <- compute_ratio(ext[[i]][nr,],phys.noerr,name=phys$name)
    x <- c(x,ratio$val)
    noerr <- ratio.noerr$dval
    err <- ratio$dval - ratio.noerr$dval
    dx <- rbind(dx,cbind(noerr,ext.sys[[i]]$serr[nr]/phys$val,err))
    mdx <- rbind(mdx,cbind(noerr,ext.sys[[i]]$mserr[nr]/phys$val,err))
  }
  list(labels=labels,x=x,dx=dx,mdx=mdx)
}

