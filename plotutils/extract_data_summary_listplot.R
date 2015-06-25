extract_data_summary_listplot <- function(ext, ext.sys, physvalues, a) {
  hbarc <- 197.3269718
  mN <- list(val=0.4398,dval=0.0038)
  # average neutron/proton mass
  mN.pdg <- list(val=938.9186795,dval=0.00000023)

  if( length(ext) != length(ext.sys) ){
    stop("extract_interpolations_listplot: ext and ext.sys should be of the same length!\n")
  }
 
  value.labels <- c("$f_K$","$f_D$","$f_{D_s}$","$M_{D_s}$")
  values <- c("f_K", "f_D", "f_Ds", "m_Ds")
  ratios <- c("m_K_ov_f_K","m_D_ov_f_D","m_Ds_ov_f_Ds",
             "f_K_ov_f_pi","f_D_ov_f_pi","f_Ds_ov_f_pi",
             "f_D_ov_f_K","f_Ds_ov_f_K",
             "f_Ds_ov_f_D")
  
  x <- NULL
  dx <- NULL
  mdx <- NULL
  labels <- NULL
  errband <- list(xleft=c(),xright=c())
  
  # in the current form the statistical and systematic errors are added in quadrature and stored in err.msys and err.sys
  # the error coming from the lattice spacing is given as a separate error for dimensional quantities
  # the error coming from the numerator is 
    
  for(name in ratios) {
    for(i in 1:length(ext)){
      if(ext[[i]]$name != name) next
      # the relevant extrapolated value is always the last one
      nr <- nrow(ext[[i]])
      phys <- physvalues[ext[[i]]$name[1]==physvalues$name,]

      err.stat <- 0 #ext[[i]]$dval[nr]/phys$val
      err.sys <- sqrt( (ext[[i]]$dval[nr]/phys$val)^2 + (ext.sys[[i]]$serr[nr]/phys$val)^2 ) #ext.sys[[i]]$serr[nr]/phys$val 
      err.msys <- sqrt( (ext[[i]]$dval[nr]/phys$val)^2 + (ext.sys[[i]]$mserr[nr]/phys$val)^2 ) #ext.sys[[i]]$mserr[nr]/phys$val
      err.a <- 0.0
      err.phys <- 0 #phys$dval*ext[[i]]$val[nr]/phys$val^2

      x <- c(x,ext[[i]]$val[nr]/phys$val)
      dx <- rbind(dx,cbind(err.stat,err.sys,err.a,err.phys))
      mdx <- rbind(mdx,cbind(err.stat,err.msys,err.a,err.phys))
      errband <- rbind(errband,data.frame( xleft=(1-(phys$dval*ext[[i]]$val[nr]/phys$val^2)), xright=(1+(phys$dval*ext[[i]]$val[nr]/phys$val^2)) ) ) 
      labels <- c(labels,ext[[i]]$texlabel[nr])
    }
  }
  
  for(name in values) {
    for(i in 1:length(ext)){
      if(ext[[i]]$name != name) next
      # the relevant extrapolated value is always the last one
      nr <- nrow(ext[[i]])
      phys <- physvalues[ext[[i]]$name[1]==physvalues$name,]
     
      err.stat <- 0 #hbarc*ext[[i]]$dval[nr]/phys$val/a$val
      err.sys <- sqrt( (hbarc*ext[[i]]$dval[nr]/phys$val/a$val)^2 + (hbarc*ext.sys[[i]]$serr[nr]/phys$val/a$val)^2 ) #hbarc*ext.sys[[i]]$serr[nr]/phys$val/a$val
      err.msys <- sqrt( (hbarc*ext[[i]]$dval[nr]/phys$val/a$val)^2 + (hbarc*ext.sys[[i]]$mserr[nr]/phys$val/a$val)^2 ) #hbarc*ext.sys[[i]]$mserr[nr]/phys$val/a$val
      err.a <- hbarc*a$dval*ext[[i]]$val[nr]/phys$val/(a$val^2)
      err.phys <- 0 #hbarc*phys$dval*ext[[i]]$val[nr]/(phys$val^2)/a$val
      
      x <- c(x,hbarc*ext[[i]]$val[nr]/a$val/phys$val)
      dx <- rbind(dx,cbind(err.stat,err.sys,err.a,err.phys))
      mdx <- rbind(mdx,cbind(err.stat,err.msys,err.a,err.phys))
      val.idx <- which(values==ext[[i]]$name[1])
      errband <- rbind(errband,data.frame(xleft=(1-(hbarc*phys$dval*ext[[i]]$val[nr]/(phys$val^2)/a$val)),xright=(1+(hbarc*phys$dval*ext[[i]]$val[nr]/(phys$val^2)/a$val)) ) )
      labels <- c(labels,value.labels[val.idx])
       
    } 
  }
  # add the nucleon
  err.stat <- hbarc*mN$dval/mN.pdg$val/a$val
  err.sys <- 0.0
  err.msys <- 0.0
  err.a <- hbarc*mN$val*a$dval/mN.pdg$val/(a$val^2)
  err.phys <- 0 #hbarc*mN.pdg$dval*mN$val/phys$val^2/a$val

  x <- c(x,hbarc*mN$val/mN.pdg$val/a$val)
  dx <- rbind(dx,cbind(err.stat,err.sys,err.a,err.phys))
  mdx <- rbind(mdx,cbind(err.stat,err.msys,err.a,err.phys))
  errband <- rbind(errband, data.frame(xleft=(1-(hbarc*mN.pdg$dval*mN$val/phys$val^2/a$val)),xright=(1+(hbarc*mN.pdg$dval*mN$val/phys$val^2/a$val)) ) )
  labels <- c(labels,"$M_N$")
   
  print(data.frame(xleft=x-err.phys,xright=x+err.phys))
   
  list(labels=labels,x=x,dx=dx,mdx=mdx,errband=errband )
}

