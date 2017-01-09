extract_data_summary_listplot <- function(ext, ext.sys, physvalues, a) {
  hbarc <- 197.3269718
  mN <- list(val=0.440,dval=0.004)
  # average neutron/proton mass
  mN.pdg <- list(val=938.9186795,dval=0.00000023)

  if( length(ext) != length(ext.sys) ){
    stop("extract_interpolations_listplot: ext and ext.sys should be of the same length!\n")
  }
 
  #value.labels <- c("$f_K$","$f_D$","$f_{D_s}$","$M_{D_s}$")
  #values <- c("f_K", "f_D", "f_Ds", "m_Ds")
  #ratios <- c("m_K_ov_f_K","m_D_ov_f_D","m_Ds_ov_f_Ds",
  #           "f_K_ov_f_pi","f_D_ov_f_pi","f_Ds_ov_f_pi",
  #           "f_D_ov_f_K","f_Ds_ov_f_K",
  #           "f_Ds_ov_f_D")
  
  value.labels <- c("$f_K$","$f_{D_s}$","$M_{D_s}$")
  values <- c("f_K","f_D","f_Ds")
  ratios <- c("f_K_ov_f_pi","f_D_ov_f_pi","f_Ds_ov_f_pi","f_Ds_ov_f_K","f_Ds_ov_f_D")
  ratios <- ratios[c(5,4,3,2,1)]
  

  # these will hold relative quantities for the summary plot with partly combined errors
  x <- NULL
  dx <- NULL
  mdx <- NULL
  labels <- NULL
  errband <- list(xleft=c(),xright=c())

  # these the quantities themselves with individual errors
  X <- NULL
  DX <- NULL
  MDX <- NULL
  
  # in the current form the statistical and systematic errors are added in quadrature for the summary plot
  # the error coming from the lattice spacing is given as a separate error for dimensional quantities
  # the error coming from the numerator is 
    
  for(name in ratios) {
    for(i in 1:length(ext)){
      if(ext[[i]]$name != name) next
      # the relevant extrapolated value is always the last one
      nr <- nrow(ext[[i]])
      phys <- physvalues[ext[[i]]$name[1]==physvalues$name,]

      err.stat <- ext[[i]]$dval[nr]/phys$val
      err.sys <- ext.sys[[i]]$serr[nr]/phys$val 
      err.msys <- ext.sys[[i]]$mserr[nr]/phys$val
      err.a <- 0.0
      err.phys <- phys$dval*ext[[i]]$val[nr]/phys$val^2

      x <- c(x,ext[[i]]$val[nr]/phys$val)
      X <- c(X,ext[[i]]$val[nr])
      dx <- rbind(dx,cbind( 0.0, sqrt(err.stat^2 + err.sys^2), err.a, 0.0 ))
      DX <- rbind(DX,phys$val*cbind( err.stat, err.sys, err.a ))
      mdx <- rbind(mdx,cbind( 0.0, sqrt(err.stat^2 + err.msys^2), err.a, 0.0 ))
      MDX <- rbind(MDX,phys$val*cbind( err.stat, err.msys, err.a))
      errband <- rbind(errband,data.frame( xleft=(1-err.phys), xright=(1+err.phys) ) ) 
      labels <- c(labels,ext[[i]]$texlabel[nr])
    }
  }
  
  for(name in values) {
    for(i in 1:length(ext)){
      if(ext[[i]]$name != name) next
      # the relevant extrapolated value is always the last one
      nr <- nrow(ext[[i]])
      phys <- physvalues[ext[[i]]$name[1]==physvalues$name,]
     
      err.stat <- hbarc*ext[[i]]$dval[nr]/phys$val/a$val
      err.sys <- hbarc*ext.sys[[i]]$serr[nr]/phys$val/a$val
      err.msys <- hbarc*ext.sys[[i]]$mserr[nr]/phys$val/a$val
      err.a <- hbarc*a$dval*ext[[i]]$val[nr]/phys$val/(a$val^2)
      err.phys <- hbarc*phys$dval*ext[[i]]$val[nr]/(phys$val^2)/a$val
      
      x <- c(x,hbarc*ext[[i]]$val[nr]/a$val/phys$val)
      X <- c(X,hbarc*ext[[i]]$val[nr]/a$val)
      dx <- rbind(dx,cbind( 0.0, sqrt(err.stat^2 + err.sys^2), err.a, 0.0 ))
      DX <- rbind(DX,phys$val*cbind( err.stat, err.sys, err.a ))
      mdx <- rbind(mdx,cbind( 0.0, sqrt(err.stat^2 + err.msys^2), err.a, 0.0 ))
      MDX <- rbind(MDX,phys$val*cbind( err.stat, err.msys, err.a))
      val.idx <- which(values==ext[[i]]$name[1])
      errband <- rbind(errband,data.frame(xleft=(1-err.phys),xright=(1+err.phys) ) )
      labels <- c(labels,value.labels[val.idx])
       
    } 
  }
  # add the nucleon
  err.stat <- hbarc*mN$dval/mN.pdg$val/a$val
  err.sys <- 0.0
  err.msys <- 0.0
  err.a <- hbarc*mN$val*a$dval/mN.pdg$val/(a$val^2)
  err.phys <- hbarc*mN.pdg$dval*mN$val/phys$val^2/a$val

  x <- c(x,hbarc*mN$val/mN.pdg$val/a$val)
  X <- c(X,hbarc*mN$val/a$val)
  dx <- rbind(dx,cbind(err.stat,err.sys,err.a,0.0))
  DX <- rbind(DX,mN.pdg$val*cbind(err.stat,0.0,err.a))
  mdx <- rbind(mdx,cbind(err.stat,err.msys,err.a,0.0))
  MDX <- rbind(MDX,mN.pdg$val*cbind(err.stat,0.0,err.a))
  errband <- rbind(errband, data.frame(xleft=(1-err.phys),xright=(1+err.phys) ) )
  labels <- c(labels,"$M_N$")
   
  list(labels=labels,x=x,dx=dx,mdx=mdx,errband=errband,summary=data.frame(cbind(labels,X=X,DX=DX,MDX=MDX)) )
}

