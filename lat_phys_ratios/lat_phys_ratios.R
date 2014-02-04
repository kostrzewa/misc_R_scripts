compute_ratio <- function(dividend,divisor,name,debug=FALSE) {
  ratio <- list( value=dividend$value / divisor$value, 
                 error=sqrt( (dividend$error/divisor$value)^2 + (divisor$error*dividend$value/divisor$value^2)^2 ), 
                 name=name )
  if(debug) {
    print(sprintf("compute_ratio: %s",as.character(name)))
    print(ratio)
  }
  return(ratio)
}

# extract an observable from table produced by read.table from the format agree upon
# ....

extract_observable <- function(datatable,name,debug=FALSE) {
  if(debug) {
    print(sprintf("Extract_observable %s",name))
    print(datatable[ datatable$V1==name, ])
  }

  if( nrow(datatable[ datatable$V1==name, ]) == 0 ) {
    stop(paste("aborting: observable could not be found in data file:",name))
  }

  lat <- list(value=datatable[ datatable$V1==name, ][2], error=datatable[ datatable$V1==name, ][3])

  # if there exists a FLAG value for the physical value of an observable, we use this
  # otherwise we fall back to the PDG value
  offset <- 0
  if( !is.na(datatable[ datatable$V1==name, ][6]) ) {
    offset <- 2
  }
  phys <- list(value=datatable[ datatable$V1==name, ][(4+offset)],error=datatable[ datatable$V1==name, ][(5+offset)])
  return( list(lat=lat,phys=phys) )
}

lat_phys_ratios <- function(filename,ratios,xlim=c(0.8,1.2),labelpos=1.15,debug=FALSE) {
  par(family="Times")

  datapts <- read.table(filename)
  if(debug) {
    print(datapts)
  }

  if(debug) {
    print(ratios)
  }

  num_observables <- nrow(ratios)
# m_pi / f_pi
  
  rep <- FALSE
  for(i in seq(1,num_observables)) {
    dividend_name <- as.character(ratios[i,]$dividend_name)
    divisor_name <- as.character(ratios[i,]$divisor_name)
    
    if(debug) {
      print(sprintf("Observable %d %s / %s" ,i,dividend_name,divisor_name))
    }

    dividend <- extract_observable(datapts,dividend_name,debug)
    divisor <- extract_observable(datapts,divisor_name,debug)
   
    name <- sprintf("%s / %s (lat)",dividend_name,divisor_name) 
    ratio_lat <- compute_ratio(dividend$lat,divisor$lat,name=name,debug)
    name <- sprintf("%s / %s (phys)",dividend_name,divisor_name) 
    ratio_phys <- compute_ratio(dividend$phys,divisor$phys,name=name,debug)

    ratio_lat_phys <- compute_ratio(ratio_lat,ratio_phys,name=as.character(ratios[i,]$ratio_name),debug)

    x <- as.numeric( ratio_lat_phys$value )
    y <- 1+num_observables-i
    dx <- as.numeric( ratio_lat_phys$error )
    name <- as.character(ratio_lat_phys$name)
    
    if(i > 1) {
      rep <- TRUE
    } 
    
    if( i %% 2 == 0 ) {
      col <- rgb(0.8,0.8,0.8,0.4)
      rect(xlim[1],y-0.5,xlim[2],y+0.5,col=col,border=NA)
    }
    
    plotwitherror(x=x, y=y, dx=dx,
      yaxt='n',
      ylab="",xlab=expression(Q[lat]/Q[phys]),
      rep=rep,xlim=xlim,ylim=c(0.5,num_observables+0.5),lwd=2)

    
    par(xpd=NA)
    text(y=y,x=labelpos,label=name)
    par(xpd=FALSE)
  }
  abline(v=1,lwd=1,lty=2)
  # don't really like the look of this, even though it makes it more consistent
  #rect(xlim[1],1-0.5,xlim[2],num_observables+0.5,col=NA)
}


plot_lat_phys_ratios <- function(filename,debug=FALSE) {
  #pdf("light_strange_ratios.pdf",width=2.7,height=5)
  require(tikzDevice)
  tikz('light_strange_ratios.tex', standAlone = TRUE, width=2.36, height=4)
  ratios <- data.frame(dividend_name="m_pi",divisor_name="f_pi",ratio_name="$m_\\pi/f_\\pi$")
  
  #ratios <- rbind(ratios,data.frame(dividend_name="m_K(a)",divisor_name="m_pi",ratio_name="m[K]/m[Pi] (a)"))
  #ratios <- rbind(ratios,data.frame(dividend_name="f_K(a)",divisor_name="f_pi",ratio_name="f[K]/f[Pi] (a)"))
 
  ratios <- rbind(ratios,data.frame(dividend_name="m_K(b)",divisor_name="m_pi",ratio_name="$m_K/m_\\pi$"))
  ratios <- rbind(ratios,data.frame(dividend_name="f_K(b)",divisor_name="f_pi",ratio_name="$f_K/f_\\pi$"))

  #ratios <- rbind(ratios,data.frame(dividend_name="m_N",divisor_name="m_pi",ratio_name="m[N]/m[Pi]"))
  #ratios <- rbind(ratios,data.frame(dividend_name="g_A",divisor_name="unity",ratio_name="g[A]"))
  #ratios <- rbind(ratios,data.frame(dividend_name="<x>",divisor_name="unity",ratio_name="x[ud]"))

  lat_phys_ratios(filename,ratios,xlim=c(0.995,1.032),labelpos=1.026,debug=debug)

  dev.off()
  tools::texi2dvi("light_strange_ratios.tex",pdf=T)
  system("pdfcrop light_strange_ratios.pdf light_strange_ratios.pdf")

  #pdf("heavylight_charm_ratios.pdf",height=5,width=5)
  tikz('heavylight_charm_ratios.tex', standAlone = TRUE, width=4.1, height=4)
  #ratios <- data.frame(dividend_name="m_D(a)",divisor_name="m_pi",ratio_name="m[D]/m[Pi] (a)")
  #ratios <- rbind(ratios,data.frame(dividend_name="f_D(a)",divisor_name="f_pi",ratio_name="f[D]/f[Pi] (a)"))

  ratios <- data.frame(dividend_name="m_D(b)",divisor_name="f_pi",ratio_name="$m_D/f_\\pi$")
  ratios <- rbind(ratios,data.frame(dividend_name="f_D(b)",divisor_name="f_pi",ratio_name="$f_D/f_\\pi$"))

  #ratios <- rbind(ratios,data.frame(dividend_name="m_Ds(a,a)",divisor_name="m_pi",ratio_name="m[Ds]/m[Pi] (a)"))
  #ratios <- rbind(ratios,data.frame(dividend_name="f_Ds(a,a)",divisor_name="f_pi",ratio_name="f[Ds]/f[Pi] (a)"))
  #ratios <- rbind(ratios,data.frame(dividend_name="f_Ds(a,a)",divisor_name="f_D(a)",ratio_name="f[Ds]/f[D] (a,a)"))
   
  ratios <- rbind(ratios,data.frame(dividend_name="m_Ds(b,b)",divisor_name="f_pi",ratio_name="$m_{D_s}/f_\\pi$")) 
  ratios <- rbind(ratios,data.frame(dividend_name="f_Ds(b,b)",divisor_name="f_pi",ratio_name="$f_{D_s}/f_\\pi$")) 
  ratios <- rbind(ratios,data.frame(dividend_name="f_Ds(b,b)",divisor_name="f_D(b)",ratio_name="$f_{D_s}/f_D$"))                     
  
  # ratios for the nucleon sector  
  #ratios <- rbind(ratios,data.frame(dividend_name="m_N",divisor_name="m_pi",ratio_name="m[N]/m[Pi]"))
  #ratios <- rbind(ratios,data.frame(dividend_name="m_N",divisor_name="f_pi",ratio_name="m[N]/f[Pi]"))
  #ratios <- rbind(ratios,data.frame(dividend_name="g_A",divisor_name="unity",ratio_name="g[A]"))
  #ratios <- rbind(ratios,data.frame(dividend_name="<x>",divisor_name="unity",ratio_name="x[ud]"))

  lat_phys_ratios(filename,ratios,xlim=c(0.9,1.13),labelpos=1.11,debug=debug)

  dev.off()
  tools::texi2dvi("heavylight_charm_ratios.tex",pdf=T)
  system("pdfcrop heavylight_charm_ratios.pdf heavylight_charm_ratios.pdf")
}

