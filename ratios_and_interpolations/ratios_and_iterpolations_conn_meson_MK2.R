library('propagate')

source("/home/bartek/code/R/misc_R_scripts/fit_extrapolate_solve_MK2.R")
source("~/code/R/misc_R_scripts/hadron_obs.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/utils.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/utils_MK2.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/mK_ov_fK.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/match_mu_s_mu_c.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/match_mu.R") 


# TODO: move entire toolset from m1/m2 mass specification to 

ratios_and_iterpolations_conn_meson_mk2 <- function(debug=F,recompute=T,loadraw=T,overview=T) {
  # certain functionality relies on stuff being strings
  options(stringsAsFactors = FALSE)
  # masses to be used in this analysis
  strange_masses <- c(0.0238,0.0245,0.0252,0.0259)
  charm_masses <- c(0.2822,0.294,0.3058,0.3176)
  light_masses <- c(0.0009)
  m.sea <- light_masses
  
  phys_ratios <- read.csv("phys_ratios.csv",as.is=c("type","name"))
  
  pheno.col <- rgb(green=1.0,red=0.0,blue=0.0,alpha=0.3)
  pheno.pch <- 15
  
  # combinations of these masses
  mass_comb <- list( ll=data.frame( m1=light_masses, m2=light_masses ),
                     ls=expand.grid( m1=light_masses, m2=strange_masses),
                     lc=expand.grid( m1=light_masses, m2=charm_masses),
                     sc=expand.grid( m1=strange_masses, m2=charm_masses) )               
  
  # File (and object) names of the data to be loaded
  names <- list(ll_c=sprintf("ll_c.m%g.m%g.matrixfit",mass_comb$ll$m1,mass_comb$ll$m2),
                #ll_na=sprintf("lln_u_%g-d_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                #ll_nb=sprintf("lln_d_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ls_c=sprintf("ls_c.m%g.m%g.matrixfit",mass_comb$ls$m1,mass_comb$ls$m2),
                lc_c=sprintf("lc_c.m%g.m%g.matrixfit",mass_comb$lc$m1,mass_comb$lc$m2),
                sc_c=sprintf("sc_c.m%g.m%g.matrixfit",mass_comb$sc$m1,mass_comb$sc$m2) )

  # single quantities
  quants <- NULL
  quants[["m_pi"]] <- list(name="m_pi", texlabel="$m_\\pi$", m1=light_masses, m2=light_masses, datanames=names$ll_c, type="mps" )
  quants[["m_K"]] <- list(name="m_K", texlabel="$m_K$", m1=light_masses, m2=strange_masses, datanames=names$ls_c, type="mps" )
  quants[["m_D"]] <- list(name="m_D", texlabel="$m_D$", m1=light_masses, m2=charm_masses, datanames=names$lc_c, type="mps" )
  quants[["m_Ds"]] <- list(name="m_Ds", texlabel="$m_{D_s}$", m1=strange_masses, m2=charm_masses, datanames=names$sc_c, type="mps" )
  quants[["f_pi"]] <- list(name="f_pi", texlabel="$f_\\pi$", m1=light_masses, m2=light_masses, datanames=names$ll_c, type="fps" )
  quants[["f_K"]] <- list(name="f_K", texlabel="$f_K$", m1=light_masses, m2=strange_masses, datanames=names$ls_c, type="fps" )
  quants[["f_D"]] <- list(name="f_D", texlabel="$f_D$", m1=light_masses, m2=charm_masses, datanames=names$lc_c, type="fps" )
  quants[["f_Ds"]] <- list(name="f_Ds", texlabel="$f_{D_s}$", m1=strange_masses, m2=charm_masses, datanames=names$sc_c, type="fps" )

  # the different ratios that we would like to compute
  ratios <- NULL
  ratios[["m_pi_ov_f_pi"]] <- list( name="m_pi_ov_f_pi", texlabel="$m_\\pi/f_\\pi$", dividend=quants[["m_pi"]], divisor=quants[["f_pi"]])
  ratios[["m_K_ov_f_K"]] <- list( name="m_K_ov_f_K", texlabel="$m_K/f_K$", dividend=quants[["m_K"]], divisor=quants[["f_K"]])
  ratios[["m_D_ov_f_D"]] <- list( name="m_D_ov_f_D", texlabel="$m_D/f_D$", dividend=quants[["m_D"]], divisor=quants[["f_D"]])
  ratios[["m_Ds_ov_f_Ds"]] <- list( name="m_Ds_ov_f_Ds", texlabel="$m_{D_s}/f_{D_s}$", dividend=quants[["m_Ds"]], divisor=quants[["f_Ds"]])

  ratios[["m_K_ov_f_pi"]] <- list( name="m_K_ov_f_pi", texlabel="$m_K/f_\\pi$", dividend=quants[["m_K"]], divisor=quants[["f_pi"]])
  ratios[["m_D_ov_f_pi"]] <- list( name="m_D_ov_f_pi", texlabel="$m_D/f_\\pi$", dividend=quants[["m_D"]], divisor=quants[["f_pi"]])
  ratios[["m_Ds_ov_f_pi"]] <- list( name="m_Ds_ov_f_pi", texlabel="$m_{D_s}/f_\\pi$", dividend=quants[["m_Ds"]], divisor=quants[["f_pi"]])

  ratios[["m_Ds_ov_m_D"]] <- list( name="m_Ds_ov_m_D", texlabel="$m_{D_s}/m_D$", dividend=quants[["m_Ds"]], divisor=quants[["m_D"]])
  ratios[["m_Ds_ov_m_K"]] <- list( name="m_Ds_ov_m_K", texlabel="$m_{D_s}/m_K$", dividend=quants[["m_Ds"]], divisor=quants[["m_K"]])
  ratios[["m_Ds_ov_m_pi"]] <- list( name="m_Ds_ov_m_pi", texlabel="$m_{D_s}/m_\\pi$", dividend=quants[["m_Ds"]], divisor=quants[["m_pi"]])

  ratios[["m_D_ov_m_K"]] <- list( name="m_D_ov_m_K", texlabel="$m_{D}/m_K$", dividend=quants[["m_D"]], divisor=quants[["m_K"]])
  ratios[["m_D_ov_m_pi"]] <- list( name="m_D_ov_m_pi", texlabel="$m_{D}/m_\\pi$", dividend=quants[["m_D"]], divisor=quants[["m_pi"]])

  ratios[["m_K_ov_m_pi"]] <- list( name="m_K_ov_m_pi", texlabel="$m_K/m_\\pi$", dividend=quants[["m_K"]], divisor=quants[["m_pi"]])

  ratios[["f_K_ov_f_pi"]] <- list( name="f_K_ov_f_pi", texlabel="$f_K/f_\\pi$", dividend=quants[["f_K"]], divisor=quants[["f_pi"]])
  ratios[["f_D_ov_f_pi"]] <- list( name="f_D_ov_f_pi", texlabel="$f_D/f_\\pi$", dividend=quants[["f_D"]], divisor=quants[["f_pi"]])
  ratios[["f_Ds_ov_f_pi"]] <- list( name="f_Ds_ov_f_pi", texlabel="$f_{D_s}/f_\\pi$", dividend=quants[["f_Ds"]], divisor=quants[["f_pi"]])

  ratios[["f_D_ov_f_K"]] <- list( name="f_D_ov_f_K", texlabel="$f_D/f_K$", dividend=quants[["f_D"]], divisor=quants[["f_K"]])
  ratios[["f_Ds_ov_f_K"]] <- list( name="f_Ds_ov_f_K", texlabel="$f_{D_s}/f_K$", dividend=quants[["f_Ds"]], divisor=quants[["f_K"]])

  ratios[["f_Ds_ov_f_D"]] <- list( name="f_Ds_ov_f_D", texlabel="$f_{D_s}/f_D$", dividend=quants[["f_Ds"]], divisor=quants[["f_D"]])

  if(loadraw || recompute) {
    # read all the data files
    for( n in names ) {
      for( file in n ) {
        file <- sprintf("%s.Rdata",file)
        if(debug) cat("Loading" ,file,"\n")
        load(file)
      }
    }
  }

  if(recompute) {
  
    hadron_obs <- construct.empty.hadron_obs()
    
    for( qty in quants ) {
      for( objname in qty$datanames ) {
        dat <- get.obs.tsboot_mk2(obj=get(objname), m.sea=m.sea, type=qty$type)
        hadron_obs[[length(hadron_obs)+1]] <- compute.quant_mk2( name=qty$name, texlabel=qty$texlabel, dat=dat )
      }
    }
    
    for( r in ratios ) {
      # there is a special case if the data of dividend and divisor are exactly the same,
      # we don't need x*x cases but only x! This would be the case for a meson mass
      # and its decay constant, for example
      if(debug) cat("Computing ratio", r$name, "\n")
      if( all(r$dividend$datanames == r$divisor$datanames) ) {
        for( objname in r$dividend$datanames ) {
          dividend <- get.obs.tsboot_mk2(obj=get(objname), type=r$dividend$type, m.sea=m.sea )
          divisor <- get.obs.tsboot_mk2(obj=get(objname), type=r$divisor$type, m.sea=m.sea )
          hadron_obs[[length(hadron_obs)+1]] <- compute.ratio_mk2(name=r$name,texlabel=r$texlabel,dividend=dividend,divisor=divisor)
        }
      } else {

        # when any pair of mass vectors from the two observables are the same
        # we need to ensure that "off-diagonal" elements, where say the strange
        # mass for the divisor and dividend are different, are skipped!
        # this would be the case, for instance for a ratio involving the kaon
        # and the D_s meson
         m1m1 <- FALSE
         m1m2 <- FALSE
         m2m1 <- FALSE
         m2m2 <- FALSE
         if( all(r$dividend$m1 == r$divisor$m1) ) {
           m1m1 <- TRUE
         } else if( all(r$dividend$m1 == r$divisor$m2) ) {
          m1m2 <- TRUE
         } else if( all(r$dividend$m2 == r$divisor$m1) ) {
           m2m1 <- TRUE
         } else if( all(r$dividend$m2 == r$divisor$m2) ) {
           m2m2 <- TRUE
         }

        for( dividend.objname in r$dividend$datanames ) {
          for( divisor.objname in r$divisor$datanames ) {
            dividend <- get.obs.tsboot_mk2(obj=get(dividend.objname), type=r$dividend$type, m.sea=m.sea )
            divisor <- get.obs.tsboot_mk2(obj=get(divisor.objname), type=r$divisor$type, m.sea=m.sea )
            
            # note the second condition: if no two mass vectors are the same
            if( any(as.vector(outer(divisor$m.val,dividend$m.val,'=='))) ||
               !(m1m1 || m1m2 || m2m1 || m2m2) )  {
              hadron_obs[[length(hadron_obs)+1]] <- compute.ratio_mk2(name=r$name,texlabel=r$texlabel,dividend=dividend,divisor=divisor)
            }
          }
        }
      }
    }
    
#    # ratios from expressions
#    cat("Computing ratio f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD\n")
#    for( i in 1:length(mass_comb$sc[,1]) ) {
#      fDs_indices <- which(hadron_obs$val.tsboot$name == "f_Ds" & 
#                          hadron_obs$val.tsboot$m11 == mass_comb$sc[i,]$m1 & 
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      mDs_indices <- which(hadron_obs$val.tsboot$name == "m_Ds" & 
#                          hadron_obs$val.tsboot$m11 == mass_comb$sc[i,]$m1 & 
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      fD_indices <- which(hadron_obs$val.tsboot$name == "f_D" & 
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      mD_indices <- which(hadron_obs$val.tsboot$name == "m_D" & 
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#                          
#      dividend <- compute.expression(expr=expression(a*sqrt(b)),
#                    envir=data.frame( a=hadron_obs$val.tsboot$val[fDs_indices], b=hadron_obs$val.tsboot$val[mDs_indices] ),
#                    m1=mass_comb$sc[i,]$m1, m2=mass_comb$sc[i,]$m2 )
#      divisor <- compute.expression(expr=expression(a*sqrt(b)),
#                    envir=data.frame( a=hadron_obs$val.tsboot$val[fD_indices], b=hadron_obs$val.tsboot$val[mD_indices] ),
#                    m1=0.0009, m2=mass_comb$sc[i,]$m2 )
#                    
#      hadron_obs <- compute.ratio(name="f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD",texlabel="$f_{D_s}\\sqrt{m_{D_s}}/f_D\\sqrt{m_D}$",
#                               dividend=dividend, divisor=divisor, res=hadron_obs )
#    }

#    cat("Computing ratio f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd\n")
#    for( i in 1:length(mass_comb$sc[,1]) ) {
#      fDs_indices <- which(hadron_obs$val.tsboot$name == "f_Ds" & 
#                          hadron_obs$val.tsboot$m11 == mass_comb$sc[i,]$m1 & 
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      mDs_indices <- which(hadron_obs$val.tsboot$name == "m_Ds" & 
#                          hadron_obs$val.tsboot$m11 == mass_comb$sc[i,]$m1 & 
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      fD_indices <- which(hadron_obs$val.tsboot$name == "f_D" & 
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#      mD_indices <- which(hadron_obs$val.tsboot$name == "m_D" & 
#                          hadron_obs$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
#                          
#      dividend <- compute.expression(expr=expression(a*b^2),
#                    envir=data.frame( a=hadron_obs$val.tsboot$val[fDs_indices], b=hadron_obs$val.tsboot$val[mDs_indices] ),
#                    m1=mass_comb$sc[i,]$m1, m2=mass_comb$sc[i,]$m2 )
#      divisor <- compute.expression(expr=expression(a*b^2),
#                    envir=data.frame( a=hadron_obs$val.tsboot$val[fD_indices], b=hadron_obs$val.tsboot$val[mD_indices] ),
#                    m1=0.0009, m2=mass_comb$sc[i,]$m2 )
#                    
#      hadron_obs <- compute.ratio(name="f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd",texlabel="$f_{D_s}m_{D_s}^2/f_D{m_D}^2$",
#                               dividend=dividend, divisor=divisor, res=hadron_obs )
#    }
    
    save(hadron_obs,file="hadron_obs.Rdata")
  } else {
    load("hadron_obs.Rdata")
  } # if(recompute)
  
#  print(select.hadron_obs(hadron_obs,by='name',filter='f_Ds'))
#  stop("Computed all quantities and ratios!")

  
  if(overview) {
    # overview plots of data
    require(tikzDevice)
    filebase <- "data_overview"
    texfile <- sprintf("%s.tex",filebase) 
    pdffile <- sprintf("%s.pdf",filebase)
    tikz(texfile, standAlone = TRUE, width=5, height=5)
    for(name in c(names(quants),names(ratios))) {
      if(debug) print(name)
      #ts.indices <- which( hadron_obs$res.tsboot$name == name )
      indices <- which( hadron_obs$val$name == name )
      obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
      print(length(obs))
      obs.df <- extract.for.plot(hadron_obs=obs,x.name='m.val',x.idx=length(obs[[1]]$m.val))
      plotwitherror(y=obs.df$y, dy=obs.df$dy, x=obs.df$x,
                    main=obs[[1]]$texlabel, xlab="$\\mu$",ylab=obs[[1]]$texlabel)
    }
    
#    if(debug) print("f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD")
#    indices <- which( hadron_obs$val$name == "f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD" )
#    plotwitherror(y=hadron_obs$val$val[indices], dy=hadron_obs$val$dval[indices], x=hadron_obs$val$m12[indices],
#                  main=hadron_obs$val$texlabel[indices[1]], xlab="$\\mu$",ylab=hadron_obs$val$texlabel[indices[1]])

#    if(debug) print("f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd")
#    indices <- which( hadron_obs$val$name == "f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd" )
#    plotwitherror(y=hadron_obs$val$val[indices], dy=hadron_obs$val$dval[indices], x=hadron_obs$val$m12[indices],
#                  main=hadron_obs$val$texlabel[indices[1]], xlab="$\\mu$",ylab=hadron_obs$val$texlabel[indices[1]])
    
    dev.off()
    tools::texi2dvi(texfile,pdf=T)
    command <- sprintf("pdfcrop %s %s",pdffile,pdffile)
    system(command)
  } # if(overview)
  
  # mu_s from FLAG ratio
  mu_s <- data.frame( val=c(0.0009*27.46), dval=c(0.44*0.0009), name="FLAG")
  
  extrapolations <- data.frame(name=c(),val=c(),dval=c(),x=c(),dx=c())
  
  #print(phys_ratios)  
  #readline("press key")
  
  # m_K_ov_f_K
  name <- "m_K_ov_f_K"
  pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.0215, y=3.2 )
  # asemble the data into the correct format
  obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
  pred.idx <- list(m.val=c(2),m.sea=vector())
  dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
  fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
  pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval)
  fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)
  
  extrapolations <- rbind(extrapolations, 
                          data.frame(name=name,
                                     val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                     dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ),
                                     x=mu_s[1,]$val, dx=mu_s[1,]$dval
                                     )
                         )
  
  fes.solve <- fes_solve(fesfit=fes.fit,unknown='x1',known=c(),
                         y=phys_ratios[phys_ratios$name == name,]$val,
                         dy=phys_ratios[phys_ratios$name == name,]$dval )
  
  solution <- data.frame(val=mean(fes.solve[,1]), dval=sd(fes.solve[,1]), name="m_K_ov_f_K")
                         
  df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
  
  #print(df)
  plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[extrapolations$name==name,],solutions=solution,
                  xlab="$\\mu_s$",ylab=obs[[1]]$texlabel)
  
#   stop()
    
  mu_s <- rbind( mu_s, solution )
  
  cols <- c('black','red','forestgreen')
  syms <- c(1,15:(15+length(mu_s$val)-1))
  
  legend.mu_s <- list( labels=c("Data","$\\mu_s$ from FLAG ratio","$\\mu_s$ from $m_K/f_K$","$\\mu_s$ from $m_K/m_\\pi$"),
                       pch=c(syms,18), col=c(cols,'blue') )
  legend.mu_c <- list( labels=c("Data","$\\mu_c$ from FLAG$\\cdot$HPQCD ratios",
                                "$\\mu_c=$ HPQCD ratio $\\cdot\\mu_s$ from $m_K/f_K$","$\\mu_c$ from $m_D/m_\\pi$"),
                        pch=c(syms,18), col=c(cols,'blue') )
  
  name <- "m_K_ov_m_pi"
  pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  # asemble the data into the correct format
  obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
  pred.idx <- list(m.val=c(2),m.sea=vector())
  dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
  fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
  pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval)
  
  fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)
    
  extrapolations <- rbind(extrapolations, 
                          data.frame(name=name,
                                     val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                     dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ),
                                     x=mu_s$val, dx=mu_s$dval
                                     )
                         )  
  
  fes.solve <- fes_solve(fesfit=fes.fit,unknown='x1',known=c(),
                         y=phys_ratios[phys_ratios$name == name,]$val,
                         dy=phys_ratios[phys_ratios$name == name,]$dval )

  df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
  plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[extrapolations$name==name,],solutions=solution,
                  xlab="$\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_s)                       

  stop()                
                  
  mu_s <- rbind( mu_s, data.frame(val=mean(fes.solve[,1]), dval=sd(fes.solve[,1]), name="m_K_ov_m_pi" ) )
  
  lg.coords <- data.frame( x=0.021, y=3.77 )
  # plot
  
#   print(extrapolations)
#   print(mu_s)
#   readline("press key")

  # mu_c from HPQCD ratio, this will now contains three values of mu_c to which we will add a fourth by solving
  mu_c <- data.frame( val=mu_s$val*11.85, 
                      dval= c( 0.0009*sqrt( (11.85*0.44)^2 + (27.46*0.16)^2 ), 
              sqrt( (mu_s$dval[2:3]*11.85)^2 + (mu_s$val[2:3]*0.16)^2 ) ),
                      name=sprintf("%s/%s",mu_s$name,"HPQCD"))
  
  name <- "m_D_ov_m_pi"
#   lg.coords <- data.frame( x=0.27, y=15.2 )
  obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
  pred.idx <- list(m.val=c(2),m.sea=vector())
  dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
  fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
  pred <- data.frame(x1=mu_c$val,dx1=mu_c$dval)
  fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

  extrapolations <- rbind(extrapolations, 
                          data.frame(name=name,
                                     val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                     dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ),
                                     x=mu_c[1:3,]$val, dx=mu_c[1:3,]$dval
                                     )
                         )  

  fes.solve <- fes_solve(fesfit=fes.fit,unknown='x1',known=c(),
                         y=phys_ratios[phys_ratios$name == name,]$val,
                         dy=phys_ratios[phys_ratios$name == name,]$dval )
  mu_c <- rbind( mu_c, data.frame(val=mean(fes.solve[,1]), dval=sd(fes.solve[,1]), name="m_D_ov_m_pi") )

  readline("press key")
  cat("Extrapolations and quark masses\n")
  print(extrapolations)
  print(mu_s)
  print(mu_c)
  readline("press key")
  
#   pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
#   mu_c_from_m_D <- match_mu.1D(name=name,alldat=hadron_obs,masses=charm_masses,
#               lg=legend.mu_c,lg.coords=lg.coords,
#               pheno=pheno,mu=mu_c, xlab="$a\\mu_c$", solve=T,
#               xlim=c(0.27,0.36),ylim=c(12.8,15.2))

  cols <- c(cols,'blue')
  syms <- c(1,15:(15+length(mu_c$val)-1))
  
  legend.mu_s <- list( labels=c("Data","$\\mu_s$ from FLAG ratio","$\\mu_s$ from $m_K/f_K$","$\\mu_s$ from $m_K/m_\\pi$"),
                       pch=syms, col=cols )
  legend.mu_c <- list( labels=c("Data","$\\mu_c$ from FLAG$\\cdot$HPQCD ratios",
                                "$\\mu_c=$ HPQCD ratio $\\cdot\\mu_s$ from $m_K/f_K$","$\\mu_c$ from $m_D/m_\\pi$"),
                        pch=syms, col=cols )
  legend.mu_sc <- list( labels=c("Data", "$\\mu_s$ and $\\mu_c$ from FLAG/HPQCD ratios",
                                 "$\\mu_c=$ HPQCD ratio $\\cdot\\mu_s$ from $m_K/f_K$",
                                 "$\\mu_c$ from $m_D/m_\\pi$, $\\mu_s$ from $m_K/m_\\pi$"), 
                        pch=syms, col=cols )

  pred <- list()
  
  name <- "m_D_ov_f_D"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.3, y=8.25 )
  pred[[name]] <- match_mu.1D(name=name,alldat=hadron_obs,masses=charm_masses,
              pheno=pheno,mu=mu_c,
              lg=legend.mu_c,lg.coords=lg.coords,xlab="$a\\mu_c$", solve=F,
              xlim=c(0.28,0.36))
              
  name <- "m_Ds_ov_f_Ds"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.0217, y=8.2 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=hadron_obs,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c, lg=legend.mu_sc, lg.coords=lg.coords,ylim=c(7.2,8.2))              
  
  name <- "m_Ds_ov_m_K"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.0209, y=4.31 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=hadron_obs,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c, lg=legend.mu_sc, lg.coords=lg.coords,ylim=c(3.7,4.31),xlim=c(0.021,0.027))  
  
  name <- "m_Ds_ov_m_D"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.0209, y=1.07 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=hadron_obs,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c, lg=legend.mu_sc, lg.coords=lg.coords,
                     ylim=c(1.045,1.07),xlim=c(0.021,0.027))

  name <- "f_K_ov_f_pi"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.023, y=1.22 )
  pred[[name]] <- match_mu.1D(name=name,alldat=hadron_obs,masses=strange_masses,
              pheno=pheno,mu=mu_s,
              lg=legend.mu_s,lg.coords=lg.coords,xlab="$a\\mu_s$",solve=F,
              ylim=c(1.19,1.22),xlim=c(0.023,0.027))

  name <- "f_D_ov_f_pi"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.255, y=1.8 )
  pred[[name]] <- match_mu.1D(name=name,alldat=hadron_obs,masses=charm_masses,
              pheno=pheno,mu=mu_c,
              lg=legend.mu_c,lg.coords=lg.coords,xlab="$a\\mu_c$",
              ylim=c(1.5,1.8))                         

  name <- "m_D_ov_m_K"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.021, y=3.5 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=hadron_obs,mass_comb=mass_comb,pheno=pheno,
                     mu_s,mu_c,
                     m11=0.0009,m12=charm_masses,
                     m21=0.0009,m22=strange_masses,
                     lg=legend.mu_sc, lg.coords=lg.coords,
                     xval="m22", debug=T,
                     xlim=c(0.021,0.027),
                     ylim=c(3.3,3.9) )
  
  name <- "m_Ds_ov_m_pi"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.023, y=15.4 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=hadron_obs,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c, lg=legend.mu_sc, lg.coords=lg.coords,
                     xlim=c(0.023,0.027),
                     ylim=c(13.8,15.4))

  name <- "f_Ds_ov_f_pi"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.023, y=2.04 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=hadron_obs,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c,lg=legend.mu_sc, lg.coords=lg.coords,
                     xlim=c(0.023,0.027),ylim=c(1.94,2.04))

  name <- "f_Ds_ov_f_D"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.023, y=1.3 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=hadron_obs,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c,lg=legend.mu_sc,lg.coords=lg.coords,
                     xlim=c(0.023,0.027),
                     ylim=c(1.16,1.3))

  name <- "m_K"
  pred[[name]] <- match_mu.1D(name=name,alldat=hadron_obs,masses=strange_masses,
              mu=mu_s)
                     
  name <- "f_K"
  pred[[name]] <- match_mu.1D(name=name,alldat=hadron_obs,masses=strange_masses,
              mu=mu_s)
                     
  name <- "m_D"
  pred[[name]] <- match_mu.1D(name=name,alldat=hadron_obs,masses=charm_masses,
              mu=mu_c)
                     
  name <- "f_D"
  pred[[name]] <- match_mu.1D(name=name,alldat=hadron_obs,masses=charm_masses,
              mu=mu_c)
                                   
  name <- "m_Ds"
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=hadron_obs,mass_comb=mass_comb,
                  mu_s=mu_s, mu_c=mu_c)
  
  name <- "f_Ds"
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=hadron_obs,mass_comb=mass_comb,
                  mu_s=mu_s, mu_c=mu_c)
                  
  print(pred)
  
  write.csv(pred,file="preds.csv")
  
  options(stringsAsFactors = TRUE)
}
