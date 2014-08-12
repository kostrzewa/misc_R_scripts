library('propagate')

source("/home/bartek/code/R/misc_R_scripts/fit_extrapolate_solve.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/utils.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/mK_ov_fK.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/match_mu_s_mu_c.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/match_mu.R")

ratios_and_iterpolations_conn_meson <- function(debug=F,recompute=T,loadraw=T,overview=T) {
  # certain functionality relies on stuff being strings
  options(stringsAsFactors = FALSE)
  # masses to be used in this analysis
  strange_masses <- c(0.0238,0.0245,0.0252,0.0259)
  charm_masses <- c(0.2822,0.294,0.3058,0.3176)
  light_masses <- c(0.0009)
  
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
  
    results <- NULL
    
    for( qty in quants ) {
      for( objname in qty$datanames ) {
        dat <- get.obs.tsboot(obj=get(objname), type=qty$type)
        results <- compute.quant( name=qty$name, texlabel=qty$texlabel, dat=dat, res=results )
      }
    }
    
    for( r in ratios ) {
      # there is a special case if the data of dividend and divisor are exactly the same,
      # we don't need x*x cases but only x!
      if(debug) cat("Computing ratio", r$name, "\n")
      if( all(r$dividend$datanames == r$divisor$datanames) ) {
        for( objname in r$dividend$datanames ) {
          dividend <- get.obs.tsboot(obj=get(objname), type=r$dividend$type )
          divisor <- get.obs.tsboot(obj=get(objname), type=r$divisor$type )
          results <- compute.ratio(name=r$name,texlabel=r$texlabel,dividend=dividend,divisor=divisor,res=results)
        }
      } else {

        # when any pair of mass vectors from the two observables are the same
        # we need to ensure that "off-diagonal" elements, where say the strange
        # mass for the divisor and dividend are different, are skipped!
        # this would be the case, for instance for a ratio involving the kaon
        # and the D_s meson or the kaon mass and its decay constant
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
            dividend <- get.obs.tsboot(obj=get(dividend.objname), type=r$dividend$type )
            divisor <- get.obs.tsboot(obj=get(divisor.objname), type=r$divisor$type )
            if( ( m1m1 && dividend$m1 == divisor$m1 ) ||
                ( m1m2 && dividend$m1 == divisor$m2 ) ||
                ( m2m1 && dividend$m2 == divisor$m1 ) ||
                ( m2m2 && dividend$m2 == divisor$m2 ) ||
               !(m1m1 || m1m2 || m2m1 || m2m2) )  {
              results <- compute.ratio(name=r$name,texlabel=r$texlabel, dividend=dividend,divisor=divisor,res=results)
            }
          }
        }
      }
    }
    
    # ratios from expressions
    cat("Computing ratio f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD\n")
    for( i in 1:length(mass_comb$sc[,1]) ) {
      fDs_indices <- which(results$val.tsboot$name == "f_Ds" & 
                          results$val.tsboot$m11 == mass_comb$sc[i,]$m1 & 
                          results$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
      mDs_indices <- which(results$val.tsboot$name == "m_Ds" & 
                          results$val.tsboot$m11 == mass_comb$sc[i,]$m1 & 
                          results$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
      fD_indices <- which(results$val.tsboot$name == "f_D" & 
                          results$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
      mD_indices <- which(results$val.tsboot$name == "m_D" & 
                          results$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
                          
      dividend <- compute.expression(expr=expression(a*sqrt(b)),
                    envir=data.frame( a=results$val.tsboot$val[fDs_indices], b=results$val.tsboot$val[mDs_indices] ),
                    m1=mass_comb$sc[i,]$m1, m2=mass_comb$sc[i,]$m2 )
      divisor <- compute.expression(expr=expression(a*sqrt(b)),
                    envir=data.frame( a=results$val.tsboot$val[fD_indices], b=results$val.tsboot$val[mD_indices] ),
                    m1=0.0009, m2=mass_comb$sc[i,]$m2 )
                    
      results <- compute.ratio(name="f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD",texlabel="$f_{D_s}\\sqrt{m_{D_s}}/f_D\\sqrt{m_D}$",
                               dividend=dividend, divisor=divisor, res=results )
    }

    cat("Computing ratio f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd\n")
    for( i in 1:length(mass_comb$sc[,1]) ) {
      fDs_indices <- which(results$val.tsboot$name == "f_Ds" & 
                          results$val.tsboot$m11 == mass_comb$sc[i,]$m1 & 
                          results$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
      mDs_indices <- which(results$val.tsboot$name == "m_Ds" & 
                          results$val.tsboot$m11 == mass_comb$sc[i,]$m1 & 
                          results$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
      fD_indices <- which(results$val.tsboot$name == "f_D" & 
                          results$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
      mD_indices <- which(results$val.tsboot$name == "m_D" & 
                          results$val.tsboot$m12 == mass_comb$sc[i,]$m2 )
                          
      dividend <- compute.expression(expr=expression(a*b^2),
                    envir=data.frame( a=results$val.tsboot$val[fDs_indices], b=results$val.tsboot$val[mDs_indices] ),
                    m1=mass_comb$sc[i,]$m1, m2=mass_comb$sc[i,]$m2 )
      divisor <- compute.expression(expr=expression(a*b^2),
                    envir=data.frame( a=results$val.tsboot$val[fD_indices], b=results$val.tsboot$val[mD_indices] ),
                    m1=0.0009, m2=mass_comb$sc[i,]$m2 )
                    
      results <- compute.ratio(name="f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd",texlabel="$f_{D_s}m_{D_s}^2/f_D{m_D}^2$",
                               dividend=dividend, divisor=divisor, res=results )
    }
    
    
    save(results,file="results.Rdata")
  } else {
    load("results.Rdata")
  } # if(recompute)
  
  if(overview) {
    # overview plots of data
    require(tikzDevice)
    filebase <- "data_overview"
    texfile <- sprintf("%s.tex",filebase) 
    pdffile <- sprintf("%s.pdf",filebase)
    tikz(texfile, standAlone = TRUE, width=5, height=5)
    for(name in names(quants)) {
      if(debug) print(name)
      #ts.indices <- which( results$res.tsboot$name == name )
      indices <- which( results$val$name == name )
      plotwitherror(y=results$val$val[indices], dy=results$val$dval[indices], x=results$val$m12[indices],
                    main=results$val$texlabel[indices[1]], xlab="$\\mu$",ylab=results$val$texlabel[indices[1]])
    }
    for(name in names(ratios)) {
      if(debug) print(name)
      #ts.indices <- which( results$res.tsboot$name == name )
      indices <- which( results$val$name == name )
      plotwitherror(y=results$val$val[indices], dy=results$val$dval[indices], x=results$val$m12[indices],
                    main=results$val$texlabel[indices[1]], xlab="$\\mu$",ylab=results$val$texlabel[indices[1]])
    }
    
    if(debug) print("f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD")
    indices <- which( results$val$name == "f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD" )
    plotwitherror(y=results$val$val[indices], dy=results$val$dval[indices], x=results$val$m12[indices],
                  main=results$val$texlabel[indices[1]], xlab="$\\mu$",ylab=results$val$texlabel[indices[1]])

    if(debug) print("f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd")
    indices <- which( results$val$name == "f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd" )
    plotwitherror(y=results$val$val[indices], dy=results$val$dval[indices], x=results$val$m12[indices],
                  main=results$val$texlabel[indices[1]], xlab="$\\mu$",ylab=results$val$texlabel[indices[1]])
    
    dev.off()
    tools::texi2dvi(texfile,pdf=T)
    command <- sprintf("pdfcrop %s %s",pdffile,pdffile)
    system(command)
  } # if(overview)
  
  # mu_s from FLAG ratio
  mu_s <- data.frame( val=c(0.0009*27.46), dval=c(0.44*0.0009) )
  
  # m_K_ov_f_K
  name <- "m_K_ov_f_K"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.0215, y=3.2 )
  mu_s_from_f_K <- match_mu.1D(name=name,alldat=results,masses=strange_masses,
              pheno=pheno,mu=mu_s,lg.coords=lg.coords,xlab="$\\mu_s$", solve=T)
  
  mu_s <- rbind( mu_s, data.frame(val=mu_s_from_f_K$predmu$mu, dval=mu_s_from_f_K$predmu$dmu) )
  
  # mu_c from HPQCD ratio
  mu_c <- data.frame( val=mu_s$val*11.85, 
                      dval= c( 0.0009*sqrt( (11.85*0.44)^2 + (27.46*0.16)^2 ), 
              sqrt( (mu_s$dval[2]*11.85)^2 + (mu_s$val[2]*0.16)^2 ) ) )

  cols <- c('black','red','forestgreen')
  syms <- c(1,15:(15+length(mu_c$val)-1))
  
  legend.mu_s <- list( labels=c("Data","$\\mu_s$ from FLAG ratio","$\\mu_s$ from $m_K/f_K$","$\\mu_s$ from $m_K/m_\\pi$"),
                       pch=c(syms,18), col=c(cols,'blue') )
  legend.mu_c <- list( labels=c("Data","$\\mu_c$ from FLAG$\\cdot$HPQCD ratios",
                                "$\\mu_c=$ HPQCD ratio $\\cdot\\mu_s$ from $m_K/f_K$","$\\mu_c$ from $m_D/m_\\pi$"),
                        pch=c(syms,18), col=c(cols,'blue') )
  
  name <- "m_K_ov_m_pi"
  lg.coords <- data.frame( x=0.021, y=3.77 )
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  mu_s_from_m_K <- match_mu.1D(name=name,alldat=results,masses=strange_masses,
              pheno=pheno,mu=mu_s,xlab="$a\\mu_s$", 
              lg=legend.mu_s,lg.coords=lg.coords,solve=T,ylim=c(3.56,3.77),
              xlim=c(0.021,0.027))
              
  mu_s <- rbind( mu_s, data.frame(val=mu_s_from_m_K$predmu$mu, dval=mu_s_from_m_K$predmu$dmu) )
              
  name <- "m_D_ov_m_pi"
  lg.coords <- data.frame( x=0.27, y=15.2 )
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  mu_c_from_m_D <- match_mu.1D(name=name,alldat=results,masses=charm_masses,
              lg=legend.mu_c,lg.coords=lg.coords,
              pheno=pheno,mu=mu_c, xlab="$a\\mu_c$", solve=T,
              xlim=c(0.27,0.36),ylim=c(12.8,15.2))
  
  mu_c <- rbind( mu_c, data.frame(val=mu_c_from_m_D$predmu$mu, dval=mu_c_from_m_D$predmu$dmu) )

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
  pred[[name]] <- match_mu.1D(name=name,alldat=results,masses=charm_masses,
              pheno=pheno,mu=mu_c,
              lg=legend.mu_c,lg.coords=lg.coords,xlab="$a\\mu_c$", solve=F,
              xlim=c(0.28,0.36))
              
  name <- "m_Ds_ov_f_Ds"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.0217, y=8.2 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=results,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c, lg=legend.mu_sc, lg.coords=lg.coords,ylim=c(7.2,8.2))              
  
  name <- "m_Ds_ov_m_K"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.0209, y=4.31 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=results,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c, lg=legend.mu_sc, lg.coords=lg.coords,ylim=c(3.7,4.31),xlim=c(0.021,0.027))  
  
  name <- "m_Ds_ov_m_D"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.0209, y=1.07 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=results,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c, lg=legend.mu_sc, lg.coords=lg.coords,
                     ylim=c(1.045,1.07),xlim=c(0.021,0.027))

  name <- "f_K_ov_f_pi"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.023, y=1.22 )
  pred[[name]] <- match_mu.1D(name=name,alldat=results,masses=strange_masses,
              pheno=pheno,mu=mu_s,
              lg=legend.mu_s,lg.coords=lg.coords,xlab="$a\\mu_s$",solve=F,
              ylim=c(1.19,1.22),xlim=c(0.023,0.027))

  name <- "f_D_ov_f_pi"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.255, y=1.8 )
  pred[[name]] <- match_mu.1D(name=name,alldat=results,masses=charm_masses,
              pheno=pheno,mu=mu_c,
              lg=legend.mu_c,lg.coords=lg.coords,xlab="$a\\mu_c$",
              ylim=c(1.5,1.8))                         

  name <- "m_D_ov_m_K"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.021, y=3.5 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=results,mass_comb=mass_comb,pheno=pheno,
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
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=results,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c, lg=legend.mu_sc, lg.coords=lg.coords,
                     xlim=c(0.023,0.027),
                     ylim=c(13.8,15.4))

  name <- "f_Ds_ov_f_pi"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.023, y=2.04 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=results,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c,lg=legend.mu_sc, lg.coords=lg.coords,
                     xlim=c(0.023,0.027),ylim=c(1.94,2.04))

  name <- "f_Ds_ov_f_D"
  pheno <- cbind(phys_ratios[ phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
  lg.coords <- data.frame( x=0.023, y=1.3 )
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=results,mass_comb=mass_comb,pheno=pheno,
                     mu_s, mu_c,lg=legend.mu_sc,lg.coords=lg.coords,
                     xlim=c(0.023,0.027),
                     ylim=c(1.16,1.3))

  name <- "m_K"
  pred[[name]] <- match_mu.1D(name=name,alldat=results,masses=strange_masses,
              mu=mu_s)
                     
  name <- "f_K"
  pred[[name]] <- match_mu.1D(name=name,alldat=results,masses=strange_masses,
              mu=mu_s)
                     
  name <- "m_D"
  pred[[name]] <- match_mu.1D(name=name,alldat=results,masses=charm_masses,
              mu=mu_c)
                     
  name <- "f_D"
  pred[[name]] <- match_mu.1D(name=name,alldat=results,masses=charm_masses,
              mu=mu_c)
                                   
  name <- "m_Ds"
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=results,mass_comb=mass_comb,
                  mu_s=mu_s, mu_c=mu_c)
  
  name <- "f_Ds"
  pred[[name]] <- match_mu_s_mu_c.2D(name=name,alldat=results,mass_comb=mass_comb,
                  mu_s=mu_s, mu_c=mu_c)
                  
  print(pred)
  
  write.csv(pred,file="preds.csv")
  
  options(stringsAsFactors = TRUE)
}
