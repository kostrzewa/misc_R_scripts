library('propagate')

source("/home/bartek/code/R/misc_R_scripts/fit_extrapolate_solve.R")
source("~/code/R/misc_R_scripts/hadron_obs.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/utils.R")
source("/home/bartek/code/R/misc_R_scripts/ratios_and_interpolations/definitions.R")

# this is a template for a driver script for the analysis involving interpolations
# in the strange and charm quark masses and subsequent extra-/interpolations
# for any other quantities

# in practice one would modify the sections between "EDIT FROM HERE" to "TO HERE" markers
# in order to customize the analysis to a different ensemble or different observables

# the individual steps of the analysis are documented further below in the analysis of
# m_K_ov_f_K

ratios_and_iterpolations_conn_meson <- function(analyses="all",debug=F,recompute=T,loadraw=T,overview=T) {
  # certain functionality relies on stuff being strings
  options(stringsAsFactors = FALSE)


  ### EDIT FROM HERE
  # masses to be used in this analysis
  analysis_name <- "iwa_b2.1-L48T96-k0.13729-mul0.0009"
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
  datanames <- list(ll_c=sprintf("ll_c.m%g.m%g.matrixfit",mass_comb$ll$m1,mass_comb$ll$m2),
                #ll_na=sprintf("lln_u_%g-d_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                #ll_nb=sprintf("lln_d_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ls_c=sprintf("ls_c.m%g.m%g.matrixfit",mass_comb$ls$m1,mass_comb$ls$m2),
                lc_c=sprintf("lc_c.m%g.m%g.matrixfit",mass_comb$lc$m1,mass_comb$lc$m2),
                sc_c=sprintf("sc_c.m%g.m%g.matrixfit",mass_comb$sc$m1,mass_comb$sc$m2) )

  ### TO HERE
  
  # single quantities
  quants <- define.meson.quants(datanames=datanames,light_masses=light_masses,strange_masses=strange_masses,charm_masses=charm_masses)

  # the different ratios that we would like to compute
  ratios <- define.meson.ratios(quants=quants)

  # read all the data files
  if(loadraw || recompute) {
    for( n in datanames ) {
      for( file in n ) {
        file <- sprintf("%s.Rdata",file)
        if(debug) cat("Loading" ,file,"\n")
        load(file)
      }
    }
  }

  # now we compute expectation values and errors from the bootstrap samples for the "quants" and "ratios"
  if(recompute) {
    if(debug) {
      cat("ratios_and_interpolations_conn_meson: computing hadron_obs!\n")
    }
    # all the data has been loaded in this environment, so we pass this along to the function
    hadron_obs <- compute.hadron_obs(envir=environment(),quants=quants,ratios=ratios,m.sea=m.sea,debug=debug)
    save(hadron_obs,file="hadron_obs.Rdata")
  } else {
    load("hadron_obs.Rdata")
  } # if(recompute)

  if(overview) {
    # overview plots of data
    require(tikzDevice)
    filebase <- "data_overview"
    temp <- sprintf("%s.%s",filebase,c("tex","pdf","aux","log"))
    tikzfiles <- list(tex=temp[1],pdf=temp[2],aux=temp[3],log=temp[4])
    rm(temp)
    tikz(tikzfiles$tex, standAlone = TRUE, width=5, height=5)

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

# the expressions are not quite ready yet,,,
#    if(debug) print("f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD")
#    indices <- which( hadron_obs$val$name == "f_Ds_sqrt_m_Ds_ov_f_D_sqrt_mD" )
#    plotwitherror(y=hadron_obs$val$val[indices], dy=hadron_obs$val$dval[indices], x=hadron_obs$val$m12[indices],
#                  main=hadron_obs$val$texlabel[indices[1]], xlab="$\\mu$",ylab=hadron_obs$val$texlabel[indices[1]])

#    if(debug) print("f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd")
#    indices <- which( hadron_obs$val$name == "f_Ds_m_Ds_sqrd_ov_f_D_mD_sqrd" )
#    plotwitherror(y=hadron_obs$val$val[indices], dy=hadron_obs$val$dval[indices], x=hadron_obs$val$m12[indices],
#                  main=hadron_obs$val$texlabel[indices[1]], xlab="$\\mu$",ylab=hadron_obs$val$texlabel[indices[1]])
    
    dev.off()
    tools::texi2dvi(tikzfiles$tex,pdf=T)                                                                                                                                                                           # use pdfcrop tool to remove plot borders
    command <- sprintf("pdfcrop %s %s",tikzfiles$pdf,tikzfiles$pdf)
    system(command)
    # remove temporary files 
    command <- sprintf("rm %s %s %s", tikzfiles$tex, tikzfiles$log, tikzfiles$aux)
    system(command)
  } # if(overview)

  ### EDIT FROM HERE
  # mu_s from FLAG ratio
  mu_s <- data.frame( val=c(0.0009*27.46), dval=c(0.44*0.0009), name="FLAG")
  
  #extrapolations <- data.frame(name=c(),val=c(),dval=c(),x=c(),dx=c())
  extrapolations <- list()
  
  #print(phys_ratios)  
  #readline("press key")

  # description of the most general type of analysis done here for a given quantity
  # we begin by enclosing the analysis in wavy brackets to provide a local namespace
  # and we define the name of the observable: m_K_ov_f_K in this case
  
  # m_K_ov_f_K
  name <- "m_K_ov_f_K"
  {
    # from the "phys_ratios" table, we load the physical values of the ratio corresponding to this "name"
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)

    # the bootstrap samples, expectation values and errors related to 'name' are stored in hadron_obs (see documentation
    # for  description of the data encapsulation format)
    # here we extract from the collection of observables the 
    
    # asemble the data into the correct format
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    # the observable depends on a number of parameters (masses), pred.idx is an index which determines which
    # mass value the observable is made to depend on, in this case m.val=2 specifies dependence on the second
    # valence mass, which is the strange mass for this observable since hadron_obs stores the mass values
    # in order of increasing size (this observable depends on the light mass and the strange mass)
    
    pred.idx <- list(m.val=c(2),m.sea=vector())

    # the extrapolation functions from the FES toolkit require the data in a data frame which is extracted
    # via extract.for.fes_fit
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)

    # do the linear interpolation of this data (in this case 1D, but in principle of arbitrary
    # dimensionality
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)

    # we are going to extrapolate now, so we choose some parameter values that we want to extrapolate towards
    # in this case the strange masses defined above
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval)
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    # add the result of this extrapolation to the list of extrapolations
    extrapolations[[length(extrapolations)+1]] <- cbind( 
                                data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                           dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                           pred,plot.x.idx='x1',plot.dx.idx='dx1')

    # we can also find the parameter values of the observable corresponding to the physical
    # value saved in pheno above, this is done here
    fes.solve <- fes_solve(fesfit=fes.fit,unknown='x1',known=c(),
                          y=phys_ratios[phys_ratios$name == name,]$val,
                          dy=phys_ratios[phys_ratios$name == name,]$dval )

    
    # and we store the result (which is a bootstrap sampling) in a data frame
    solution <- data.frame(val=mean(fes.solve[,1]), dval=sd(fes.solve[,1]), name=name)
    print(solution)

    # just as we extracted data for the fit above, here we extract it for plotting
    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))

    # and plot th data points with the phenomenological value indicated by a green band, th extrapolations indicated by coloured points
    # and the original data as black points
    plot.hadron_obs(df=df,name=name,pheno=pheno,solutions=solution,extrapolations=extrapolations[[length(extrapolations)]],
                    debug=debug, #labelx="$a\\mu_s$",
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel)

    # finally we add the resulting interpolated strange mass to the vector of strange masses
    mu_s <- rbind( mu_s, solution )
  }
    
  # some definitions required by many of the plots below
  cols <- c('black','red','forestgreen')
  syms <- c(1,15:(15+length(mu_s$val)))
  
  legend.mu_s <- list( labels=c("Data",
                                "$\\mu_s$ from FLAG ratio",
                                "$\\mu_s$ from $M_K/f_K$",
                                "$\\mu_s$ from $M_K/M_\\pi$"),
                       pch=syms, col=c(cols,'blue') )
  legend.mu_c <- list( labels=c("Data","$\\mu_c$ from FLAG$\\cdot$HPQCD ratios",
                                "$\\mu_c=$ HPQCD ratio $\\cdot\\mu_s$ from $M_K/f_K$",
                                "$\\mu_c=$ HPQCD ratio $\\cdot\\mu_s$ from $M_K/M_\\pi$",
                                "$\\mu_c$ from $M_D/M_\\pi$"),
                        pch=c(syms,18), col=c(cols,'cyan','blue') )
  
  legend.mu_sc <- list( labels=c("Data", "$\\mu_s$ and $\\mu_c$ from FLAG/HPQCD ratios",
                                 "$\\mu_c=$ HPQCD ratio $\\cdot\\mu_s$ from $M_K/f_K$",
                                 "$\\mu_c=$ HPQCD ratio $\\cdot\\mu_s$ from $M_K/M_\\pi$",
                                 "$\\mu_c$ from $M_D/M_\\pi$, $\\mu_s$ from $M_K/M_\\pi$"), 
                        pch=c(syms,18), col=c(cols,'cyan','blue') )
  
  name <- "m_K_ov_m_pi"
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    
    # asemble the data into the correct format
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)

    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval)
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)
    extrapolations[[length(extrapolations)+1]] <- cbind( 
                                data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                           dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                           pred,plot.x.idx='x1',plot.dx.idx='dx1')
    
    fes.solve <- fes_solve(fesfit=fes.fit,unknown='x1',known=c(),
                          y=phys_ratios[phys_ratios$name == name,]$val,
                          dy=phys_ratios[phys_ratios$name == name,]$dval )
    solution <- data.frame(val=mean(fes.solve[,1]), dval=sd(fes.solve[,1]), name=name)

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],solutions=solution,
                    #labelx="$a\\mu_s$",
                    lg=legend.mu_s,debug=debug,
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, ylim=c(3.58,3.75))
                    
    mu_s <- rbind( mu_s, solution )
  }

  # mu_c from HPQCD ratio, this will now contains three values of mu_c to which we will add a fourth by solving
  mu_c <- data.frame( val=mu_s$val*11.85, 
                      dval= c( 0.0009*sqrt( (11.85*0.44)^2 + (27.46*0.16)^2 ), 
              sqrt( (mu_s$dval[2:3]*11.85)^2 + (mu_s$val[2:3]*0.16)^2 ) ),
                      name=sprintf("%s/%s",mu_s$name,"HPQCD"))
  
  name <- "m_D_ov_m_pi"
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val,dx1=mu_c$dval)
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    fes.solve <- fes_solve(fesfit=fes.fit,unknown='x1',known=c(),
                          y=phys_ratios[phys_ratios$name == name,]$val,
                          dy=phys_ratios[phys_ratios$name == name,]$dval )

    solution <- data.frame(val=mean(fes.solve[,1]), dval=sd(fes.solve[,1]), name=name)

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],solutions=solution,
                    #labelx="$a\\mu_c$",
                    lg=legend.mu_c,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel,ylim=c(13,15))                

    mu_c <- rbind( mu_c, solution )
  }
  
  pred <- list()

  # redo m_K_ov_f_K for the purpose of extrapolating to the pair of quark masses from m_K and m_D
  name <- "m_K_ov_f_K"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval)
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')
                                       
    fes.solve <- fes_solve(fesfit=fes.fit,unknown='x1',known=c(),
                           y=phys_ratios[phys_ratios$name == name,]$val,
                           dy=phys_ratios[phys_ratios$name == name,]$dval )

    solution <- data.frame(val=mean(fes.solve[,1]), dval=sd(fes.solve[,1]), name=name)
                                   
    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],solutions=solution,
                    lg=legend.mu_s, debug=debug,
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, ylim=c(3.05,3.3), xlim=c(0.022,0.027))
  }
  
  name <- "m_D_ov_f_D"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val,dx1=mu_c$dval)
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')
                                       
    fes.solve <- fes_solve(fesfit=fes.fit,unknown='x1',known=c(),
                           y=phys_ratios[phys_ratios$name == name,]$val,
                           dy=phys_ratios[phys_ratios$name == name,]$dval )

    solution <- data.frame(val=mean(fes.solve[,1]), dval=sd(fes.solve[,1]), name=name)                                       
                                       
    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_c, ylim=c(7.8,10.3), xlim=c(0.26,0.40))
  }

  # we could also add this value to the list of possible charm masses
  # mu_c <- rbind( mu_c, data.frame(val=mean(fes.solve[,1]), dval=sd(fes.solve[,1]), name=name) )

  name <- "m_Ds_ov_f_Ds"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(1,2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc, ylim=c(7.1,8.2))
  }
  
  name <- "m_Ds_ov_m_K"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc,ylim=c(3.7,4.35))
  }
  
  name <- "m_Ds_ov_m_D"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)
    
  #   print(fes.extrapolate)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc,ylim=c(1.04,1.07))
  }

  name <- "m_Ds_ov_m_pi"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc,ylim=c(13.7,15.55))
  }

  name <- "m_D_ov_m_K"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc,ylim=c(3.5,4.2))
  }
    
  name <- "f_K_ov_f_pi"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    # asemble the data into the correct format
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval)
    
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)
    
    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')
    
    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_s,ylim=c(1.18,1.22))
  }

  name <- "f_D_ov_f_pi"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val[1:4],dx1=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_c, ylim=c(1.45,1.85))
  }  

  name <- "f_Ds_ov_f_pi"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc,ylim=c(1.94,2.15))
  }
  
  name <- "f_D_ov_f_K"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc,ylim=c(1.5,1.8))
  }
  
  name <- "f_Ds_ov_f_K"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc,ylim=c(1.55,1.85))
  }

  name <- "f_Ds_ov_f_D"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc,ylim=c(1.13,1.3))
  }

  name <- "f_Ds_ov_f_K"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc,ylim=c(1.15,1.3))
  }  
  
  name <- "m_K"
  if( analyses=="all" || any(analyses==name) )
  {
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval)

    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_s)
  }
                     
  name <- "f_K"
  if( analyses=="all" || any(analyses==name) )
  {
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval)

    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_s)
  }
                     
  name <- "m_D"
  if( analyses=="all" || any(analyses==name) )
  {
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val[1:4],dx1=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_c)
  }
                     
  name <- "f_D"
  if( analyses=="all" || any(analyses==name) )
  {
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val[1:4],dx1=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_c)
  }

  name <- "m_Ds"
  if( analyses=="all" || any(analyses==name) )
  {
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(1,2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)
  }
  
  name <- "f_Ds"
  if( analyses=="all" || any(analyses==name) )
  {
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(1,2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)
  }
  ### TO HERE
  
  # rename objects so that they refer to the ensemble that they are based on
  # and save to disk
  for( objname in c("mu_s","mu_c","extrapolations") ){
    savename <- gsub("-","_",sprintf("%s.%s",objname,analysis_name))
    assign(savename,get(objname))
    save(list=savename,file=sprintf("%s.Rdata",savename))
  }

  options(stringsAsFactors = TRUE)
}
