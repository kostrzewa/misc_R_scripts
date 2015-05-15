library('propagate')

source("~/code/R/misc_R_scripts/fit_extrapolate_solve.R")
source("~/code/R/misc_R_scripts/hadron_obs.R")
source("~/code/R/misc_R_scripts/ratios_and_interpolations/utils.R")
source("~/code/R/misc_R_scripts/ratios_and_interpolations/definitions.R")

# this is a template for a driver script for the analysis involving interpolations
# in the strange and charm quark masses and subsequent extra-/interpolations
# for any other quantities

# in practice one would modify 

ratios_and_iterpolations_conn_meson <- function(analyses="all",debug=FALSE,recompute=TRUE,loadraw=TRUE,overview=TRUE,
												mc=TRUE, fitrange.systematic=TRUE, propagate.systematic=TRUE) {
  # certain functionality relies on stuff being strings
  options(stringsAsFactors = FALSE)

  ### EDIT FROM HERE
  # masses to be used in this analysis
  analysis_name <- "iwa_b2.1-L24T48-csw1.57551-k0.1373-mul0.003"
  strange_masses <- c(0.0224,0.0231,0.0238,0.0245,0.0252,0.0259)
  charm_masses <- c(0.2586,0.2704,0.2822,0.294,0.3058,0.3176)
  light_masses <- c(0.003)
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
  matrixfit.datanames <- list(ll_c=sprintf("ll_c.m%g.m%g.matrixfit",mass_comb$ll$m1,mass_comb$ll$m2),
                ll_nud=sprintf("ll_n_ud.m%g.m%g.matrixfit",mass_comb$ll$m1,mass_comb$ll$m2),
                ll_ndu=sprintf("ll_n_du.m%g.m%g.matrixfit",mass_comb$ll$m1,mass_comb$ll$m2),
                ls_c=sprintf("ls_c.m%g.m%g.matrixfit",mass_comb$ls$m1,mass_comb$ls$m2),
                lc_c=sprintf("lc_c.m%g.m%g.matrixfit",mass_comb$lc$m1,mass_comb$lc$m2),
                sc_c=sprintf("sc_c.m%g.m%g.matrixfit",mass_comb$sc$m1,mass_comb$sc$m2) )
  
  fitrange.datanames <- list(ll_c=sprintf("llc_u_%g_u_%g.fitrange",mass_comb$ll$m1,mass_comb$ll$m2),
                ll_nud=sprintf("lln_u_%g_d_%g.fitrange",mass_comb$ll$m1,mass_comb$ll$m2),
                ll_ndu=sprintf("lln_d_%g_u_%g.fitrange",mass_comb$ll$m1,mass_comb$ll$m2),
                ls_c=sprintf("ls_u_%g_sp_%g.fitrange",mass_comb$ls$m1,mass_comb$ls$m2),
                lc_c=sprintf("lc_u_%g_cp_%g.fitrange",mass_comb$lc$m1,mass_comb$lc$m2),
                sc_c=sprintf("sc_sp_%g_cp_%g.fitrange",mass_comb$sc$m1,mass_comb$sc$m2) )

  ### TO HERE
  
  # single quantities
  quants <- define.meson.quants(datanames=matrixfit.datanames,light_masses=light_masses,strange_masses=strange_masses,charm_masses=charm_masses)
  fitrange.quants <- define.meson.quants(datanames=fitrange.datanames,light_masses=light_masses,strange_masses=strange_masses,charm_masses=charm_masses)

  # the different ratios that we would like to compute
  ratios <- define.meson.ratios(quants=quants)
  fitrange.ratios <- define.meson.ratios(quants=fitrange.quants)

  # read all the data files
  if(loadraw || recompute) {
    for( file in unlist(matrixfit.datanames) )  {
      file <- sprintf("%s.Rdata",file)
      if(debug) cat("Loading" ,file,"\n")
      load(file)
    }
  }
  
  fitrange.systematic.errors <- NULL
  # compute the systematic error due to fitrange choice
  if(fitrange.systematic!=FALSE){
    if(fitrange.systematic=="compute"){
      fitrange.systematic.errors <- compute.fitrange_systematic(quants=fitrange.quants,ratios=fitrange.ratios,m.sea=m.sea,debug=debug,mc=mc)
      save(fitrange.systematic.errors,file="fitrange.systematic.errors.Rdata",compress=FALSE)
    } else {
      load("fitrange.systematic.errors.Rdata")
    }
  }

  # now we compute expectation values and errors from the bootstrap samples for the "quants" and "ratios"
  if(recompute) {
    if(debug) {
      cat("ratios_and_interpolations_conn_meson: computing hadron_obs!\n")
    }
    # all the data has been loaded in the current environment, so we pass this along to the function
    hadron_obs <- compute.hadron_obs(envir=environment(),quants=quants,ratios=ratios,m.sea=m.sea,debug=debug)
    if(fitrange.systematic!=FALSE){
      hadron_obs <- add.fitrange.systematic.hadron_obs(hadron_obs=hadron_obs,fitrange.serr=fitrange.systematic.errors)
      # free lots of memory
      rm(fitrange.systematic.errors)
    }
    save(hadron_obs,file="hadron_obs.Rdata",compress=FALSE)
    # free lots more memory
    for( dataname in unlist(matrixfit.datanames) ) {
      rm(dataname)
    }  
  } else {
    load("hadron_obs.Rdata")
    # corner case, stored hadron_obs without fitrange analysis, but fitrange estimate should be performed
    if( (fitrange.systematic != FALSE || propagate.systematic==TRUE) && is.null(hadron_obs[[1]]$fr)){
      hadron_obs <- add.fitrange.systematic.hadron_obs(hadron_obs=hadron_obs,fitrange.serr=fitrange.systematic.errors)
      # free lots of memory
      rm(fitrange.systematic.errors)
      save(hadron_obs,file="hadron_obs.Rdata",compress=FALSE)
    }
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
      obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
      obs.df <- extract.for.plot(hadron_obs=obs,x.name='m.val',x.idx=length(obs[[1]]$m.val))
      plotwitherror(y=obs.df$y, dy=obs.df$dy, mdy=obs.df$mdy, x=obs.df$x,
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
  
  # mu_s from iwa2.1-L48T96-csw1.57551-k0.13729-mul0.0009 analysis
  load("mu_s.iwa_b2.1_L48T96_k0.13729_mul0.0009.Rdata")
  mu_s <- mu_s.iwa_b2.1_L48T96_k0.13729_mul0.0009
  load("mu_s.sys.iwa_b2.1_L48T96_k0.13729_mul0.0009.Rdata")
  mu_s.sys <- mu_s.sys.iwa_b2.1_L48T96_k0.13729_mul0.0009

  load("mu_c.iwa_b2.1_L48T96_k0.13729_mul0.0009.Rdata")
  mu_c <- mu_c.iwa_b2.1_L48T96_k0.13729_mul0.0009
  load("mu_c.sys.iwa_b2.1_L48T96_k0.13729_mul0.0009.Rdata")
  mu_c.sys <- mu_c.sys.iwa_b2.1_L48T96_k0.13729_mul0.0009
  
  #extrapolations <- data.frame(name=c(),val=c(),dval=c(),x=c(),dx=c())
  extrapolations <- list()
  extrapolations.sys <- list()
  
  #print(phys_ratios)  
  #readline("press key")
  
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
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    
    # asemble the data into the correct format
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)

    fes.fit <- fes_fit_linear(mc=mc,dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval, x1.from=mu_s$name)
    fes.extrapolate <- fes_extrapolate(mc=mc,fesfit=fes.fit, pred=pred)
    extrapolations[[length(extrapolations)+1]] <- cbind( 
                                data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                           dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                           pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],
                    lg=legend.mu_s,debug=debug,
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel)#, ylim=c(3.58,3.75))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_s.sys$val, dx1=mu_s.sys$dval, x1.from=mu_s.sys$name)
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }
 
  name <- "m_D_ov_m_pi"
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val,dx1=mu_c$dval,x1.from=mu_c$name)
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],
                    #labelx="$a\\mu_c$",
                    lg=legend.mu_c,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel,ylim=c(13,15))                

    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_c.sys$val, dx1=mu_c.sys$dval, x1.from=mu_c.sys$name)
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }
  
  pred <- list()

  # redo m_K_ov_f_K for the purpose of extrapolating to the pair of quark masses from m_K and m_D
  name <- "m_K_ov_f_K"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval,x1.from=mu_s$name)
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')
                                       
    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],
                    lg=legend.mu_s, debug=debug,
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, ylim=c(3.05,3.3), xlim=c(0.022,0.027))
    
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_s.sys$val, dx1=mu_s.sys$dval)
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }
  
  name <- "m_D_ov_f_D"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val,dx1=mu_c$dval,x1.from=mu_c$name)
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')
                                       
    fes.solve <- fes_solve(fesfit=fes.fit,unknown='x1',known=c(),
                           y=phys_ratios[phys_ratios$name == name,]$val,
                           dy=phys_ratios[phys_ratios$name == name,]$dval )
    
    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_c)#, ylim=c(7.8,10.3), xlim=c(0.26,0.40))
    
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_c.sys$val, dx1=mu_c.sys$dval, x1.from=mu_c.sys$name)
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }

  name <- "m_Ds_ov_f_Ds"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(1,2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)#, ylim=c(7.1,8.2))


    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }
  
  name <- "m_Ds_ov_m_K"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind( 
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)#,ylim=c(3.7,4.35))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }
  
  name <- "m_Ds_ov_m_D"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)#,ylim=c(1.04,1.07))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }

  name <- "m_Ds_ov_m_pi"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(mc=mc, fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)#,ylim=c(13.7,15.55))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }

  name <- "m_D_ov_m_K"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)#,ylim=c(3.5,4.2))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }
    
  name <- "f_K_ov_f_pi"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    # asemble the data into the correct format
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(mc=mc,dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval,x1.from=mu_s$name)
    
    fes.extrapolate <- fes_extrapolate(mc=mc,fesfit=fes.fit, pred=pred)
    
    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')
    
    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_s)#,ylim=c(1.18,1.22))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_s.sys$val, dx1=mu_s.sys$dval, x1.from=mu_s.sys$name)
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }

  name <- "f_D_ov_f_pi"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(mc=mc,dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val[1:4],dx1=mu_c$dval[1:4],x1.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_c)#, ylim=c(1.45,1.85))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_c.sys$val[1:4], dx1=mu_c.sys$dval[1:4], x1.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }  

  name <- "f_Ds_ov_f_pi"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)#,ylim=c(1.94,2.15))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }
  
  name <- "f_D_ov_f_K"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(mc=mc,dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(mc=mc,fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)#,ylim=c(1.5,1.8))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }

  name <- "f_Ds_ov_f_D"
  if( analyses=="all" || any(analyses==name) )
  {
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],#solutions=solution,
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)#,ylim=c(1.13,1.3))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }

  name <- "f_Ds_ov_f_K"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    pheno <- cbind(phys_ratios[phys_ratios$name==name,],col=pheno.col,pch=pheno.pch)
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    
    pred.idx <- list(m.val=c(2,3),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(3))
    plot.hadron_obs(df=df,name=name,pheno=pheno,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)#,ylim=c(1.15,1.3))
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }  
  
  name <- "m_K"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F,mc=mc)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval,x1.from=mu_s$name)

    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred, mc=mc)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_s)
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_s.sys$val, dx1=mu_s.sys$dval, x1.from=mu_s.sys$name)
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }
                     
  name <- "f_K"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_s$val,dx1=mu_s$dval,x1.from=mu_s$name)

    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_s$",ylab=obs[[1]]$texlabel, lg=legend.mu_s)
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_s.sys$val, dx1=mu_s.sys$dval, x1.from=mu_s.sys$name)
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }
                     
  name <- "m_D"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val[1:4],dx1=mu_c$dval[1:4],x1.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_c)
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_c.sys$val[1:4], dx1=mu_c.sys$dval[1:4], x1.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }
                     
  name <- "f_D"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)
    pred.idx <- list(m.val=c(2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=mu_c$val[1:4],dx1=mu_c$dval[1:4],x1.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x1',plot.dx.idx='dx1')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_c)
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=mu_c.sys$val[1:4], dx1=mu_c.sys$dval[1:4], x1.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                  data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                  pred.sys, plot.x.idx='x1', plot.dx.idx='dx1')
    }
  }

  name <- "m_Ds"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(1,2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }
  
  name <- "f_Ds"
  if( analyses=="all" || any(analyses==name) )
  {
    cat(sprintf("Doing interpolations for %s\n",name))
    obs <- select.hadron_obs(hadron_obs,by='name',filter=name)

    pred.idx <- list(m.val=c(1,2),m.sea=vector())
    dat.fes <- extract.for.fes_fit(hadron_obs=obs,pred.idx=pred.idx)
    fes.fit <- fes_fit_linear(dat=dat.fes,debug=F)
    pred <- data.frame(x1=c(mu_s$val,mu_s$val[3]),x2=mu_c$val[1:4],
                       dx1=c(mu_s$dval,mu_s$dval[3]),dx2=mu_c$dval[1:4],
                       x1.from=c(mu_s$name,mu_s$name[3]),x2.from=mu_c$name[1:4])
    fes.extrapolate <- fes_extrapolate(fesfit=fes.fit, pred=pred)

    extrapolations[[length(extrapolations)+1]] <- cbind(
                            data.frame(name=name, texlabel=obs[[1]]$texlabel,  val=apply(X=fes.extrapolate$y,MARGIN=2,FUN=mean),
                                       dval=sqrt( apply(X=fes.extrapolate$y,MARGIN=2,FUN=sd)^2 + apply(X=fes.extrapolate$dy,MARGIN=2,FUN=mean)^2 ) ),
                                       pred,plot.x.idx='x2',plot.dx.idx='dx2')

    df <- extract.for.plot(hadron_obs=obs,x.name="m.val",x.idx=c(2))
    plot.hadron_obs(df=df,name=name,extrapolations=extrapolations[[length(extrapolations)]],
                    xlab="$a\\mu_c$",ylab=obs[[1]]$texlabel, lg=legend.mu_sc)
    if(propagate.systematic){
      cat(sprintf("Estimating systematic error for %s\n",name))
      dat.fes.sys <- extract.for.fes_fit(hadron_obs=obs, pred.idx=pred.idx, fr=TRUE)
      fes.fit.sys <- fes_fit_linear(dat=dat.fes.sys, debug=FALSE, mc=mc)
      pred.sys <- data.frame(x1=c(mu_s.sys$val,mu_s.sys$val[3]),x2=mu_c.sys$val[1:4],dx1=c(mu_s.sys$dval,mu_s.sys$dval[3]),dx2=mu_c.sys$dval[1:4],
                             x1.from=c(mu_s.sys$name,mu_s.sys$name[3]),x2.from=mu_c.sys$name[1:4])
      fes.extrapolate.sys <- fes_extrapolate(fesfit=fes.fit.sys, pred=pred.sys, mc=mc)
      save(fes.extrapolate.sys,file=sprintf("fes.extrapolate.sys.%s.Rdata",name),compress=FALSE)
      qt <- t(apply(X=fes.extrapolate.sys$y, MARGIN=2, FUN=quantile, probs=c(0,0.1573,0.5,0.8427,0.99), na.rm=TRUE))
      extrapolations.sys[[length(extrapolations.sys)+1]] <- cbind(
                                 data.frame(name=name, texlabel=obs[[1]]$texlabel,  qt0=qt[,1], qt1=qt[,2], median=qt[,3], qt2=qt[,4], qt3=qt[,5], mserr=abs(qt[,3]-qt[,2]), serr=abs(qt[,3]-qt[,4])),
                                 pred.sys, plot.x.idx='x2', plot.dx.idx='dx2')
    }
  }
  ### TO HERE
  
  # rename objects so that they refer to the ensemble that they are based on
  # and save to disk
  for( objname in c("extrapolations","extrapolations.sys") ){
    savename <- gsub("-","_",sprintf("%s.%s",objname,analysis_name))
    assign(savename,get(objname))
    cat("Saving",sprintf("%s.Rdata",savename),"\n")
    save(list=savename,file=sprintf("%s.Rdata",savename),compress=FALSE)
  }

  options(stringsAsFactors = TRUE)
}
