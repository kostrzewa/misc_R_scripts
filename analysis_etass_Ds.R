source("~/code/R/misc_R_scripts/lat_phys_ratios/compute_ratio.R")
#source("~/code/R/misc_R_scripts/predictNLS.R")
library("propagate")

analysis_etass_Ds <- function(debug=F,pause=F,read_from_file=F,compute_ratio=T,analyses_to_be_done_input,skip=0) {
  strange_masses <- c(0.021,0.023,0.025,0.027,0.031)
  charm_masses <- c(0.18,0.25,0.27,0.29,0.31,0.36,0.37,0.38)
  light_masses <- c(0.003)#,0.006,0.01)
  skip_strange_masses <- NULL
  skip_charm_masses <- NULL
  #skip_strange_masses <- c(0.031)
  #skip_charm_masses <- c(0.36,0.37,0.38)
  mu_c_ov_mu_s <- list(val=11.85, err=0.16)
  a <- c(0.082,0.092)
  da <- c(0.002,0.002)
  hbarcovera <- 197.3/a
  hbarc <- 197.3
  eta_ss_mass <- 686/hbarcovera
  d_eta_ss_mass <- 6/hbarcovera
  d_scale_eta_ss_mass <- eta_ss_mass*sqrt( (6/686)^2 + (da/a)^2 )

  cat( "\nPhenomenological value of the eta_ss mass in lattice units with scale uncertainty: \n" )
  cat( "Scale:\t", a, "+-", da, "fm\n" )
  cat( "Mass^2:\t", eta_ss_mass^2, "+-", 2*eta_ss_mass*d_scale_eta_ss_mass, "\n\n") 

  # make all possible strange and charm pairs
  D_s_combinations <- expand.grid( m1=strange_masses, m2=charm_masses )
  
  # because of the way tha expand.grid works, the first four entries in 'm1' will 
  # correspond to the correct strange masses for the four mass-diagonal
  # combinations of strange with strange
  eta_ss_combinations <- expand.grid( m1=strange_masses, m2=strange_masses )

  # the same is true for the light masses
  light_combinations <- expand.grid( m1=light_masses, m2=light_masses )

  # create vectors of directory names
  #eta_ss_dirs <- sprintf("strange_%s_strange_%s",strange_masses,strange_masses)
  #D_s_dirs <- sprintf("strange_%s_charm_%s",D_s_combinations$m1, D_s_combinations$m2 )
  eta_ss_dirs <- sprintf("ss_charged_s_p_%s-s_p_%s",strange_masses,strange_masses)
  #D_s_dirs <- sprintf("sc_s_p_%s-c_p_%s",D_s_combinations$m1, D_s_combinations$m2 )
  D_s_dirs <- sprintf("cs_c_p_%s-s_p_%s",D_s_combinations$m2, D_s_combinations$m1 )
  
  light_dirs <- sprintf("p0_u_%s-d_%s",light_masses,light_masses)
  light_dirs2 <- sprintf("p0_d_%s-u_%s",light_masses,light_masses)

  analyses <- NULL

  if(debug) {
    print(D_s_dirs)
    print(eta_ss_dirs)
  }

  # fit parameters for the various matrix and effective mass fits
  analyses[[1]] <- list( dirs=eta_ss_dirs, name="eta_ss_charged", 
                    mass_diagonal=T, q_masses=eta_ss_combinations, 
                    t1=12, t2=23, t1_plot=11, t2_plot=24, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[2]] <- list( dirs=D_s_dirs, name="D_s", 
                      mass_diagonal=F, q_masses=D_s_combinations, 
                      t1=15, t2=23, t1_plot=14, t2_plot=24, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[3]] <- list( dirs=light_dirs, name="ll_charged",
                         mass_diagonal=T, q_masses=light_combinations,
                         t1=10, t2=23, t1_plot=9, t2_plot=24, basename="outprcv.", observable=c(1) , sign=c(1) )
 
  analyses[[4]] <- list( dirs=light_dirs, name="ll_g5_neutral",
                        mass_diagonal=T, q_masses=light_combinations,
                        t1=12, t2=23, t1_plot=11, t2_plot=24, basename="outprcvn.", observable=c(1), sign=c(1) ) 

  analyses[[5]] <- list( dirs=light_dirs2, name="ll_g5_neutral",
                        mass_diagonal=T, q_masses=light_combinations,
                        t1=12, t2=23, t1_plot=11, t2_plot=24, basename="outprcvn.", observable=c(1), sign=c(1) ) 
 
  analyses[[6]] <- list( dirs=light_dirs, name="ll_11_neutral",
                         mass_diagonal=T, q_masses=light_combinations,
                         t1=12, t2=23, t1_plot=11, t2_plot=24, basename="outprcvn.", observable=c(5) , sign=c(-1) ) 

  # we will collect the meson masses here
  analysis_results <- NULL

  analyses_to_be_done <- NULL
  if( missing(analyses_to_be_done_input) ) {
    cat("Analyses to be done missing, doing all!\n")
    analyses_to_be_done <- seq(1,length(analyses))
  } else {
    analyses_to_be_done <- analyses_to_be_done_input
  }

  if( read_from_file == FALSE ) {
    for( ctr_analyses in analyses_to_be_done ) {
      masses <- data.frame(NULL)

      for( ctr_dirs in seq(1,length( analyses[[ctr_analyses]]$dirs ) ) ) {
        if(!file.exists( analyses[[ctr_analyses]]$dirs[ctr_dirs] )) {
          cat("## Skipping", analyses[[ctr_analyses]]$dirs[ctr_dirs], "because it doesn't exist!\n")
          next
        }
        fitresult <- do_meson_analysis(directory=analyses[[ctr_analyses]]$dirs[ctr_dirs], name=analyses[[ctr_analyses]]$name,debug=debug,
                                        basename=analyses[[ctr_analyses]]$basename, t1=analyses[[ctr_analyses]]$t1, t2=analyses[[ctr_analyses]]$t2,
                                        t1_plot=analyses[[ctr_analyses]]$t1_plot,t2_plot=analyses[[ctr_analyses]]$t2_plot,
                                        observable=analyses[[ctr_analyses]]$observable , sign=analyses[[ctr_analyses]]$sign,
                                        skip=skip )

        # interactive pausing between analyses to make debugging easier
        if(pause) {
          cat( "Press any key for next fit! \n" )
          readline()
        }

        masses <- rbind( masses, data.frame( m1=analyses[[ctr_analyses]]$q_masses$m1[ctr_dirs], 
                                             m2=analyses[[ctr_analyses]]$q_masses$m2[ctr_dirs],
                                             mass=fitresult$opt.res$par[1], dmass=sd(fitresult$massfit.tsboot[,1])))
      } # for(ctr_dirs)

      # add a new element to "analysis_results"
      analysis_results[[ctr_analyses]] <- list( res=masses, name=analyses[[ctr_analyses]]$name )

      pdf( sprintf("plots/%s.pdf",analyses[[ctr_analyses]]$name) , title=analyses[[ctr_analyses]]$name )
      par(family="Palatino")
      
      # for the eta_ss_charged we need a little more space on the plot
      ylim <- NULL
      if(analyses[[ctr_analyses]]$name == "eta_ss_charged"){
        ylim <- c(0.075,max(masses$mass^2)+0.05*max(masses$mass^2))
      } 
      plotwitherror(x=masses$m1,
                    y=(masses$mass)^2,
                    dy=2*masses$mass*masses$dmass,
                    xlab=expression(amu[s]),
                    ylab="(aM)^2",main=analyses[[ctr_analyses]]$name,
                    ylim=ylim,
                    xaxp=c(min(strange_masses)-2*0.002,max(strange_masses)+2*0.002,length(strange_masses)+4) )
      
      # in the case of the eta_ss_charged we would like to add the phenomenological value to the plot
      if(analyses[[ctr_analyses]]$name == "eta_ss_charged"){
        abline(h=eta_ss_mass^2)
        ytop_dm <- eta_ss_mass^2 + 2*eta_ss_mass*d_eta_ss_mass
        ybottom_dm <- eta_ss_mass^2 - 2*eta_ss_mass*d_eta_ss_mass
        top_ytop_dm_da <- eta_ss_mass^2 + 2*eta_ss_mass*d_scale_eta_ss_mass
        top_ybottom_dm_da <- ytop_dm
        bottom_ytop_dm_da <- ybottom_dm
        bottom_ybottom_dm_da <- eta_ss_mass^2 - 2*eta_ss_mass*d_scale_eta_ss_mass 
        col_dm <- c(rgb(0.6,0.0,0.0,0.3),rgb(0.0,0.6,0.0,0.3))
        col_dm_da <- c(rgb(0.0,0.0,0.6,0.3),rgb(0.6,0.0,0.6,0.3))
        rect(xleft=min(strange_masses)-0.003,xright=max(strange_masses)+0.003,ybottom=ybottom_dm,ytop=ytop_dm,col=col_dm,border=NA)
        rect(xleft=min(strange_masses)-0.003,xright=max(strange_masses)+0.003,ybottom=top_ybottom_dm_da,ytop=top_ytop_dm_da,col=col_dm_da,border=NA)
        rect(xleft=min(strange_masses)-0.003,xright=max(strange_masses)+0.003,ybottom=bottom_ybottom_dm_da,ytop=bottom_ytop_dm_da,col=col_dm_da,border=NA)
        legendlabels <- c(sprintf("Error on phenomenological value of eta_ss mass @ a=%s fm",a[2]), 
                          sprintf("Additional error from scale uncertainty a=%s +- %s fm ",a[2],da[2]),
                          sprintf("Error on phenomenological value of eta_ss mass @ a=%s fm",a[1]),                                                    
                          sprintf("Additional error from scale uncertainty a=%s +- %s fm",a[1],da[1]) )
        legend(x=(min(strange_masses)-0.002),y=0.0985,legend=legendlabels,col=c(col_dm[2],col_dm_da[2],col_dm[1],col_dm_da[1]),pch=c(15,15),bty="n")
      }
      dev.off()
    } # for(ctr_analyses)
    save( analysis_results, file="analysis_results.Rdata", list="analysis_results" )
  } else { # if(read_from_file)
    load(file="analysis_results.Rdata" , verbose=F)
  } # if(read_from_file)


  # compute a first ratio M(eta_ss)/M(D_s) without taking into account the correlations
  # between the two measurements
  if( compute_ratio == TRUE ) {
    # this will hold all computed ratios in a list
    ratios <- NULL
    ctr_ratios <- 0

    for( s_mass in strange_masses ) {
      # select eta_ss mass measurements for the given strange mass
      dividend <- list( val=analysis_results[[1]]$res$mass[ analysis_results[[1]]$res$m1==s_mass ]^2 , 
                        dval=2*analysis_results[[1]]$res$mass[ analysis_results[[1]]$res$m1==s_mass ]*
                                analysis_results[[1]]$res$dmass[ analysis_results[[1]]$res$m1==s_mass ] )
      for( c_mass in charm_masses ) {
        # select D_s mass measurements for the given strange and charm masses
        divisor <- list( val=analysis_results[[2]]$res$mass[ ( analysis_results[[2]]$res$m1==s_mass & analysis_results[[2]]$res$m2==c_mass ) ]^2,
                         dval=2*analysis_results[[2]]$res$mass[ ( analysis_results[[2]]$res$m1==s_mass & analysis_results[[2]]$res$m2==c_mass ) ]*
                                analysis_results[[2]]$res$dmass[ ( analysis_results[[2]]$res$m1==s_mass & analysis_results[[2]]$res$m2==c_mass ) ] )
        if(debug) {
          #cat("Building ratio for mus_s =", s_mass,", mu_c =", c_mass, "\n")
          #print(dividend)
          #print(divisor)
        }
        
        # if the given ratio does not exist, we skip this pair
        if( length(divisor$val) == 0 | length(dividend$val) == 0 ) {
          next
        }
        
        ctr_ratios <- ctr_ratios+1
        ratio <- compute_ratio(dividend=dividend, divisor=divisor, name="(M(eta_ss)/M(D_s))^2",debug=F)
        ratios <- rbind( ratios, data.frame(val=ratio$val, dval=ratio$dval, s_mass=s_mass, c_mass=c_mass ) )
      }
    }
    
    # if there are any data points which we don't want to include in the fit
    fitdata <- ratios[ !( (ratios$s_mass %in% skip_strange_masses) | (ratios$c_mass %in% skip_charm_masses) ), ] 
    if(debug) {
      cat("\nRatios selected for mu_s/mu_c fit\n")
      print(fitdata)
    }

    fitexpr <- expression( a*s_mass/(c_mass^2+b)+c )
    print(fitexpr)
    # fit the squared mass ratio in 2D
    fitfun <- function(s_mass,c_mass,a,b,c) {
      #return( a*s_mass/(c_mass^2+b) + c  )
      eval(fitexpr,envir=list(s_mass=s_mass,c_mass=c_mass,a=a,b=b,c=c) )
    }
    
    chiSqr <- function(par,ratios) {
      return( sum( ((ratios$val - fitfun(ratios$s_mass, ratios$c_mass, par[1], par[2], par[3]))/ratios$dval)^2 ) )
    }
    opt.res <- optim(par=c(0.1,0.1,0.1), fn=chiSqr, ratios=ratios)
   
    if(debug) {
      # TODO print ChiSqr + d.o.f. info
      print(opt.res)
    }
    
    # fit with 3 parameters and a ChiPT/HQET form
    fit_r_v_sc <- nls( #as.formula(value ~ eval(fitexpr,envir=list(s_mass=s_mass,c_mass=c_mass,a=a,b=b,c=c)))  
                      val ~ a*s_mass/(c_mass^2+b)+c #fitfun(s_mass,c_mass,a,b,c) 
                      , start=list(a=opt.res$par[1],b=opt.res$par[2],c=opt.res$par[3]), 
                      data=fitdata, weights=1/(fitdata$dval^2), control=list(tol=1e-6,minFactor=1e-10,maxiter=20000) )
    
    # linear approximation fit
    # fit_r_v_sc <- lm( value ~ s_mass + c_mass , data=fitdata, weights=1/(fitdata$error^2) )


    if(debug) {
      cat("\nChiSqr and parameter values from optim:\n", 
              "ChiSqr: ", opt.res$value, "Pars:", opt.res$par, "\n")
      print(summary(fit_r_v_sc))
      if( class(fit_r_v_sc) == 'lm' ) {
        cat("\nConfidence Intervals\n")
        print(confint(fit_r_v_sc))
      }
    }

    fit_r_v_sc_coefs <- coef(fit_r_v_sc)
    fit_r_v_sc_dcoefs <- coef(summary(fit_r_v_sc))[, "Std. Error"]

    cat("Pheno (M_ss/M_sc)^2",0.348^2," +- ",2*0.348*0.004,"\n") 
 
    library(tikzDevice)
    texfile <- "plots/M_eta_ss_over_M_D_s.tex"
    tikz(texfile, standAlone = TRUE, width=7, height=5)
  
    #pdf("plots/M_eta_ss_over_M_D_s.pdf",title="(M_ss/M_sc)^2",width=9,height=7)
    par(family="Palatino")
    # plot only data that was used for the fit
    plotwitherror(x=fitdata$s_mass,y=fitdata$val,dy=ratios$dval,main="$(M_{ss}/M_{sc})^2$",ylab="$(M(\\eta_{ss})/M(D_s))^2$",xlab="$a\\mu_s$",
                  #xlim=c(min(strange_masses)-0.01,max(strange_masses)+0.008),
                  xlim=c(0.012,max(strange_masses)+0.008),
                  ylim=c(min(fitdata$val)-0.01,max(fitdata$val)+0.07),                   
                  xaxp=c(min(strange_masses)-5*0.002,max(strange_masses)+4*0.002,length(strange_masses)+9) )

    # if not all data was used for the fit, add remaining points in red
    newdata <- ratios[ (ratios$s_mass %in% skip_strange_masses) | (ratios$c_mass %in% skip_charm_masses), ]
    print(newdata) 
    if(length(newdata[,1])>0){
      cat("Adding new data\n")
      plotwitherror(x=newdata$s_mass,y=newdata$val,dy=newdata$dval,rep=T,col="red",pch=16)
    }
   
        
    #let's add some extra charm masses
    pred_charm_masses <- c( 0.16, charm_masses)

    line_colours <- rainbow(length(pred_charm_masses)+2)
    # indicate fit lines and predictions on the plot with confidence intervals
    for( i in 1:length(pred_charm_masses) ) {
      pred_points <- data.frame( s_mass=seq(min(strange_masses)-0.01, max(strange_masses)+0.01, length.out=40), 
                                 c_mass=rep(pred_charm_masses[i],40) )
      preds <- NULL
      if( class(fit_r_v_sc) == 'lm' ) {
        cat("Using predict to make predictions\n")
        preds <- predict(fit_r_v_sc, newdata = pred_points, interval='confidence')
      } else if ( class(fit_r_v_sc) == 'nls' ) {
        cat("running predictNLS\n")
        preds <- predictNLS(fit_r_v_sc, newdata = pred_points, interval='confidence', do.sim=F)
      }

      # for the linear model we can do confidence intervals
      prediction <- NULL
      conf_band <- NULL
      if( class(fit_r_v_sc) == 'lm' ) {
        prediction <- preds[,1]
        conf_band <- c(rev(preds[,3]),preds[,2])
      # otherwise we have to approximate them by using the errors in the parameters  
      } else if ( class(fit_r_v_sc) == 'nls' ) {
        prediction <- preds$summary[,1]
        conf_band <- c(rev(preds$summary[,5]),preds$summary[,6])
      }
      lines(x=pred_points$s_mass,y=prediction,col=line_colours[i])
      polycol <- col2rgb(line_colours[i])/255
      polygon(c(rev(pred_points$s_mass), pred_points$s_mass), 
             conf_band, col = rgb(t(polycol),alpha=0.2), border=NA)
    }

    # compute a line and confidence bands along the line mu_c = mu_c_ov_mu_s*mu_s
    pred_s_masses <- seq(min(strange_masses)-0.015, max(strange_masses)+0.01, length.out=40)
    pred_points <- data.frame( s_mass=pred_s_masses, c_mass=pred_s_masses*mu_c_ov_mu_s$val
                              ,ds_mass=rep(0,length(pred_s_masses)), dc_mass=pred_s_masses*mu_c_ov_mu_s$err )
  
    if( class(fit_r_v_sc) == 'lm' ) {
      preds <- predict(fit_r_v_sc, newdata = pred_points, interval = 'confidence')
    } else if ( class(fit_r_v_sc) == 'nls' ) {
      preds <- predictNLS(fit_r_v_sc, newdata = pred_points, interval = 'confidence', do.sim=F)
    }
    
    prediction <- NULL
    conf_band <- NULL
    if( class(fit_r_v_sc) == 'lm' ){
      prediction <- preds[,1]
      conf_band <- c(rev(preds[,3]),preds[,2])
    } else if ( class(fit_r_v_sc) == 'nls' ) {
      prediction <- preds$summary[,1]
      conf_band <- c(rev(preds$summary[,6]),preds$summary[,5]) 
    }
    
    polycol <- col2rgb(line_colours[length(pred_charm_masses)+1])/255
    lines(x=pred_points$s_mass, y=prediction,col=line_colours[length(pred_charm_masses)+1]) 
    polygon(c(rev(pred_points$s_mass), pred_points$s_mass), conf_band, col = rgb(t(polycol),alpha=0.2))

    legend.xpos = min(strange_masses)-0.01
    legend.ypos = max(ratios$val)+0.08
    legend.labels = c( paste( paste("$a\\mu_c =",pred_charm_masses), "$") , "$a\\mu_c = 11.85(16) * a\\mu_s$", "$(M_{ss}/M_{sc})^2 = 0.121(3)$" )

    legend(x=legend.xpos,y=legend.ypos,legend=legend.labels,col=line_colours,lty=rep(1,length(pred_charm_masses)+2),bty="n")

    abline(h=0.348^2,col=line_colours[length(pred_charm_masses)+2])
    rect(xleft=min(strange_masses)-0.015,xright=max(strange_masses)+0.01,ytop=0.348^2+2*0.348*0.004,ybottom=0.348^2-2*0.348*0.004,
         col=rgb(t(col2rgb(line_colours[length(pred_charm_masses)+2])/255),alpha=0.3),border=NA )

    dev.off()
    tools::texi2dvi(texfile,pdf=T)

  }
  
}

do_meson_analysis <- function(directory,name,t1,t2,t1_plot,t2_plot,debug=F,basename="outprcv.",observable=c(1),sign=c(+1),skip=0) {
  if(debug) {
    cat(sprintf("Doing %s meson analysis in %s, observable %d\n",name,directory,observable))
  }

  files <- getorderedfilelist(basename=basename,path=directory)
  # prepend directory name
  if(debug) {
    cat("Reading", basename, "correlators from",directory,"\n")
  }
  cmicor <- readcmidatafiles(files[skip+1:length(files)])
  
  # when the data sample is still small, need to adjust boot.l
  boot.l <- 3
  if( length(files) < 100 ) {
    boot.l <- 2
  }

  if(debug) {
    cat("Performing matrixfit\n")
  }
  meson.cor <- extract.obs(cmicor, vec.obs=observable,sign.vec=sign)
  meson.cor <- bootstrap.cf(meson.cor,boot.R=400,boot.l=boot.l,seed=232415)
  meson.cor.matrixfit <- matrixfit(meson.cor, t1=t1, t2=t2, symmetrise=T, boot.R=400, boot.l=boot.l, useCov=F)
  summary(meson.cor.matrixfit)
  

  if(debug) {
    cat("Performing bootstrapped effective mass fit\n")
  }
  meson.cor.effectivemass <- bootstrap.effectivemass(meson.cor, boot.R=400, boot.l=boot.l, type="acosh")
  meson.cor.effectivemass <- fit.effectivemass(meson.cor.effectivemass, t1=t1, t2=t2, useCov=F)
  summary(meson.cor.effectivemass)

  ywidth <- 5*max( meson.cor.effectivemass$deffMass[t1:t2] )
  ymid <- mean( meson.cor.effectivemass$effMass[t1:t2] )

  plotfilename <- sprintf("plots/%s_%s.pdf",name,directory)
  pdf(plotfilename,onefile=T,title=paste(name, directory, sep=" "))
  plot(meson.cor.effectivemass, xlab=c("t/a"), ylab="aM", main=paste(name,directory,sep=" "), 
       xlim=c(t1_plot,t2_plot), ylim=c(ymid-ywidth,ymid+ywidth) )
  plot(meson.cor.matrixfit, xlab="t/a", ylab="C(t)",main=paste(name, directory,sep=" "))     
  dev.off();
  
  cat("\n")

  return(invisible( meson.cor.effectivemass ))
}

