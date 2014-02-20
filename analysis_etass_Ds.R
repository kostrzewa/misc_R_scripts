source("~/code/R/misc_R_scripts/lat_phys_ratios/compute_ratio.R")

analysis_etass_Ds <- function(debug=F,t1,t2,t1_plot,t2_plot,pause=F,read_from_file=F) {
  strange_masses <- c(0.021,0.023,0.025,0.027)
  charm_masses <- c(0.27,0.29,0.31)
  #strange_masses <- c(0.021, 0.023)
  #charm_masses <- c(0.27)
  a <- 0.095
  hbarcovera <- 197.3/a
  eta_ss_mass <- 686/hbarcovera
  d_eta_ss_mass <- 6/hbarcovera

  local_t1_plot <- 0
  local_t2_plot <- 0

  if(missing(t1_plot)){
    local_t1_plot <- t1-1
  } else {
    local_t1_plot <- t1_plot
  }

  if(missing(t2_plot)){
    local_t2_plot <- t2+1
  } else {
    local_t2_plot <- t2_plot
  }

  # make all possible strange and charm pairs
  D_s_combinations <- expand.grid( m1=strange_masses, m2=charm_masses )
  
  # because of the way tha expand.grid works, the first four entries in 'm1' will 
  # correspond to the correct strange masses for the four mass-diagonal
  # combinations of strange with strange
  eta_ss_combinations <- expand.grid( m1=strange_masses, m2=strange_masses )

  # create vectors of directory names
  eta_ss_dirs <- sprintf("strange_%s_strange_%s",strange_masses,strange_masses)
  D_s_dirs <- sprintf("strange_%s_charm_%s",D_s_combinations$m1, D_s_combinations$m2 )

  analyses <- NULL

  analyses[[1]] <- list( dirs=eta_ss_dirs, name="eta_ss_charged", 
                    mass_diagonal=T, q_masses=eta_ss_combinations, 
                    t1=12, t2=23, t1_plot=11, t2_plot=24, basename="outprcv." )
  
  analyses[[2]] <- list( dirs=D_s_dirs, name="D_s", 
                      mass_diagonal=F, q_masses=D_s_combinations, 
                      t1=15, t2=23, t1_plot=14, t2_plot=24, basename="outprcvn." ) 

# analyses[[2]] <- list( dirs=eta_ss_dirs, name="eta_ss_neutral", 
#                   mass_diagonal=T, q_masses=eta_ss_combinations, 
#                   t1=12, t2=23, t1_plot=11, t2_plot=24, basename="outprcvn." )

  # we will collect the meson masses here
  analysis_results <- NULL

  if( read_from_file == FALSE ) {
    for( ctr_analyses in seq(1,length(analyses)) ) {
      masses <- data.frame(NULL)
      for( ctr_dirs in seq(1,length( analyses[[ctr_analyses]]$dirs ) ) ) {
        effmassfit <- do_meson_analysis(directory=analyses[[ctr_analyses]]$dirs[ctr_dirs], name=analyses[[ctr_analyses]]$name,debug=debug,
                                        basename=analyses[[ctr_analyses]]$basename, t1=analyses[[ctr_analyses]]$t1, t2=analyses[[ctr_analyses]]$t2,
                                        t1_plot=analyses[[ctr_analyses]]$t1_plot,t2_plot=analyses[[ctr_analyses]]$t2_plot)

        # interactive pausing between analyses to make debugging easier
        if(pause) {
          cat( "Press any key for next fit! \n" )
          readline()
        }

        masses <- rbind( masses, data.frame( m1=analyses[[ctr_analyses]]$q_masses$m1[ctr_dirs], 
                                             m2=analyses[[ctr_analyses]]$q_masses$m2[ctr_dirs],
                                             mass=effmassfit$opt.res$par[1], dmass=sd(effmassfit$massfit.tsboot[,1]) ) )
      } # for(ctr_dirs)

      analysis_results[[ctr_analyses]] <- list( res=masses, name=analyses[[ctr_analyses]]$name )

      pdf( sprintf("plots/%s.pdf",analyses[[ctr_analyses]]$name) , title=analyses[[ctr_analyses]]$name )
      plotwitherror(x=masses$m1,y=(masses$mass)^2,dy=2*masses$mass*masses$dmass,xlab=expression(amu[s]),ylab="(aM)^2",main=analyses[[ctr_analyses]]$name )
      dev.off()
    } # for(ctr_analyses)
    save( analysis_results, file="analysis_results.Rdata", list="analysis_results" )
  } else { # if(read_from_file == TRUE)
    load(file="analysis_results.Rdata" , verbose=F)
  }

  # compute a first ratio M(eta_ss)/M(D_s) without taking into account the correlations
  # between the two measurements

  # this will hold all computed ratios in a list
  ratios <- NULL
  ctr_ratios <- 0

  for( s_mass in strange_masses ) {
    # select eta_ss mass measurements for the given strange mass
    dividend <- list( value=analysis_results[[1]]$res$mass[ analysis_results[[1]]$res$m1==s_mass ] , 
                      error=analysis_results[[1]]$res$dmass[ analysis_results[[1]]$res$m1==s_mass ] )
    for( c_mass in charm_masses ) {
      ctr_ratios <- ctr_ratios+1    
     
      # select D_s mass measurements for the given strange and charm masses
      divisor <- list( value=analysis_results[[2]]$res$mass[ ( analysis_results[[2]]$res$m1==s_mass & analysis_results[[2]]$res$m2==c_mass ) ],
                       error=analysis_results[[2]]$res$dmass[ ( analysis_results[[2]]$res$m1==s_mass & analysis_results[[2]]$res$m2==c_mass ) ] )
      
      ratios <- rbind( ratios,  data.frame( compute_ratio(dividend=dividend, divisor=divisor, name="M(eta_ss)/M(D_s)",debug=F), 
                                        s_mass=s_mass, c_mass=c_mass ) )
    }
  }

  pdf("plots/M_eta_ss_over_M_D_s.pdf",title=as.character(ratios$name[1]))
  plotwitherror(x=ratios$s_mass,y=ratios$value,dy=ratios$error,main=ratios$name[1],ylab="M(eta_ss)/M(D_s)",xlab=expression(amu[s]))
  abline(h=0.348)
  rect(xleft=min(strange_masses),xright=max(strange_masses),ytop=0.352,ybottom=0.344,col=rgb(0.6,0.0,0.0,0.3) )
  dev.off()

}

do_meson_analysis <- function(directory,name,t1,t2,t1_plot,t2_plot,debug=F,basename="outprcv.",observable=c(1)) {
  if(debug) {
    cat(sprintf("Doing %s meson analysis for %s\n",name,directory))
  }

  plotfilename <- sprintf("plots/%s_%s.pdf",name,directory)
  pdf(plotfilename,onefile=T,title=paste(name, directory, sep=" "))
  
  files <- getorderedfilelist(basename=basename,path=directory)
  # prepend directory name
  if(debug) {
    cat("Reading", basename, "correlators from",directory,"\n")
  }
  cmicor <- readcmidatafiles(files)

  if(debug) {
    cat("Performing matrixfit\n")
  }
  meson.cor <- extract.obs(cmicor, vec.obs=c(1))
  meson.cor.matrixfit <- matrixfit(meson.cor, t1=t1, t2=t2, symmetrise=T, boot.R=400, boot.l=1, useCov=F)
  plot(meson.cor.matrixfit, xlab="t/a", ylab="C(t)",main=paste(name, directory,sep=" "))

  meson.cor.effectivemass <- bootstrap.effectivemass(meson.cor, boot.R=400, boot.l=1, type="acosh")
  meson.cor.effectivemass <- fit.effectivemass(meson.cor.effectivemass, t1=t1, t2=t2, useCov=F)

  ywidth <- 5*max( meson.cor.effectivemass$deffMass[t1:t2] )
  ymid <- mean( meson.cor.effectivemass$effMass[t1:t2] )

  plot(meson.cor.effectivemass, xlab=c("t/a"), ylab="aM", main=paste(name,directory,sep=" "), 
       xlim=c(t1_plot,t2_plot), ylim=c(ymid-ywidth,ymid+ywidth) )
  dev.off();
  
  summary(meson.cor.matrixfit)
  summary(meson.cor.effectivemass)
  cat("\n")

  return(invisible(meson.cor.effectivemass))
}

