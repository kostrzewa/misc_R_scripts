# function for the automated analysis of connected meson two-point functions
# with the aim of extracting the (ground-state) masses and decay constants

# to adjust it for a given purpose the 
#   *_masses
#   mass_comb
#   dirs
#   analyses
# objects must be edited

source("~/code/R/misc_R_scripts/analysis_conn_meson_2pt/do_conn_meson_2pt_analysis.R")

analysis_conn_meson_2pt <- function(analyses_to_be_done_input,kappa,boot.R=400,boot.l=20,debug=F,pause=F,skip=0,seed=12345,useCov=F,read.cor=T,study.fitrange=F) {
  # masses to be used in this analysis
  strange_masses <- c(0.0224,0.0231,0.0238,0.0245,0.0252,0.0259)
  charm_masses <- c(0.2586,0.2704,0.2822)#,0.294,0.3058,0.3176)
  light_masses <- c(0.006)
  
  # combinations of these masses
  mass_comb <- list( ll=data.frame( m1=light_masses, m2=light_masses ),
                     ls=expand.grid( m1=light_masses, m2=strange_masses),
                     lc=expand.grid( m1=light_masses, m2=charm_masses),
                     sc=expand.grid( m1=strange_masses, m2=charm_masses) )
  
  correlators_dir <- "correlators/"

  # (sub-)directories where the correlators will be loaded from
  dirs <- list( ll_c=sprintf("llc_u_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ls_c=sprintf("ls_u_%g-sp_%g",mass_comb$ls$m1,mass_comb$ls$m2),
                lc_c=sprintf("lc_u_%g-cp_%g",mass_comb$lc$m1,mass_comb$lc$m2),
                sc_c=sprintf("sc_sp_%g-cp_%g",mass_comb$sc$m1,mass_comb$sc$m2),
                ll_n_ud=sprintf("lln_u_%g-d_%g",mass_comb$ll$m1,mass_comb$ll$m2),
                ll_n_du=sprintf("lln_d_%g-u_%g",mass_comb$ll$m1,mass_comb$ll$m2) )

  analyses <- NULL
  # the "analyses" list drives the automation process, the first four elements are self-explanatory
  #   t1,t2 are the fit range, t1_plot,t2_plot the x plotting limits
  #   basename is the correlator filename "prefix"
  #   observable is a numerical vector identifying which "gamma combinations" will be fitted (in the CMI format)
  #   sign is a numerical vector indicating whether the correlator is of "cosh" (+1) or "sinh" (-1) form
  analyses[[1]] <- list( dirs=dirs$ll_c, name="ll_c", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=12, t2=25, t1_plot=5, t2_plot=32, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[2]] <- list( dirs=dirs$ls_c, name="ls_c", mass_diagonal=F, q_masses=mass_comb$ls,
                         t1=10, t2=31, t1_plot=8, t2_plot=32, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[3]] <- list( dirs=dirs$lc_c, name="lc_c", mass_diagonal=F, q_masses=mass_comb$lc,
                         t1=10, t2=31, t1_plot=8, t2_plot=32, basename="outprcv.", observable=c(1), sign=c(1) )

  analyses[[4]] <- list( dirs=dirs$sc_c, name="sc_c", mass_diagonal=F, q_masses=mass_comb$sc,
                         t1=10, t2=31, t1_plot=8, t2_plot=32, basename="outprcv.", observables=c(1), sign=c(1) )
  
  analyses[[5]] <- list( dirs=dirs$ll_n_ud, name="ll_n_ud", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=10, t2=31, t1_plot=5, t2_plot=32, basename="outprcvn.", observable=c(5), sign=c(-1) )
  
  analyses[[6]] <- list( dirs=dirs$ll_n_du, name="ll_n_du", mass_diagonal=T, q_masses=mass_comb$ll,
                         t1=10, t2=31, t1_plot=5, t2_plot=32, basename="outprcvn.", observable=c(5), sign=c(-1) )

  analyses_to_be_done <- NULL
  if( missing(analyses_to_be_done_input) ) {
    cat("Analyses to be done missing, doing all!\n")
    analyses_to_be_done <- seq(1,length(analyses))
  } else {
    analyses_to_be_done <- analyses_to_be_done_input
  }

  analysis_results <- NULL
  for( ctr_analyses in analyses_to_be_done ) {
    for( ctr_dirs in seq(1,length( analyses[[ctr_analyses]]$dirs ) ) ) {
      if(!file.exists( analyses[[ctr_analyses]]$dirs[ctr_dirs] )) {
        cat("## Skipping", analyses[[ctr_analyses]]$dirs[ctr_dirs], "because it doesn't exist!\n")
        next
      }

      result <- do_conn_meson_2pt_analysis(directory=analyses[[ctr_analyses]]$dirs[ctr_dirs], name=analyses[[ctr_analyses]]$name,debug=debug, pause=pause,
                                  basename=analyses[[ctr_analyses]]$basename, t1=analyses[[ctr_analyses]]$t1, t2=analyses[[ctr_analyses]]$t2,
                                  t1_plot=analyses[[ctr_analyses]]$t1_plot,t2_plot=analyses[[ctr_analyses]]$t2_plot,
                                  observable=analyses[[ctr_analyses]]$observable , sign=analyses[[ctr_analyses]]$sign,
                                  skip=skip, kappa=kappa, q_masses=analyses[[ctr_analyses]]$q_masses[ctr_dirs,],
                                  boot.R=boot.R, boot.l=boot.l, seed=seed, useCov=useCov, read.cor=read.cor, study.fitrange=study.fitrange
                                 )
          
      analysis_results <- rbind(analysis_results, result)
    }
  }

  write.csv(analysis_results, file="meson_2pt.csv")
}
