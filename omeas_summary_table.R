omeas_summary_table <- function(filename = "omeas.summary.RData", 
                                outfile = "mpcac_v_kappa.dat"){
  load(filename)

  dat <- NULL
  for( ename in names(resultsum) ){
    ens <- resultsum[[ename]]

    dat <- rbind(dat,
                 data.frame(name = ename,
                            beta = ens$params$beta,
                            csw  = ens$params$csw,
                            kappa = ens$params$kappa,
                            L = ens$params$L,
                            T = ens$params$T,
                            mu = ens$params$mul,
                            muh = ens$params$muh,
                            musigma = ens$params$musigma,
                            mudelta = ens$params$mudelta,
                            Ntot    = ens$params$N.online + ens$params$skip,
                            Ntherm  = ens$params$N.online,
                            skip    = ens$params$skip,

                                  mpcac = ens$obs$mpcac_mc["val",1],
                                 dmpcac = ens$obs$mpcac_mc["dval",1],
                              tau_mpcac = ens$obs$mpcac_mc["tauint",1],
                             dtau_mpcac = ens$obs$mpcac_mc["dtauint",1],
                            
                                    mpi = ens$obs$mpi["val",1],
                                   dmpi = ens$obs$mpi["dval",1],
                                tau_mpi = ens$obs$mpi["tauint",1],
                               dtau_mpi = ens$obs$mpi["dtauint",1],
                            
                                    fpi = ens$obs$fpi["val",1],
                                   dfpi = ens$obs$fpi["dval",1],
                                tau_fpi = ens$obs$fpi["tauint",1],
                               dtau_fpi = ens$obs$fpi["dtauint",1],
                            
                                      P = ens$obs$P["val",1],
                                     dP = ens$obs$P["dval",1],
                                  tau_P = ens$obs$P["tauint",1],
                                 dtau_P = ens$obs$P["dtauint",1],

                                 deltaH = ens$obs$dH["val",1],
                                ddeltaH = ens$obs$dH["dval",1],
                             tau_deltaH = ens$obs$dH["tauint",1],
                            dtau_deltaH = ens$obs$dH["dtauint",1],

                                 exp_dH = ens$obs$expdH["val",1],
                                dexp_dH = ens$obs$expdH["dval",1],
                             tau_exp_dH = ens$obs$expdH["tauint",1],
                            dtau_exp_dH = ens$obs$expdH["dtauint",1]
                            )
                 )
  }
  
  Ls <- unique(dat$L)
  mus <- unique(dat$mu)
  kappas <- unique(dat$kappa)

  require("RColorBrewer")
  allclrs <- brewer.pal(n = 8,
                        name = "Dark2")
  allpchs <- c(0:8,15:16)

  dat <- cbind(dat, 
               data.frame(colour = allclrs[ match( dat$mu, mus ) ],
                          kappacolour = allclrs[ match( dat$kappa, kappas ) ], 
                          pch    = allpchs[ match( dat$L, Ls ) ])
               )

  write.table(x = dat,
              file = outfile,
              quote = TRUE,
              row.names = FALSE)
}

