pi0 <- function(pion0.disc_path, pion0.conn_path, pionc.conn_path, pionc.con.fitrange, pion0.con.fitrange, pion0.full.fitrange, 
                kappa, mul, skip.conf, skipfiles=0, boot.R=1800, boot.l=15, read.data=TRUE, seed=c(34872), useCov=FALSE) {
  pionc.con <- NULL
  pion0.con <- NULL
  pion0.disc <- NULL

  if(read.data) {
    files <- getorderedfilelist(path=pionc.conn_path,basename="outprcv.")
    if(!missing(skip.conf))
      files <- files[-skip.conf]
    cmicor <- readcmidatafiles(files[(skipfiles+1):length(files)])
    pionc.con <- extract.obs(cmicor,  vec.obs=c(1), sign=c(1))

    files <- getorderedfilelist(path=pion0.conn_path,basename="outprcvn.")
    if(!missing(skip.conf))
      files <- files[-skip.conf]
    cmicor <- readcmidatafiles(files[(skipfiles+1):length(files)])
    pion0.con <- extract.obs(cmicor,  vec.obs=c(5), sign=c(-1))

    files <- getorderedfilelist(path=pion0.disc_path,basename=sprintf("disc.%s.%s.k0v4.",kappa,mul), last.digits=4)
    if(!missing(skip.conf))
      files <- files[-skip.conf]
    print("Reading loops files")
    vdata <- readcmiloopfiles(files[(skipfiles+1):length(files)])
    print("Extracting loops")
    pi0loops <- extract.loop(vdata, obs=9)
    
    print("Computedisc")
    pi0disc.ll <- computeDisc(pi0loops, smeared=FALSE, real=TRUE, subtract.vev=TRUE)
    pi0disc.lf <- computeDisc(cf = pi0loops, cf2 = pi0loops, smeared=FALSE, smeared2=TRUE, real=TRUE, subtract.vev=TRUE)
    pi0disc.fl <- computeDisc(cf = pi0loops, cf2 = pi0loops, smeared=TRUE, smeared2=FALSE, real=TRUE, subtract.vev=TRUE)
    pi0disc.ff <- computeDisc(pi0loops, smeared=TRUE, real=TRUE, subtract.vev=TRUE)
    print("Concatenate")
    pion0.disc <- c(pi0disc.ll, pi0disc.lf, pi0disc.fl, pi0disc.ff)

    save(pionc.con,file="pionc.con.Rdata")
    save(pion0.con,file="pion0.con.Rdata")
    save(pion0.disc, file="pion0.disc.Rdata")
  } else {
    cat("Warning: reading data from .Rdata files, if you updated anything, set read.data=TRUE!\n")
    load("pionc.con.Rdata")
    load("pion0.con.Rdata")
    load("pion0.disc.Rdata")
  }

  pionc.con <- bootstrap.cf(cf=pionc.con,boot.R=boot.R,boot.l=boot.l,seed=seed)
  cat("Conn. PiC")
  pionc.con.matrixfit <<- matrixfit(cf=pionc.con,t1=pionc.con.fitrange[1],t2=pionc.con.fitrange[2],boot.R=boot.R,boot.l=boot.l,useCov=useCov)
  summary(pionc.con.matrixfit)
  save(pionc.con.matrixfit,file="pionc.con.matrixfit.Rdata")

  pionc.con.effMass <- bootstrap.effectivemass(cf=pionc.con,boot.R=boot.R,boot.l=boot.l,type='solve')
  pionc.con.effMass.fit <<- fit.effectivemass(cf=pionc.con.effMass,t1=pionc.con.fitrange[1],t2=pionc.con.fitrange[2],replace.na=T,useCov=useCov)
  summary(pionc.con.effMass.fit)
  save(pionc.con.effMass.fit,file="pionc.con.effMass.fit")

  pion0.full <- add.cf(pion0.con, pion0.disc, a=0.5, b=2.)
  pion0.con <- bootstrap.cf(cf=pion0.con,boot.R=boot.R,boot.l=boot.l,seed=seed)

  print("Conn. Pi0")
  pion0.con.matrixfit <<- matrixfit(cf=pion0.con,t1=pion0.con.fitrange[1],t2=pion0.con.fitrange[2],boot.R=boot.R,boot.l=boot.l,useCov=useCov)
  summary(pion0.con.matrixfit)
  save(pion0.con.matrixfit,file="pion0.com.matrixfit.Rdata")
  
  pion0.con.effMass <- bootstrap.effectivemass(cf=pion0.con,boot.R=boot.R,boot.l=boot.l,type='solve')
  pion0.con.effMass.fit <<- fit.effectivemass(cf=pion0.con.effMass,t1=pion0.con.fitrange[1],t2=pion0.con.fitrange[2],replace.na=T,useCov=useCov)
  summary(pion0.con.effMass.fit)
  save(pion0.con.effMass.fit,file="pion0.con.effMass.fit.Rdata")
  
  pion0.full <- bootstrap.cf(cf=pion0.full,boot.R=boot.R,boot.l=boot.l,seed=seed)

  print("Full Pi0")
  pion0.full.matrixfit <<- matrixfit(cf=pion0.full,t1=pion0.full.fitrange[1],t2=pion0.full.fitrange[2],boot.R=boot.R,boot.l=boot.l,useCov=useCov)
  summary(pion0.full.matrixfit)
  save(pion0.full.matrixfit,file="pion0.full.matrixfit.Rdata")

  pion0.full.effMass <- bootstrap.effectivemass(cf=pion0.full,boot.R=boot.R,boot.l=boot.l,type='solve')
  pion0.full.effMass.fit <<- fit.effectivemass(cf=pion0.full.effMass,t1=pion0.full.fitrange[1],t2=pion0.full.fitrange[2],replace.na=T,useCov=useCov)
  summary(pion0.full.effMass.fit)
  save(pion0.full.effMass.fit,file="pion0.effMass.fit.Rdata")

  print("Disc Pi0")
  pion0.disc <- mul.cf(cf=pion0.disc,a=4)
  pion0.disc <- bootstrap.cf(cf=pion0.disc,boot.R=boot.R,boot.l=boot.l,seed=seed)
  pion0.disc.matrixfit <<- matrixfit(cf=pion0.disc,t1=pion0.full.fitrange[1],t2=pion0.full.fitrange[2],boot.R=boot.R,boot.l=boot.l,useCov=useCov)
  summary(pion0.disc.matrixfit)
  save(pion0.disc.matrixfit,file="pion0.disc.matrixfit.Rdata")

  pion0.disc.effMass <- bootstrap.effectivemass(cf=pion0.disc,boot.R=boot.R,boot.l=boot.l,type='solve')
  pion0.disc.effMass.fit <<- fit.effectivemass(cf=pion0.disc.effMass,t1=pion0.full.fitrange[1],t2=pion0.full.fitrange[2],replace.na=T,useCov=useCov)
  summary(pion0.disc.effMass.fit)
  save(pion0.disc.effMass.fit,file="pion0.disc.effMass.fit.Rdata")

}
