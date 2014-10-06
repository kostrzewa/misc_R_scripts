pi0 <- function(disc_path,conn_path) {
  files <- getorderedfilelist(path=conn_path,basename="outprcvn.")
  cmicor <- readcmidatafiles(files)
  pion0.con <<- extract.obs(cmicor,  vec.obs=c(5), sign=c(-1))

  files <- getorderedfilelist(path=disc_path,basename="disc.0.1373.0.006.k0v4.", last.digits=4)
  print("Reading loops files")
  vdata <- readcmiloopfiles(files)
  print("Extracting loops")
  pi0loops <- extract.loop(vdata, obs=9)

  print("Computedisc")
  pi0disc.ll <- computeDisc(pi0loops, smeared=FALSE, real=TRUE, subtract.vev=TRUE)
  pi0disc.lf <- computeDisc(cf = pi0loops, cf2 = pi0loops, smeared=FALSE, smeared2=TRUE, real=TRUE, subtract.vev=TRUE)
  pi0disc.fl <- computeDisc(cf = pi0loops, cf2 = pi0loops, smeared=TRUE, smeared2=FALSE, real=TRUE, subtract.vev=TRUE)
  pi0disc.ff <- computeDisc(pi0loops, smeared=TRUE, real=TRUE, subtract.vev=TRUE)
  print("Concatenate")
  pion0.disc <<- c(pi0disc.ll, pi0disc.lf, pi0disc.fl, pi0disc.ff)
  print("Write disc to file")
  write.csv(pion0.disc$cf,file="pion0.disc.csv")
  pion0.full <<- add.cf(pion0.con, pion0.disc, a=0.5, b=2.)
  print("Write full to file")
  write.csv(pion0.full$cf,file="pion0.full.csv")
  print("Write conn to file")
  write.csv(pion0.con$cf,file="pion0.conn.csv")

  pion0.con <<- bootstrap.cf(cf=pion0.con,boot.R=400,boot.l=20,seed=21313)

  print("Connected matrixfit")
  pion0.con.matrixfit <<- matrixfit(cf=pion0.con,t1=11,t2=23,useCov=FALSE)
  summary(pion0.con.matrixfit)
  
  pion0.full <<- bootstrap.cf(cf=pion0.full,boot.R=400,boot.l=20,seed=12398)

  print("Full matrixfit")
  pion0.full.matrixfit <<- matrixfit(cf=pion0.full,t1=6,t2=20,useCov=FALSE)
  summary(pion0.full.matrixfit)

  pion0.full.effMass <<- bootstrap.effectivemass(cf=pion0.full,boot.R=400,boot.l=20,type='solve')
  pion0.full.effMass.fit <<- fit.effectivemass(cf=pion0.full.effMass,t1=5,t2=20,replace.na=T)
  summary(pion0.full.effMass.fit)
}



