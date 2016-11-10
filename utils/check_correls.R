check_correls <- function(ens,dirglob="ll"){
  for(E in ens){
    cat(sprintf("Ensemble %s\n",E))
    dirs <- Sys.glob(sprintf("%s/%s*",E,dirglob))
    for(dir in dirs){
      print(dir)
      fname <- sprintf("%s-%s",E,strsplit(dir,split='/')[[1]][2])
      files <- getorderedfilelist(path=dir,basename="outpr")
      cmicor <- readcmidatafiles(files=files,obs=c(1:20))
      cmicf <- extract.obs(cmicor=cmicor,vec.obs=c(1:20),symmetrise=FALSE)
      uwcf <- uwerr.cf(cf=cmicf)
      tikzfiles <- tikz.init(basename=fname,width=20,height=5)
      plotwitherror(y=uwcf$tauint,dy=uwcf$dtauint,x=c(1:length(uwcf$tauint)),
                    ylab="$\\tau_{\\mathrm{int}}(C_i(t))$")
      tikz.finalize(tikzfiles)
    }
  }
}
