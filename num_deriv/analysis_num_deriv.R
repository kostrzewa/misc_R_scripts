analysis_num_deriv <- function(dir,volume=2^4,int.steps=101,trajs=1,endian="little") {
  pat <- "*f_numerical.bin"
  #anpattern <- "*f_analytical.bin"
  # extract precisions from filenames 
  files <- list.files(dir,pattern=pat,full.names=FALSE)
  tokens <- strsplit(x=files,split="_")
  tokens <- array(data=unlist(tokens),dim=c(3,length(tokens)))
  precs <- tokens[1,]

  print(precs)

  df_a <- list()
  df_n <- list()

  for( prec in precs ) {
    print(prec)
    n_filename <- sprintf("%s_f_numerical.bin",prec)
    a_filename <- sprintf("%s_f_analytical.bin",prec)
    
    n_file <- file(n_filename,"rb")
    a_file <- file(a_filename,"rb")

    df_n[[prec]] <- t(array(data=readBin(n_file,double(),n=trajs*int.steps*volume*4*8,endian=endian),dim=c(volume*4*8,int.steps*trajs)))
    df_a[[prec]] <- t(array(data=readBin(a_file,double(),n=trajs*int.steps*volume*4*8,endian=endian),dim=c(volume*4*8,int.steps*trajs)))

    close(n_file)
    close(a_file)
  }
  
  pdf(file="diff_a.pdf",onefile=T,width=6,height=6)
  for( i in 1:22 ) {
    plot(df_a[[length(df_a)]][,1]-df_a[[i]][,1])
  }
  dev.off()

  pdf(file="diff_n.pdf",onefile=T,width=6,height=6)
  for( i in 1:22 ) {
    plot(df_n[[length(df_n)]][,1]-df_n[[i]][,1])
  }
  dev.off()
  
  pdf(file="diff_a_n.pdf",onefile=T,width=6,height=6)
  for( i in 1:22 ) {
    plot(df_a[[length(df_a)]][,1]-df_n[[i]][,1])
  }
  dev.off()
 
}
