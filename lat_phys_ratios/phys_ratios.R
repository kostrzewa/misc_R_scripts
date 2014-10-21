source("~/code/R/misc_R_scripts/lat_phys_ratios/compute_ratio.R")

phys_ratios <- function(file) {
  dat <- read.table(file,header=T)
  combinations <- expand.grid( dividend=dat$name, divisor=dat$name )
  ratios <- NULL
  
  
  for( i in 1:length(combinations[,1]) ) {
    # skip combinations which are 1
    if( combinations[i,]$dividend == combinations[i,]$divisor ) {
      next()
    }
    
    name <- sprintf("%s_ov_%s", combinations[i,]$dividend, combinations[i,]$divisor )
    
    dividend.idx <- which( dat$name == combinations[i,]$dividend )
    divisor.idx <- which( dat$name == combinations[i,]$divisor )
    
    # extract data from dat
    dividend <- NULL
    divisor <- NULL
    
    # prefer isospin symmetric, then PDG, then FLAG
    if( !is.na( dat[dividend.idx,]$iso ) ) {
      dividend <- data.frame( val=dat[dividend.idx,]$iso, dval=dat[dividend.idx,]$diso, type="ISOSYM" )
    } else if( !is.na( dat[dividend.idx,]$pdg ) ) {
      dividend <- data.frame( val=dat[dividend.idx,]$pdg, dval=dat[dividend.idx,]$dpdg, type="PDG" )
    } else {
      dividend <- data.frame( val=dat[dividend.idx,]$flag, dval=dat[dividend.idx,]$dflag, type="FLAG" )
    }

    if( !is.na( dat[divisor.idx,]$iso ) ) {
      divisor <- data.frame( val=dat[divisor.idx,]$iso, dval=dat[divisor.idx,]$diso, type="ISOSYM" )
    } else if( !is.na( dat[divisor.idx,]$pdg ) ) {
      divisor <- data.frame( val=dat[divisor.idx,]$pdg, dval=dat[divisor.idx,]$dpdg, type="PDG" )
    } else {
      divisor <- data.frame( val=dat[divisor.idx,]$flag, dval=dat[divisor.idx,]$dflag, type="FLAG" )
    }
    
    if( divisor$type == dividend$type ) {
      type <- divisor$type
    } else {
      type <- sprintf("%s/%s",dividend$type,divisor$type)
    }
    
    ratios <- rbind(ratios,c(compute_ratio(name=name,dividend=dividend,divisor=divisor),type=type))
  }
  print(ratios)
  write.csv(ratios,file="phys_ratios.csv")
}
