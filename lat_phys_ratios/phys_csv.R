source("~/code/R/misc_R_scripts/lat_phys_ratios/compute_ratio.R")

phys_csv <- function(file,debug=FALSE,all=FALSE) {
  # certain functionality relies on stuff being strings
  options(stringsAsFactors = FALSE) 
  dat <- read.table(file,header=T,stringsAsFactors=FALSE)
  if(debug) print(dat)
  combinations <- expand.grid( dividend=dat$name, divisor=dat$name )
  values <- NULL
  
  pdgcol <- 2
  flagcol <- 4
  isocol <- 6
  
  columns <- c(pdgcol,flagcol,isocol)
  types <- c("PDG","FLAG","ISOSYM")
  
  savename <- "phys_values.csv"
  if(all) savename <- sprintf("all_%s",savename) 
  
  for( i in 1:length(combinations[,1]) ) {
    # when the dividend and divisor are the same, we transfer the value of the quantitiy to the table in physical units
    if( combinations[i,]$dividend == combinations[i,]$divisor ) {
      # prefer isospin symmetric, then PDG, then FLAG
      val.idx <- which( dat$name == combinations[i,]$dividend )
      if( !is.na( dat[val.idx,]$iso ) ) {
        values <- rbind(values,data.frame( name=dat[val.idx,]$name, val=dat[val.idx,]$iso, dval=dat[val.idx,]$diso, type="ISOSYM" ))
      } else if( !is.na( dat[val.idx,]$pdg ) ) {
        values <- rbind(values,data.frame( name=dat[val.idx,]$name, val=dat[val.idx,]$pdg, dval=dat[val.idx,]$dpdg, type="PDG" ))
      } else {
        values <- rbind(values,data.frame( name=dat[val.idx,]$name, val=dat[val.idx,]$flag, dval=dat[val.idx,]$dflag, type="FLAG" ))
      } 
    } else {
    # otherwise compute ratios with error propagation 
      name <- sprintf("%s_ov_%s", combinations[i,]$dividend, combinations[i,]$divisor )
    
      dividend.idx <- which( dat$name == combinations[i,]$dividend )
      divisor.idx <- which( dat$name == combinations[i,]$divisor )
      
      # extract data from dat
      dividend <- NULL
      divisor <- NULL
      
      # compute all possible ratios taking data from all possible sources
      if( all ) {
        for( i in seq(1,length(types))) {
          dividend <- NULL
          if(!is.na(dat[dividend.idx,columns[i]])) {
            dividend <- data.frame( val=dat[dividend.idx,columns[i]], dval=dat[dividend.idx,(columns[i]+1)], type=types[i] )
          } else {
            next()
          }
          for( j in seq(1,length(types))) {
            divisor <- NULL
            if(!is.na(dat[divisor.idx,columns[j]])) {
              divisor <- data.frame( val=dat[divisor.idx,columns[j]], dval=dat[divisor.idx,(columns[j]+1)], type=types[j] )
            } else {
              next()
            }        
            if( divisor$type == dividend$type ) {
              type <- divisor$type
            } else {
              type <- sprintf("%s/%s",dividend$type,divisor$type)
            }
            values <- rbind(values,c(compute_ratio(name=name,dividend=dividend,divisor=divisor),type=type))
          }
        }
      } else {
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
        values <- rbind(values,c(compute_ratio(name=name,dividend=dividend,divisor=divisor),type=type))
      }
    } # if(dividend==divisor)
  }
  write.csv(values,file=savename)
  return(as.data.frame(values))
}
