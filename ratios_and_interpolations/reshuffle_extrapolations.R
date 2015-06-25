# function designed to reshuffle a list of extrapolation results and their 
# errors into a particular ordering
reshuffle_ext_table <- function(ext) {
  rval <- ext
  ordering <- c("$M_K/M_\\pi$", "$M_{D}/M_\\pi$", "$M_{D_s}/M_\\pi$", "$M_{D_s}/M_K$", "$M_{D}/M_K$", "$M_{D_s}/M_D$", 
                "$M_K/f_K$", "$M_D/f_D$", "$M_{D_s}/f_{D_s}$",
                "$f_K/f_\\pi$", "$f_D/f_\\pi$", "$f_{D_s}/f_\\pi$", "$f_D/f_K$", "$f_{D_s}/f_K$", "$f_{D_s}/f_D$",
                "$aM_K$", "$aM_D$", "$aM_{D_s}$",
                "$af_K$", "$af_D$", "$af_{D_s}$")

  if(length(ordering)!=length(ext)){
    stop("reshuffule_extrapolations: ordering vector and extrapolations list do not have the same length!")
  }

  for( i in 1:length(ordering) ){
    for( j in 1:length(ext$labels) ){
      if( ext$labels[j] != ordering[i] ) next
      for( k in 1:length(rval) ) {
        if(class(rval[[k]])=="numeric" || class(rval[[k]])=="character"){
          rval[[k]][i] <- ext[[k]][j]
        } else {
          rval[[k]][i,] <- ext[[k]][j,]
        }
      }
    }
  }
  rval
}
  

