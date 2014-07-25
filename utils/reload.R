# convenience functions for re/un-loading a package including any dynamic
# libraries belonging to the package

unload <- function(pkgName) {
  detachName <- sprintf("package:%s",pkgName)
  if( any(search() == detachName) ) {
    detach(name=detachName, unload = TRUE, character.only=TRUE)
    library.dynam.unload(pkgName, system.file(package = pkgName))
  } else {
    cat(sprintf("No package '%s' in search()\n",pkgName))
  }
}

reload <- function(pkgName) {
  unload(pkgName)
  require(package=pkgName,character.only=TRUE)                                                                                                               
}

