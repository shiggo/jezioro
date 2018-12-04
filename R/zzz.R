.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("This is jezioro v", packageVersion("jezioro"), sep=""))
}
