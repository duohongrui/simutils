.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to simutils")
  suppressPackageStartupMessages({
    require(powsimR)
    require(rtracklayer)
  })
}
