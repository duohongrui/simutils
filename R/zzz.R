.onAttach <- function(libname, pkgname) {
  suppressPackageStartupMessages({
    require(powsimR)
    require(rtracklayer)
  })
}
