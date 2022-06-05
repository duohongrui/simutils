#' Fix Paths on Windows
#'
#' @param path A path under Windows
#'
#' @export
#'
fix_path <- function(path){
  stringr::str_replace_all(path, pattern = "\\\\", replacement = "/")
}
