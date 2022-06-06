#' A String of Time
#'
#' @export
#'
#' @examples
#' time_string <- time_string()
time_string <- function(){
  format(Sys.time(), format = "%Y%m%d%H%M%S")
}
