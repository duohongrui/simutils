#' A Random Seed
#'
#' Get a random seed
#'
#' @export
#' @examples
#' seed <- random_seed()
random_seed <- function(){
  seed <- .Machine$integer.max
}
