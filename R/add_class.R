#' Add A New Class
#'
#' @param x An object.
#' @param define_class A new self-defined class
#' @export
#'
#' @examples
#' a <- matrix(rnorm(1000), nrow = 100)
#' a <- add_class(a, 'my_matrix')
add_class <- function(x, define_class){
  base::class(x) <- c(base::class(x), define_class)
  x
}
