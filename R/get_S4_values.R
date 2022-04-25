#' Get Values From S4 Object
#'
#' @param S4Object A S4object
#'
#' @return A list
#' @export
#'
get_S4_values <- function(S4Object){
  purrr::map(slotNames(S4Object), function(x){
    exp <- paste0(quote(S4Object),"@", x)
    eval(parse(text = exp))
  })
}
