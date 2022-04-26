#' Get Values From S4 Object
#'
#' This function is used to get the slot values from S4 objects conveniently. In
#' addition, a list only containing one or more S4 subjects is supported.
#' @param S4Object A S4 object or a list only containing S4 objects.
#'
#' @return A list
#' @importFrom stats setNames
#' @importFrom methods slotNames
#' @export
#'
get_S4_values <- function(S4Object){
  if(is.list(S4Object)){
    # Set list names if they are NULL.
    if(is.null(names(S4Object))){
      stats::setNames(S4Object, paste0("S4_object_", base::seq_len(length(S4Object))))
    }
    # A new variable
    S4_list <- S4Object
    # Get the S4 objects from the list and run get_S4_values function.
    result <- purrr::map(names(S4_list), function(i){
      sub_S4 <- S4_list[[i]]
      get_S4_values(sub_S4)
    })
    # Return
    result
  }
  # If S4Object is a S4 object and then get the slot values from it.
  if(!is.list(S4Object) & base::isS4(S4Object)){
    purrr::map(methods::slotNames(S4Object), function(x){
      exp <- paste0(quote(S4Object),"@", x)
      eval(parse(text = exp))
    })%>% stats::setNames(methods::slotNames(S4Object))
  }
}
