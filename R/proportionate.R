#' A Number Multiplies Proportions Under Limitations
#'
#' @param number A number
#' @param result_sum_strict The limitation of the sum of results
#' @param prop The proportions to be multiplied by the number
#' @param prop_sum_strict The limitation of the sum of proportions
#' @param digits The number of decimal places
#'
#' @export
#'
#' @examples
#' ## Save 1 decimal place
#' a <- proportionate(number = 355,
#'                    prop = c(0.2, 0.6, 0.15, 0.36),
#'                    digits = 1)
#' ## The sum of the proportions is 1
#' b <- proportionate(number = 355,
#'                    prop = c(0.2, 0.6, 0.15, 0.05),
#'                    prop_sum_strict = 1,
#'                    digits = 1)
#' ## Save 0 decimal place
#' c <- proportionate(number = 355,
#'                    prop = c(0.2, 0.6, 0.15, 0.05),
#'                    prop_sum_strict = 1,
#'                    digits = 0)
#' ## The sum of the results is 355
#' d <- proportionate(number = 355,
#'                    result_sum_strict = 355,
#'                    prop = c(0.2, 0.6, 0.15, 0.05),
#'                    prop_sum_strict = 1,
#'                    digits = 0)
proportionate <- function(
    number,
    result_sum_strict = NULL,
    prop,
    prop_sum_strict = NULL,
    digits = 0
){
  # Check
  if(!is.null(prop_sum_strict)){
    if(sum(prop) != prop_sum_strict) stop("The sum of proportions is not equal to the specific value")
  }
  # Assign
  len_prop <- length(prop)
  assign_num <- number * prop
  assign_num <- round(assign_num, digits = digits)
  if(!is.null(result_sum_strict)){
    if(sum(assign_num) != result_sum_strict){
      assign_num <- c(assign_num[1:(len_prop-1)],
                      result_sum_strict - sum(assign_num[1:(len_prop-1)]))
      assertthat::assert_that(sum(assign_num) == result_sum_strict)
    }
  }
  return(assign_num)
}
