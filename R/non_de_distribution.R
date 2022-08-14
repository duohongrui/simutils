#' Test if the p-value is Uniformly Distributed
#'
#' @param PValue The p-value vector.
#' @export
#' @examples
#' pvalue <- runif(1000)
#' result <- test_uni_distribution(pvalue)
test_uni_distribution <- function(
  PValue
){
  if(!requireNamespace("spgs", quietly = TRUE)){
    stop("Please install \"spgs\" package using \"install.packages('spgs')\" command.")
  }
  result <- spgs::chisq.unif.test(PValue, bins = 20)
  if(result$p.value < 0.05){
    score <- 0
  }else{
    score <- 1
  }
  return(dplyr::lst(result,
                    score))
}
