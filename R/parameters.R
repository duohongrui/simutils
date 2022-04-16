#' Define Numeric Parameters
#'
#' @param id A parameter name.
#' @param type The type of parameter. Must be a numeric.
#' @param default The default value.
#' @param lower The lowest range of this parameter.
#' @param upper The highest range of this parameter.
#' @param border Boolean. Whether the default value can get to the range border.
#' @param description The description information of this parameter.
#'
#' @return A list.
#' @export
#'
#' @examples
#' probability <- param_numeric(id = "probability",
#'                              default = 0.5,
#'                              lower = 0,
#'                              upper = 1,
#'                              border = TRUE,
#'                              description = "A numeric value to specify the probability to randomly select the right number")
param_numeric <- function(id,
                          type = "numeric",
                          default,
                          lower = NULL,
                          upper = NULL,
                          border = TRUE,
                          description){

  assertthat::assert_that(is.character(id),
                          is.numeric(default),
                          is.character(description))
  if(type != "numeric"){
    stop("This is not fit for other types except for numeric values")
  }

  if(is.null(lower)){
    cat("The parameter lower is not defined and will be set to -Inf")
    lower <- -Inf
  }
  if(is.null(upper)){
    cat("The parameter upper is not defined and will be set to Inf")
    upper <- Inf
  }
  assertthat::assert_that(is.numeric(lower),
                          is.numeric(upper))
  if(lower >= upper){
    stop('Lower value must be smaller than upper value')
  }

  if(default < lower | default > upper) stop("The defaul value is not between the lower and the upper")

  if(!border){
    if(default == lower | default == upper) stop("The defaul value can not euqal to the lower or the upper\nPlease reset!")
  }

  tibble::lst(id,
              type,
              default,
              lower,
              upper,
              border,
              description)

}





#' Define Integer Parameters
#'
#' @param id A parameter name.
#' @param type The type of parameter. Must be an integer.
#' @param default The default value.
#' @param lower The lowest range of this parameter.
#' @param upper The highest range of this parameter.
#' @param border Boolean. Whether the default value can get to the range border.
#' @param description The description information of this parameter.
#'
#' @return A list.
#' @export
#'
#' @examples
#' group_num <- param_integer(id = "group_num",
#'                            default = 1L,
#'                            lower = 1L,
#'                            upper = 50L,
#'                            border = TRUE,
#'                            description = "An integer to specify how many groups to be simulated")
param_integer <- function(id,
                          type = "integer",
                          default,
                          lower = NULL,
                          upper = NULL,
                          border = TRUE,
                          description){

  assertthat::assert_that(is.character(id),
                          is.integer(default),
                          is.character(description))
  if(type != "integer"){
    stop("This is not fit for other types except for integers")
  }

  if(is.null(lower)){
    cat("The parameter lower is not defined and will be set to -Inf")
    lower <- -Inf
  }
  if(is.null(upper)){
    cat("The parameter upper is not defined and will be set to Inf")
    upper <- Inf
  }
  assertthat::assert_that(is.integer(lower),
                          is.integer(upper))
  if(lower >= upper){
    stop('Lower value must be smaller than upper value')
  }

  if(default < lower | default > upper) stop("The defaul value is not between the lower and the upper")

  if(!border){
    if(default == lower | default == upper) stop("The defaul value can not euqal to the lower or the upper\nPlease reset!")
  }

  tibble::lst(id,
              type,
              default,
              lower,
              upper,
              border,
              description)

}





#' Define Logical Parameters
#'
#' @param id A parameter name.
#' @param type The type of parameter. Must be a logical value.
#' @param default The default value. TRUE or FALSE.
#' @param description The description information of this parameter.
#'
#' @return A list.
#' @export
#'
#' @examples
#' verbose <- param_Boolean(id = "verbose",
#'                          default = TRUE,
#'                          description = "Whether to return the information during the process")
param_Boolean <- function(id,
                          type = "logical",
                          default,
                          description){

  assertthat::assert_that(is.character(id),
                          is.logical(default),
                          is.character(description))
  if(type != "logical"){
    stop("This is not fit for other types except for logical values")
  }

  tibble::lst(id,
              type,
              default,
              description)

}





#' Define Character Parameters
#'
#' @param id A parameter name.
#' @param type The type of parameter. Must be characters.
#' @param default The default values. Length to one or more.
#' @param alternatives Options that can be chosen.
#' @param description The description information of this parameter.
#'
#' @return A list.
#' @export
#'
#' @examples
#' dimension_method <- param_character(id = 'dimension_method',
#'                                     default = "UMAP",
#'                                     alternatives = c("UMAP", "TSNE", "PCA"),
#'                                     description = " A character string specifying the algorithm to use for dimensionality reduction.")
param_character <- function(id,
                            type = "character",
                            default,
                            alternatives,
                            description){

  assertthat::assert_that(is.character(id),
                          is.character(default),
                          is.character(description))
  if(type != "character"){
    stop("This is not fit for other types except for characters")
  }

  if(!default %in% alternatives){
    stop("The default value is not in your alternatives. Please check the spelling and input")
  }

  tibble::lst(id,
              type,
              default,
              alternatives,
              description)

}


