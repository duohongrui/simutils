#' Method information and introduction
#'
#' @param method Method name.
#' @param programming Programming language.
#' @param url Relative website URL.
#' @param authors Methods creator.
#' @param description Descriptions to the method.
#' @param manuscript Manuscript information.
#'
#' @return A list.
#' @export
#' @seealso [authors_definition()], [manuscript_definition()]
#'
#' @examples
#' a <- method_definition(
#' method = "A",
#' programming = "Python",
#' url = "www.xxxx.com",
#' authors = authors_definition(first = "Duo",
#'                              last = "hongrui",
#'                              email = "duohongrui@cqnu.edu.cn",
#'                              github = "https://github.com/duohongrui",
#'                              orcid = "0000-0001-8683-015X"),
#' manuscript = manuscript_definition(title = "xxxxx",
#'                                    doi = "xxxxxx",
#'                                    journal = "xxxxx",
#'                                    date = "2022",
#'                                    peer_review = TRUE),
#' description = "xxxxxx")
method_definition <- function(
  method,
  programming = NULL,
  url = NULL,
  authors = list(),
  manuscript = list(),
  description = NULL){
  as.list(environment())
}


#' Authors information
#'
#' @param first First name.
#' @param last Last name.
#' @param email Email address.
#' @param github Github URL.
#' @param orcid Orcid URL.
#'
#' @return A list.
#' @export
#'
authors_definition <- function(
  first,
  last,
  email = NULL,
  github = NULL,
  orcid = NULL){
  as.list(environment())
}


#' Manuscript Information
#'
#' @param title Manuscript title.
#' @param doi Doi number.
#' @param journal Journal name.
#' @param date Publication date.
#' @param peer_review Has been peer reviewed?
#'
#' @return A list.
#' @export
#'
manuscript_definition <- function(
  title,
  doi,
  journal = NULL,
  date = NULL,
  peer_review = TRUE){
  as.list(environment())
}
