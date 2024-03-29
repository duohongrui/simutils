% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method_definition.R
\name{method_definition}
\alias{method_definition}
\title{Method information and introduction}
\usage{
method_definition(
  method,
  programming = NULL,
  url = NULL,
  authors = list(),
  manuscript = list(),
  description = NULL,
  vignette = NULL
)
}
\arguments{
\item{method}{Method name.}

\item{programming}{Programming language.}

\item{url}{Relative website URL.}

\item{authors}{Methods creator.}

\item{manuscript}{Manuscript information.}

\item{description}{Descriptions to the method.}

\item{vignette}{The URL of method vignette}
}
\value{
A list.
}
\description{
Method information and introduction
}
\examples{
a <- method_definition(
method = "A",
programming = "Python",
url = "www.xxxx.com",
authors = authors_definition(first = "Hongrui",
                             last = "Duo",
                             email = "duohongrui@cqnu.edu.cn",
                             github = "https://github.com/duohongrui",
                             orcid = "0000-0001-8683-015X"),
manuscript = manuscript_definition(title = "xxxxx",
                                   doi = "xxxxxx",
                                   journal = "xxxxx",
                                   date = "2022",
                                   peer_review = TRUE),
description = "xxxxxx")
}
\seealso{
\code{\link[=authors_definition]{authors_definition()}}, \code{\link[=manuscript_definition]{manuscript_definition()}}
}
