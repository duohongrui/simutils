% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters.R
\name{param_dataframe}
\alias{param_dataframe}
\title{Define Dataframe Object Parameters}
\usage{
param_dataframe(
  id = NULL,
  type = c("data.frame", "tibble"),
  description = NULL,
  force = FALSE,
  process = "simulation",
  function_name
)
}
\arguments{
\item{id}{A parameter name. Default is NULL, which means no dataframe object input.}

\item{type}{The type of parameter. Must be a data.frame or tibble.}

\item{description}{The description information of this parameter.}

\item{force}{Logical, whether this parameter must be needed or not.}

\item{process}{Two options, simulation or estimation. Which process does this
parameter involved in.}

\item{function_name}{The function name where the parameter exists.}
}
\value{
A list.
}
\description{
Define Dataframe Object Parameters
}
\examples{
network_input <- param_dataframe(
id = 'network',
process = "simulation",
description = "The network information is used to show the relationship of nodes."
)
}
