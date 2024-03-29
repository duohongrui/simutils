% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEGs_manipulation.R
\name{calculate_DEGs_properties}
\alias{calculate_DEGs_properties}
\title{Summarize the Ability of Simulating DEGs}
\usage{
calculate_DEGs_properties(
  count_matrix,
  group,
  DEGs,
  DEA_method,
  model_method,
  verbose = TRUE
)
}
\arguments{
\item{count_matrix}{A matrix or a list with names of gene expression data.}

\item{group}{A vector of characters which indicate which group a cell belongs to. A list is also supported when many matrices of data is input in \code{count_matrix} parameter.}

\item{DEGs}{A list of DEGs with the names of \code{xxxvsxxx}. Note that the names of DEGs are in the rownames of the matrix or the dataframe and the names of \code{xxx} is in the \code{group} characters. If you have input the lists of \code{group} and \code{count_matrix}, the every sub-list in DEGs list should match with every data in \code{count_matrix} and \code{group} at the same index or position.}

\item{DEA_method}{The DEA method to get the DEGs. Choices: edgeRQLF, edgeRQLFDetRate, MASTcpmDetRate, MASTtpmDetRate, MASTcpm, MASTtpm, limmatrend, limmavoom, ttest and wilcox.}

\item{model_method}{The method to establish the model. SVM, Decision tree or RF (Random Forest).}

\item{verbose}{Whether the process massages are returned.}
}
\value{
A list
}
\description{
Summarize the Ability of Simulating DEGs
}
