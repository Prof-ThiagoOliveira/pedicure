% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all.R
\name{amat}
\alias{amat}
\title{Relationship matrix}
\usage{
amat(ped, keep = rep(TRUE, nrow(ped)))
}
\arguments{
\item{ped}{A data frame with (at least) three columns that correspond to the
individual, male parent and female parent, respectively. Founders
use 0 (zero) or \code{NA} in the parental columns.}

\item{keep}{A logical vector identifying the rows of \code{ped} to retain.}
}
\value{
A \code{matrix} object containing the numerator
relationship matrix.
}
\description{
Calculate the numerator relationship matrix from a pedigree.
}
\section{Side-effects}{

A file \code{"A.txt"} is created
}

