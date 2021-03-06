% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all.R
\name{ped.dot}
\alias{ped.dot}
\title{Dot file representing a pedigree tree.}
\usage{
ped.dot(
  ped,
  keep = rep(TRUE, nrow(ped)),
  dotfile = "ped",
  numeric = TRUE,
  url = "",
  height = 0.5,
  width = 0.75,
  rotate = 0
)
}
\arguments{
\item{ped}{A data frame with (at least) three columns that correspond to the
individual, male parent and female parent, respectively. Founders
use 0 (zero) or \code{NA} in the parental columns.}

\item{keep}{A logical vector identifying the rows of \code{ped} to retain. The
default is \code{rep(TRUE, nrow(ped))}.}

\item{dotfile}{The output dot file primary name, default is \code{"ped"}; the
suffix \code{".dot"} is appended to the file name. Also used to
label the pedigree tree.}

\item{numeric}{Convert \code{ped} to integer values for graphing. The default is
\code{TRUE}.}

\item{url}{If not \code{""} (the default) then \code{url} is linked to the resulting graph.}

\item{height}{Node height.}

\item{width}{Node width}

\item{rotate}{If \code{rotate=90} landscape mode is selected; the default is
\code{0}.}
}
\description{
Exports the relationship structure in a pedigree as a directed
graph in a dot file suitable for plotting.
}
\details{
The resulting dot file can be edited prior to rendering with
\emph{Graphviz}, say, or conversion to a graphics file format with
the \emph{dot} application (part of the \emph{Graphviz} package). The
default draws founder individuals in rectangles, offspring in
ellipses, and offspring with a single parent in circles.
}
\section{Side effects}{


The \code{dotfile} is written to the working directory using
\code{sink()}.
}

