\name{dat}
\alias{dat}
\docType{data}
\title{ A toy RNA-Seq data. }
\description{
 This toy data contains 1000 genes, and 6 samples in each of the two classes.
 The sequencing depths of the 12 samples are rep(c(1, 2), 6).
 The first 80 genes are set to be overexpressed in Class 2, and the next 20 genes are set to be underexpressed in Class 2.
 The other 900 genes are null.
 This is a Poisson-distributed data.
}
\usage{data(dat)}
\format{
  A list.
  \describe{
		\item{\code{n}}{the count matrix}
		\item{\code{y}}{the outcome vector}
		\item{\code{type}}{the outcome type}
		\item{\code{pair}}{whether the samples are paired or not}
  }
}
\keyword{datasets}
