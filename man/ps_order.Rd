\name{Trans.Data}
\alias{Trans.Data}
\title{ Transform the data matrix. }
\description{
 Some sequencing data are largely overdispersed. We use this function to transform
 the data matrix so that they look more like Poisson. This function will determine
 the proper power for transformation according to the data automatically.
}
\usage{
 Trans.Data(n, div=10, pow.file="")
}
\arguments{
 \item{n}{ The data matrix. The rows are counts for a gene, and the columns are counts from an experiment.}
 \item{div}{ Number of groups that the data are divided into. Often the default value is a good choice.}
 \item{pow.file}{ The files that we want to store the values of pairs of log(N_{\\cdot j}) and \\theta_{-1}. This file might be used to check the performance of our transformation. Often, we expect that their relationship is monotone and smooth. The default value is ""; in this case, no file is generated.}
}
\details{
 The genes are divided into \code{div} groups according to their average expression level. In each group, a proper power is estimated. Using the fitted natural cubic spline between log(N_{\\cdot j}) and \\theta_{-1}, we can predict a proper value of theta for each gene. Please refer to the paper for details.
}
\value{
 A new data matrix. This data matrix is used for estimating false discovery rates.
}
\references{
 Li J, Witten DM, Johnstone I, Tibshirani R (2011).
 Normalization, testing, and false discovery rate estimation for RNA-sequencing data.
 Submitted.
}
\author{
 Jun Li.
}
\examples{
 data(dat)
 n.trans <- Trans.Data(dat$n)
}
\keyword{ models }
