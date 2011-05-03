\name{Est.Depth}
\alias{Est.Depth}
\title{ estimate the sequencing depths }
\description{
 Given a data matrix, this function will estimate the sequencing depths
 on the basis of the Poisson goodness-of-fit statistic.
 This estimate is applicable to data with any types of outcome, as it
 estimates under the null hypothesis. The outcome vector is not used 
 as an input parameter.
}
\usage{
 Est.Depth(n, iter=5)
}
\arguments{
 \item{n}{ The data matrix. The rows are counts for a gene, and the columns are counts from an experiment.}
 \item{iter}{ Number of iterations used. The default value is usually a good choice in most cases.}
}
\details{
 The value in the data matrix does not need to be integers.
}
\value{
 a list.
 \item{cmeans}{The estimated sequencing depths. A vector with elements summed to 1.}
 \item{keep}{The genes in set \code{S}: the half of null genes.}
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
 seq.depth <- Est.Depth(dat$n)$cmeans
}
\keyword{ models }
\keyword{ iteration }
