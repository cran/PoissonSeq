\name{PS.Est.Depth}
\alias{PS.Est.Depth}
\title{ estimate the sequencing depths }
\description{
 Given a data matrix, this function will estimate the sequencing depths
 on the basis of the Poisson goodness-of-fit statistic.
 This estimate is applicable to data with any types of outcome, as it
 estimates under the null hypothesis. \cr
 A more detailed instruction as well as sample data is available at \cr
 \code{http://www.stanford.edu/~junli07/research.html}.
}
\usage{
 PS.Est.Depth(n, iter=5, ct.sum=5, ct.mean=0.5)
}
\arguments{
 \item{n}{ The data matrix. The rows are counts for a gene, and the columns are counts from an experiment.}
 \item{iter}{ Number of iterations used. Default value: 5. The default value is usually a good choice.}
 \item{ct.sum}{ if the total number of reads of a gene across all experiments <= ct.sum, this gene will not be considered for estimating sequencing depth. Default value: 5. }
 \item{ct.mean}{if the mean number of reads of a gene across all experiments <= ct.mean, this gene will not be considered for estimating sequencing depth. Default value: 0.5.}
}
\details{
 The value in the data matrix does not need to be integers.
}
\value{
 estimated sequencing depth. a vector. their product is 1.
}
\references{
 Li J, Witten DM, Johnstone I, Tibshirani R (2011).
 Normalization, testing, and false discovery rate estimation for RNA-sequencing data.
 To appear, Biostatistics.
}
\author{
 Jun Li.
}
\examples{
 data(dat)
 seq.depth <- PS.Est.Depth(dat$n)
}
\keyword{ models }
\keyword{ iteration }
