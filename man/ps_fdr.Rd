\name{PS.FDR}
\alias{PS.FDR}
\title{ estimate FDR }
\description{
 Given the score statstics and the permuted score statistics, 
 this function estimate the false discovery rate using a modified
 version of permutation plug-in method.
 It also reports the the estimated permuation-based p-values,
 and the estimated percentage of null genes.
}
\usage{
 PS.FDR(ps.obj, nvals)
}
\arguments{
 \item{ps.obj}{ A list containing the score statistics and the permuted score statistics. 
 This list is often created by calling Score.Stat(dat, para).}
 \item{nvals}{ Number of cutoffs used. }
}
\details{
 Often, nvals is set to be the number of genes in the data.
 A larger nvals will not significantly increase the computational time of our algorithm.
 If nvals is set to be larger than the number of genes, our algorithm
 will set it as the number of genes instead.
}
\value{
 a list.
 \item{dels}{The cutoffs used. The function will choose this cutoff vector automatically according to the value of nvals.}
 \item{tt}{The score statistics of the genes.}
 \item{sig.ord}{The order of the score statistics, from most significant to most insignificant. For example, if sig.ord=(50, 10, ...), then the most significant gene is the 50th, the second most significant gene is the 10th, etc.}
 \item{nc}{A vector of number of significant genes called. If nvals is set to be the number of genes, then nc = 1 : (number of genes). Otherwise, nc is roughly an arithmetical series with length nvals, and maximum value number of genes.}
 \item{fdr}{A vector of estimated false discovery rates. fdr is of the same length as nc. fdr[i] is the estimated false discovery rate when nc[i] genes are called.}
 \item{pi0}{Estimated proportion of null genes.}
 \item{pval}{Permutation-based p-values.}
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
 data(dat)	# read in the data
 # set the parameters
 para <- list(trans=TRUE, npermu=100, nvals=nrow(dat$n), div=10, pow.file='pow.txt')

 ps.obj <- Score.Stat(dat, para=para)	# calculate the score statistics
 res <- PS.FDR(ps.obj, nrow(dat$n))	# estimate FDRs
}
\keyword{ models }
\keyword{ htest }
