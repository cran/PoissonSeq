\name{PS.Main}
\alias{PS.Main}
\title{ estimate the false discovery rates for RNA-Seq data. }
\description{
 This function is the main function of this package.
 Given the data matrix and the outcome vector, this function will return
 the estimated false discovery rates, the estimated sequencing depths,
 and the estimated permuation-based p-values.
}
\usage{
 PS.Main(dat, para=list())
}
\arguments{
 \item{dat}{ The input RNA-Seq data. It should have the following attributes: 
 
 (1) \code{n}: the data matrix
 
 (2) \code{y}: the outcome vector
 
 (3) \code{type}: 'twoclass', 'multiclass' or 'quant'
 
 (4) \code{pair}: paired data or not
 
 All the above four attributes must be given.}
 
 \item{para}{ A list of parameters. It can have the following attributes:
 
 (1) \code{trans}: whether do transformation beforehand or not. default value: TRUE
 
 (2) \code{npermu}: number of permuations. default value: 100
 
 (3) \code{nvals}: number of cutoffs. default value: number of genes.
 
 (4) \code{div}: the number of divisions of genes for estimating theta. default value: 10
 
 (5) \code{pow.file}: the file to store the power transform curve (mean(log(mu)) ~ 1/theta). default value: 'pow.txt'
 
 All the above attributes are optional.}
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
 \item{seq.depth}{Estimated sequencing depths.}
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
 res <- PS.Main(dat)
}
\keyword{ models }
