\name{PoissonSeq-package}
\alias{PoissonSeq-package}
\alias{PoissonSeq}
\docType{package}
\title{
 Significance analysis of sequencing data based on a Poisson log linear model
}
\description{
 This package implements all the methods described in the paper
 
 Li J, Witten DM, Johnstone I, Tibshirani R (2011).
 Normalization, testing, and false discovery rate estimation for RNA-sequencing data.
 Submitted.
 
 Especially, it implements a new method for estimating sequencing depths according to the
 Poisson goodness-of-fit statistic, and estimate the false discovery rates using a modified
 version of the permutation plug-in method.
}
\details{
\tabular{ll}{
Package: \tab PoissonSeq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-05-3\cr
License: \tab GPL (version 2 or newer)\cr
LazyData: \tab yes\cr
}
 This package mainly includes the following functions:\cr
 \code{Est.Depth}: Estimate the sequencing depths.\cr
 \code{Trans.Data}: Transform the data matrix to make it more like Poisson-distributed.\cr
 \code{Permu.Ind}: Get the permutation index matrix.\cr
 \code{PS.FDR}: Estimate the false discovery rates for RNA-Seq data.\cr
 \code{Score.Stat}: Calculate the score statistics and the permuted score statistics.\cr
 \code{PS.Main}: Estimate the false discovery rates for RNA-Seq data.\cr
}
\author{
Jun Li

Maintainer: Jun Li <junli07@stanford.edu>
}
\references{
 Li J, Witten DM, Johnstone I, Tibshirani R (2011).
 Normalization, testing, and false discovery rate estimation for RNA-sequencing data.
 Submitted.
}
\examples{
 data(dat)	# load the example data
  
 # for most cases, just a single sentence below is all what you need to do
 res <- PS.Main(dat)
 
 # if you want to do step by step, you can use the following sentences.
 # They are what PS.Main does.
 
 # set the parameters
 para <- list(trans=TRUE, npermu=100, nvals=nrow(dat$n), div=10, pow.file='pow.txt')

 seq.depth <- Est.Depth(dat$n)	# estimate sequencing depth
 dat$n <- Trans.Data(n=dat$n)	# transform the data
 ps.obj <- Score.Stat(dat=dat, para=para)	# calculate the score statistics
 res <- PS.FDR(ps.obj=ps.obj, nvals=para$nvals)	# estimate FDRs
}
\keyword{ package }
\keyword{ models }
