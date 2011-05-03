\name{Permu.Ind}
\alias{Permu.Ind}
\title{ Get the permuation index matrix }
\description{
 Generate a matrix of permuted indexes.
 Each row is a vector of permuted indexes.
}
\usage{
 Permu.Ind(y, npermu, pair=F)
}
\arguments{
 \item{y}{ The outcome vector. The values are real numbers for quantitative outcomes. For K-class data, values must be 1, 2, ..., K, with value k indicating Class K. For 2-class data, values must be 1 or 2. Especially, for paired 2-class data, the vector must be 1, 2, 1, 2, 1, 2, ... Each continuous pair of 1 and 2 corresponds to a pair of samples.}
 \item{npermu}{ Number of permutations.}
 \item{pair}{ Whether the data is paired or not. Can be true only for 2-class data.}
}
\details{
 When the sample size is small, the number of possible permutations can be smaller than npermu.
 In this case, this function will generate all possible permuations instead. So the number of rows in the returned matrix is less than npermu.
}
\value{
 a matrix of permuted indexes. Each row is a vector of permuted indexes.
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
 pind <- Permu.Ind(dat$y, 200, dat$pair)
}
\keyword{ models }
