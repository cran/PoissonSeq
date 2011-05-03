\name{Score.Stat}
\alias{Score.Stat}
\title{ calculate the score statistics and the permuted score statistics }
\description{
 Calculate the score statistic given the data matrix and the outcome vector.
}
\usage{
 Score.Stat(dat, para)
}
\arguments{
 \item{dat}{ The data list, including the data matrix and the outcome vector.}
 \item{para}{ The parameter list. Particularly, it sets the number of permutations used. }
}
\details{
 This function will calculate both the score statistics for the original data matrix,
 and the score statistics for the permuted data matrix.
 It deals with data with a two-class, multiple-class or quantitative outcome.
}
\value{
 a list.
 \item{tt}{The score statistics for the original data matrix. Unsigned.}
 \item{tt.signed}{The signed score statistics of for the original data matrix. This vector will not be used by our program for estimating false discovery rates.}
 \item{ttstar0}{The score statistics for the permuted data matrix. Unsigned.}
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
 para <- list(npermu=100)
 ps.obj <- Score.Stat(dat, para=para)
}
\keyword{ models }
\keyword{ univar }
