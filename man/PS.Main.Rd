\name{PS.Main}
\alias{PS.Main}
\title{ detecting differentially expressed genes from RNA-Seq data. }
\description{
 This function is the main function of this package.
 Given the data matrix and the outcome vector, this function returns
 the estimated permuation-based p-values, 
 the estimated permutation-based false discovery rates,
 et al.\cr
 A more detailed instruction as well as sample data is available at \cr
 \code{http://www.stanford.edu/~junli07/research.html}.
}
\usage{
 PS.Main(dat, para=list())
}
\arguments{
 \item{dat}{ The input RNA-Seq data. It \emph{must} have the following three attributes: 
 
 (1) \code{n}: the data matrix. Rows for genes, columns for experiments (samples).
 
 (2) \code{y}: the outcome vector
 
 (3) \code{type}: 'twoclass', 'multiclass' or 'quant'
 
 The following attributes are \emph{optional}. If not specified, the default values will be used.
 
 (4) \code{pair}: paired data or not. Default value: FALSE. Only take effect for twoclass data.
 
 (5) \code{gname}: gene names. Default value: 1 : nrow(n). That is, the i'th gene is named "i".}
 
 \item{para}{ A list of parameters. It can have the following attributes:
 
 (1) \code{trans}: to tranform the data using the order transformation or not to transform it. default value: TRUE
 
 (2) \code{npermu}: number of permuations. default value: 100
 
 (3) \code{seed}: random seed to generate the permutation indexes. default value: 10
 
 (4) \code{ct.sum}: if the total number of reads of a gene across all experiments <= ct.sum, this gene will not be considered for differential expression detection. Default value: 5.
 
 (5) \code{ct.mean}: if the mean number of reads of a gene across all experiments <= ct.mean, this gene will not be considered for differential expression detection. Default value: 0.5.

 (6) \code{div}: the number of divisions of genes for estimating theta. default value: 10
 
 (7) \code{pow.file}: the file to store the power transform curve (mean(log(mu)) ~ 1/theta). default value: 'pow.txt'
 
 All the above attributes are \emph{optional}.}
}
\value{
	a data frame (table) containing the following columns. Each row stands for a gene. The genes are sorted from the most significant to the most insignificant. 
 	\item{nc}{number of significant genes called. nc = 1 : (number of genes).}
 	\item{gname}{the sorted gene names.}
 	\item{tt}{The score statistics of the genes.}
 	\item{pval}{Permutation-based p-values of the genes.}
 	\item{fdr}{Estimated false discovery rate.}
 	\item{log.fc}{Estimated log fold change of the genes. Only available for twoclass outcomes.}
}
\references{
 Li J, Witten DM, Johnstone I, Tibshirani R (2012).
 Normalization, testing, and false discovery rate estimation for RNA-sequencing data.
 Biostatistics 13(3): 523-38.
}
\author{
 Jun Li.
}
\examples{
 data(dat)
 res <- PS.Main(dat)
}
\keyword{ models }
