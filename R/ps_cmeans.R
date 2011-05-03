###################################################
#		Esimate Cmeans using our methods
###################################################
Est.Depth <- function(n, iter=5)
{
	SMALL.VAL <- 1e-8
	cmeans <- colSums(n) / sum(n)
	keep <- NULL
	
	for (i in 1 : iter)
	{
		n0 <- rowSums(n) %*% t(cmeans)
		prop <- rowSums((n - n0) ^ 2 / (n0 + SMALL.VAL))
		qs <- quantile(prop, c(0.25, 0.75))
		keep <- (prop >= qs[1]) & (prop <= qs[2])
	
		cmeans <- colMeans(n[keep, ])
		cmeans <- cmeans / sum(cmeans)
	}
	
	return(list(cmeans=cmeans, keep=keep))
}
