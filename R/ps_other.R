############################################################
#		Filter genes with too small counts
############################################################
PS.Filter <- function(dat, ct.sum=5, ct.mean=0.5)
{
	if (is.null(dat$gname))
	{
		dat$gname <- 1 : nrow(dat$n)
	}
	
	keep <- (rowMeans(dat$n) > ct.mean) & (rowSums(dat$n) > ct.sum)
	cat(length(keep) - sum(keep), "genes has been filtered because they contains too small number of reads across the experiments.")
	
	dat$n <- dat$n[keep, ]
	dat$gname <- dat$gname[keep]
	
	return(dat)
}

############################################################
#		Summarize results
############################################################
PS.Sum <- function(dat, fdr.res, seq.depth, trans)
{
	### sort nc and fdr
	tt <- sort(abs(fdr.res$tt), decreasing=T)
	nc <- 1 : length(fdr.res$nc)
	ord.fdr <- fdr.res$fdr[order(fdr.res$nc)]
	pval <- sort(fdr.res$pval)
	
	### make sure the fdr is monotone increasing
	fdr <- rep(min(c(ord.fdr[length(ord.fdr)], 1)), length(ord.fdr))
	for (i in (length(fdr) - 1) : 1)
	{
		fdr[i] <- min(c(fdr[i + 1], ord.fdr[i], 1))
	}
	
	### other elements
	gname <- dat$gname[fdr.res$sig.ord]
	
	res.table <- data.frame(nc=nc, gname=gname, tt=tt, pval=pval, fdr=fdr)
	
	### get fold change for two class data
	if (dat$type == "twoclass")
	{
		if (trans)
		{
			exp.1 <- rowSums(dat$n.ori[, dat$y==1]) / sum(seq.depth[dat$y==1])
			exp.2 <- rowSums(dat$n.ori[, dat$y==2]) / sum(seq.depth[dat$y==2])
		} else
		{
			exp.1 <- rowSums(dat$n[, dat$y==1]) / sum(seq.depth[dat$y==1])
			exp.2 <- rowSums(dat$n[, dat$y==2]) / sum(seq.depth[dat$y==2])
		}
		log.fc <- log(exp.2 / exp.1)
		res.table <- cbind(res.table, log.fc=log.fc[fdr.res$sig.ord])
	}

	return(res.table)
}
