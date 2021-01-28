lik.arch <- function(psi,data,df)
{
	#################################
	#	lik.arch by Tucker McElroy
	#
	#	evaluates ARCH(p) lik using euclidean parametrization, p > 0
	#	input: psi is a p+1 vector of reals, data is T-vector of data
	#	output: value of -2 times log Gaussian/t lik
	###############################################

	p <- length(psi)-1
	T <- length(data)
	arch <- psi2arch(psi)

	lik <- 0
	for(t in (p+1):T)
	{
		arch.sd <- sqrt(arch[1] + sum(arch[-1]*data[(t-1):(t-p)]^2))
		if(df==Inf) { lik <- lik + (-2)*log(dnorm(data[t],sd=arch.sd)) } else {
			lik <- lik + (-2)*log(dt(data[t]/arch.sd,df=df)/arch.sd) }
	}
	return(lik)
}
