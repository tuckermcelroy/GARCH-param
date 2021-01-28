lik.tarch <- function(psi,data)
{
	#################################
	#	lik.tarch by Tucker McElroy
	#
	#	evaluates t-ARCH(p) lik using euclidean parametrization, p > 0
	#	input: psi is a p+2 vector of reals, data is T-vector of data,
	#		last parameter corresponds to degrees of freedom
	#	output: value of -2 times log t lik
	###############################################

	p <- length(psi)-2
	T <- length(data)
	arch <- psi2arch(psi[1:(p+1)])
	df <- 2 + exp(psi[p+2])

	lik <- 0
	for(t in (p+1):T)
	{
		arch.sd <- sqrt(arch[1] + sum(arch[-1]*data[(t-1):(t-p)]^2))
		lik <- lik + (-2)*log(dt(data[t]/arch.sd,df=df)/arch.sd) 
	}
	return(lik)
}
