sim.arch <- function(arch,T,init,df)
{
	#################################
	#	sim.arch by Tucker McElroy
	#
	#	simulates ARCH(p) of length T, p > 0
	#	input: coefs a0, a1, ..., ap of ARCH(p), T is desired length,
	#		init is p vector of initial values
	#	output: simulation with Gaussian or t errors (df = degrees freedom)
	###############################################

	p <- length(arch)-1
	data <- init
	for(t in (p+1):T)
	{
		if(df==Inf) { z <- rnorm(1) } else { z <- rt(1,df=df) }
		datum <- z*sqrt(arch[1] + sum(arch[-1]*data[(t-1):(t-p)]^2))
		data <- c(data,datum)
	}
	return(data)
}
