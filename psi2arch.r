psi2arch <- function(psi)
{

	#################################
	#	psi2arch by Tucker McElroy
	#
	#	gives exact parametrization of ARCH(p)
	#	input: p+1 vector of reals
	#	output: coefs a0, a1, ..., ap of ARCH(p)
	###############################################

	p <- length(psi)-1
	a.0 <- exp(psi[1])
	if(p > 0)
	{
		r <- (1 + exp(-psi[2]))^(-1)
		if(p > 1)
		{
			a.1 <- (1 + sum(exp(-psi[3:(p+1)])))^(-1)
			a.j <- a.1
			for(j in 2:p)
			{
				a.j <- c(a.j,exp(-psi[j+1])*a.1)
			}
			a.j <- r*a.j
		} else { a.j <- r }
		a.0 <- c(a.0,a.j)
	}
	return(a.0)
}



