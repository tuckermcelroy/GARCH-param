arch2psi <- function(arch)
{

	#################################
	#	arch2psi by Tucker McElroy
	#
	#	inverts exact parametrization of ARCH(p)
	#	input: coefs a0, a1, ..., ap of ARCH(p)
	#	output: p+1 vector of reals
	###############################################

	p <- length(arch)-1
	psi.0 <- log(arch[1])
	if(p > 0)
	{
		r <- sum(arch[-1])
		alpha.j <- arch[-1]/r
		psi.1 <- -log(1/r - 1)
		if(p > 1)
		{
			psi.j <- psi.1
			for(j in 2:p)
			{
				psi.j <- c(psi.j,-log(alpha.j[j]/alpha.j[1]))
			}
		} else { psi.j <- psi.1 }
		psi.0 <- c(psi.0,psi.j)
	}
	return(psi.0)
}



