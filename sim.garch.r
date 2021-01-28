sim.garch <- function(coef,T,init,df)
{
	#################################
	#	sim.garch by Tucker McElroy
	#
	#	simulates GARCH(p,q) of length T, p > 0 
	#		(but q=0 is allowed)
	#	input: coefs a0, a1, ..., ap, b1,..., bq 
	#		of GARCH(p,q), T is desired length,
	#		init is p vector of initial values
	#	output: simulation with Gaussian or t errors (df = degrees freedom)
	###############################################

	p <- length(init)
	q <- length(coef)-1-p
	acoef <- coef[1:(p+1)]
	bcoef <- NULL
	if(q > 0) { bcoef <- coef[(p+2):(p+q+1)] }
	theta <- bcoef
	phi <- acoef[-1]
	psi <- c(1,ARMAtoMA(ar=theta,ma=NULL,lag.max=T))
	psi <- polymul(psi,phi) 
	data <- init
	for(t in (p+1):T)
	{
		if(df==Inf) { z <- rnorm(1) } else { z <- rt(1,df=df) }
		new.sigt <- acoef[1]/(1-sum(theta)) + sum(psi[1:(t-1)]*data[(t-1):1]^2)
		data <- c(data,z*sqrt(new.sigt))
	} 

	return(data)
}
