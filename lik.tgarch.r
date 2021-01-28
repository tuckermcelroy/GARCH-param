lik.tgarch <- function(psi,data,p.order)
{
	#################################
	#	lik.tgarch by Tucker McElroy
	#
	#	evaluates t-GARCH(p,q) lik using euclidean parametrization, 
	#		p.order > 0
	#	input: psi is a p+q+2 vector of reals, data is T-vector of data,
	#		last parameter corresponds to degrees of freedom
	#	output: value of -2 times log t lik
	###############################################

	p <- p.order
	q <- length(psi)-p-2
	T <- length(data)
	df <- 2 + exp(psi[p+q+2])
	coef <- psi2arch(psi[1:(p+q+1)])
	acoef <- coef[1:(p+1)]
	bcoef <- NULL
	if(q > 0) { bcoef <- coef[(p+2):(p+q+1)] }
	theta <- bcoef
	phi <- acoef[-1]
	psi <- c(1,ARMAtoMA(ar=theta,ma=NULL,lag.max=T))
	psi <- polymul(psi,phi) 

	lik <- 0
	for(t in (p+1):T)
	{
		new.sigt <- sqrt(acoef[1]/(1-sum(theta)) + sum(psi[1:(t-1)]*data[(t-1):1]^2))
		lik <- lik + (-2)*log(dt(data[t]/new.sigt,df=df)/new.sigt) 
	}
	return(lik)
}
