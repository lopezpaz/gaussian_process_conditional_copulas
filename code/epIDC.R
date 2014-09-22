###############################################################################
#
# GAUSSIAN PROCESS CONDITIONAL COPULAS
# fit_copula() on vines.R:31 fits these models to data
#
###############################################################################

# Generalized FITC code
source("gFITC.R")

# RBF kernel width heuristic
sigest <- function (x, frac = 0.5) {
  x <- as.matrix(x)
  m <- dim(x)[1]
  n <- floor(frac*m)
  index <- sample(1:m, n, replace = TRUE)
  index2 <- sample(1:m, n, replace = TRUE)
  temp <- x[index,, drop=FALSE] - x[index2,,drop=FALSE]
  dist <- rowSums(temp^2)
  srange <- 1/quantile(dist[dist!=0],probs=c(0.9,0.5,0.1))
  return(mean(srange[c(1,3)]))
}

# Density of a Gaussian copula in terms of Kendall's tau.
dGAC <- function(x, y, tau) {
	rho <- sin((tau * pi) / 2)
	x <- qnorm(x)
	y <- qnorm(y)
  1/(2*pi*sqrt(1-rho^2))*exp(-1/(2*(1-rho^2))*(x^2+y^2-2
  *rho*x*y))*dnorm(x)^-1*dnorm(y)^-1
}

# Density of a Gaussian copula in terms of Kendall's tau.
logdGAC <- function(x, y, tau) {
	rho <- sin((tau * pi) / 2)
	x <- qnorm(x)
	y <- qnorm(y)
  -log(2*pi)-0.5*log(1-rho^2)-1/(2*(1-rho^2))*(x^2+y^2-2
  *rho*x*y)-dnorm(x,log=T)-dnorm(y,log=T)
}

likelihoodFactor <- function(a, x, y, copulaDensity) {

	tau <- 0.99 * (2 * pnorm(a) - 1)
	copulaDensity(x, y, tau)
}
	
epIDCinternal <- function(X, Z, Zbar, sigma, sigma0, l, 
                          m0, copulaDensity, start = NULL) {
	# We initialize the structure with the problem information
	gFITCinfo <- initGFITCinfo(Z, Zbar, sigma, sigma0, l, m0)

	# We initialize the approximate factors
	f1Hat <- list(eta1 = rep(1e-100, gFITCinfo$n), 
    eta2 = rep(1e-100, gFITCinfo$n), logZ = rep(0, gFITCinfo$n))

	f2Hat <- list(eta1 = rep(1e-100, gFITCinfo$n),
    eta2 = rep(1e-100, gFITCinfo$n))

	# We initialize the posterior approximation
	a <- list(eta1 = rep(0, gFITCinfo$n), eta2 = rep(0, gFITCinfo$n))

	# We refine the second approximate factor
	a$eta2 <- f2Hat$eta2 <- gFITCinfo$diagKnn^-1
	a$eta1 <- f2Hat$eta1 <- rep(0, gFITCinfo$n)

	# We check for an initial solution
	if (!is.null(start)) {
		a     <- start$a
		f1Hat <- start$f1Hat
		f2Hat <- start$f2Hat
	}

	# Main loop of EP
	i           <- 1
	damping     <- 1
	convergence <- FALSE

	while (!convergence && i < 50) {
		aOld <- a

		# We refine the first approximate factor
		f1HatOld <- f1Hat
		repeat {
			f1Hat <- f1HatOld
			for (j in 1 : gFITCinfo$n) {
				if (f2Hat$eta2[ j ] > 0) {
					x      <- X[ j, 1 ]
					y      <- X[ j, 2 ]
					va     <- f2Hat$eta2[ j ]^-1
					ma     <- f2Hat$eta1[ j ] * va
					maxa   <- min(ma + 10 * sqrt(va), 6 * sqrt(gFITCinfo$sigma))
					mina   <- max(ma - 10 * sqrt(va), -6 * sqrt(gFITCinfo$sigma))
					target <- function(a) likelihoodFactor(a, x, y, copulaDensity) * dnorm(a, ma, sqrt(va))

					# We obtain the normalization constant and the first and second moments
					Z <- integrateFast(target, mina, maxa)
					if (Z >= 1e-10) {
						f1Hat$logZ[ j ] <- log(Z)
						maNew           <- integrateFast(function(a) target(a) * a, mina, maxa) / Z
						vaNew           <- integrateFast(function(a) target(a) * (a - maNew)^2, mina, maxa) / Z
						eta1New         <- maNew * vaNew^-1
						eta2New         <- vaNew^-1
						eta1HatNew      <- eta1New - f2Hat$eta1[ j ]
						eta2HatNew      <- eta2New - f2Hat$eta2[ j ]
						f1Hat$eta2[ j ] <- damping * eta2HatNew + (1 - damping) * f1Hat$eta2[ j ]
						f1Hat$eta1[ j ] <- damping * eta1HatNew + (1 - damping) * f1Hat$eta1[ j ]
					}
				}
			}

			if (checkPDposterior(gFITCinfo, f1Hat$eta2))
				break
			else
				damping <- damping * 0.5
		}

		# We refine the second approximate factor.
		ret        <- computeTitledDistribution(gFITCinfo, f1Hat$eta2, f1Hat$eta1)
		eta1HatNew <- ret$mNew / ret$vNew - f1Hat$eta1
		eta2HatNew <- ret$vNew^-1 - f1Hat$eta2

		eta1HatNew[ eta2HatNew < 0 ] <- f2Hat$eta1[ eta2HatNew < 0 ]
		eta2HatNew[ eta2HatNew < 0 ] <- f2Hat$eta2[ eta2HatNew < 0 ]

		f2Hat$eta1 <- damping * eta1HatNew + (1 - damping) * f2Hat$eta1
		f2Hat$eta2 <- damping * eta2HatNew + (1 - damping) * f2Hat$eta2

		# We update the posterior approximation
		a$eta1 <- f1Hat$eta1 + f2Hat$eta1
		a$eta2 <- f1Hat$eta2 + f2Hat$eta2

		# We check for convergence
		change <- max(abs(aOld$eta1 / aOld$eta2 - a$eta1 / a$eta2))
		change <- max(change, abs(aOld$eta2^-1 - a$eta2^-1))

		if (change < 1e-4)
			convergence <- T

		# Annealed damping scheme
		damping <- damping * 0.99
		i       <- i + 1
	}

	# We compute the evidence and its gradient
	logZ         <- computeEvidence(a, f1Hat, f2Hat, gFITCinfo, X, copulaDensity)
	gradientLogZ <- computeDerivativesEvidence(gFITCinfo, f1Hat$eta2, f1Hat$eta1)

	# We are done!
	list(a = a, f1Hat = f1Hat, f2Hat = f2Hat, logZ = logZ,
		gradientLogZ = gradientLogZ, gFITCinfo = gFITCinfo,
		copulaDensity = copulaDensity)
}

computeEvidence <- function(a, f1Hat, f2Hat, gFITCinfo, X, copulaDensity) {
	logZ    <- sum(f1Hat$logZ)
	logZret <- logZ + sum(0.5 * log(2 * pi) -
		0.5 * log(f2Hat$eta2) + 0.5 * log(a$eta2) -
		0.5 * a$eta1^2 / a$eta2 + 0.5 * f1Hat$eta1^2 /
		f1Hat$eta2 + 0.5 * f2Hat$eta1^2 / f2Hat$eta2)

	logZret + getFITCevidence(gFITCinfo, f1Hat$eta2, f1Hat$eta1)
}

# Function that predicts for a new data point
predictIDC <- function(ret, zNew) {
	retInternal <- predictFITC(ret$gFITCinfo, ret$f1Hat$eta2, ret$f1Hat$eta1, zNew)
	mZnew       <- retInternal$m
	vZnew       <- retInternal$v
	t(apply(cbind(mZnew, vZnew), 1, function(x) 0.99 * (2 * pnorm(rnorm(5000,
    x[ 1 ], sqrt(x[ 2 ]))) - 1)))
}

# Function that computes an integral using the Newton-Cotes formula.
integrateFast <- function(f, minx, maxx) {
  x <- seq(minx, maxx, length.out = 101)
  h <- x[ 2 ] - x[ 1 ]

  # We evaluate the function to be integrated
  fvalues1 <- f(x)
  fvalues2 <- t.default(matrix(fvalues1[ 1 : (length(fvalues1) - 1) ], ncol = length(x) / 4, nrow = 4))
  fvalues3 <- cbind(fvalues2, c(fvalues2[ 2 : nrow(fvalues2) , 1 ], fvalues1[ length(fvalues1) ]))
  aux      <- fvalues3 %*% c(7, 32, 12, 32, 7) * 2 * h / 45
  x2       <- seq(x[ 5 ], x[ length(x) ], 4 * h)
  sum(aux)
}

epIDCexternal <- function(X, Z, copulaDensity, m = round(nrow(X) / 10),
                          tau_mle) {
	# We initialize the hyper-parameters
	m0     <- qnorm((tau_mle + 1) / 2)
	sigma  <- 0
	sigma0 <- -10
	l      <- log(sqrt(apply(Z,2,sigest)))
	Zbar   <- matrix(Z[ sample(1 : nrow(Z), m), ], m, ncol(Z))

	# We initialize the gradient optimization process
  ret         <- epIDCinternal(X, Z, Zbar, sigma, sigma0, l, m0, copulaDensity)
	bestSoFar   <- ret
	eps         <- 0.1
	convergence <- F
	iteration   <- 1

	while (!convergence && iteration < 100) {
		sigmaNew  <- sigma + eps * ret$gradientLogZ$dLogZdSigma
		sigma0New <- sigma0 + eps * ret$gradientLogZ$dLogZdSigma0
		lNew      <- l + eps * ret$gradientLogZ$dLogZdl
		ZbarNew   <- Zbar + eps * ret$gradientLogZ$dLogZdXbar
		gFITCinfo <- initGFITCinfo(Z, ZbarNew, sigmaNew, sigma0New, lNew, m0)

		while (!checkPDposterior(gFITCinfo, ret$f1Hat$eta2)) {
			eps       <- eps * 0.5
			sigmaNew  <- sigma + eps * ret$gradientLogZ$dLogZdSigma
			sigma0New <- sigma0 + eps * ret$gradientLogZ$dLogZdSigma0
			lNew      <- l + eps * ret$gradientLogZ$dLogZdl
			ZbarNew   <- Zbar + eps * ret$gradientLogZ$dLogZdXbar
			gFITCinfo <- initGFITCinfo(Z, ZbarNew, sigmaNew, sigma0New, lNew, m0)
		}

		sigma  <- sigmaNew
		sigma0 <- sigma0New
		l      <- lNew
		Zbar   <- ZbarNew
    retNew <- epIDCinternal(X, Z, Zbar, sigma, sigma0, l, m0, copulaDensity,
    ret)

		if (abs(retNew$logZ - ret$logZ) < 1e-2)
			convergence <- T

		if (retNew$logZ < ret$logZ)
			eps <- eps * 0.5
		else {
			if (retNew$logZ > bestSoFar$logZ)
				bestSoFar <- retNew
			eps <- eps * 1.1
		}

		ret       <- retNew
		iteration <- iteration + 1
	}

	bestSoFar
}
