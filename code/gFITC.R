#
# Script with the routines necessary to implement the generalized FITC approximation.
#
# Naish-Guzman, A. & Holden, S. B.
# The Generalized FITC Approximation
# Advances in Neural Information Processing Systems 20, 2007
#

##
# Function which determines if the posterior covariance
# matrix is positive definite
#

checkPDposterior <- function(gFITCinfo, tauTilde) {

	# Function which checks if the posterior covariance matrix
	# is positive definite

	M <- rot180(diag(gFITCinfo$m) + t(gFITCinfo$PRt) %*%
		(matrix(tauTilde / (1 + gFITCinfo$D *
		tauTilde), gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt))

	eigenValues <- eigen(M)$values
	eigenValues[ abs(eigenValues) < 1e-6 ] <- 0

	if (any(eigenValues <= 0))
		F
	else
		T
}

##
# Function that initializes the struture with the problem information.
#
# @param	X	n x d matrix with the data points.
# @param	Xbar	m x d matrix with the pseudo inputs.
# @param	sigma	scalar with the log-amplitude of the GP.
# @param	sigma0	scalar with the log-noise level in the GP.
# @param	l	d-dimensional vector with the log-lengthscales.
# @param	m0	scalar with the mean of the approximate prior
# 
# @return	gFITCinfo	List with the problem information 
#

initGFITCinfo <- function(X, Xbar, sigma, sigma0, l, m0 = 0) {

	# We initialize the structure with the data and the kernel
	# hyper-parameters

	gFITCinfo <- list()
	gFITCinfo$X <- X
	gFITCinfo$Xbar <- Xbar
	gFITCinfo$m <- nrow(Xbar)
	gFITCinfo$d <- ncol(Xbar)
	gFITCinfo$n <- nrow(X)
	gFITCinfo$sigma <- exp(sigma) + 1e-3
	gFITCinfo$sigma0 <- exp(sigma0)
	gFITCinfo$l <- sqrt(exp(l))

	# We compute the kernel matrices

	gFITCinfo$XbarScaled <- Xbar * matrix(gFITCinfo$l, nrow(Xbar),
		ncol(Xbar), byrow = T)
	gFITCinfo$XScaled <- X * matrix(gFITCinfo$l, nrow(X), ncol(X),
		byrow = T)
	gFITCinfo$Kmm <- computeKmm(gFITCinfo$XbarScaled, gFITCinfo$sigma,
		gFITCinfo$sigma0)
	gFITCinfo$Knm <- computeKnm(gFITCinfo$XScaled, gFITCinfo$XbarScaled,
		gFITCinfo$sigma)
	gFITCinfo$P <- gFITCinfo$Knm
	gFITCinfo$R <- cholInverse(gFITCinfo$Kmm)
	gFITCinfo$PRt <- gFITCinfo$P %*% t(gFITCinfo$R)

	# We compute the diagonal matrices

	gFITCinfo$diagKnn <- computeDiagKnn(X, gFITCinfo$sigma,
		gFITCinfo$sigma0)
	gFITCinfo$D <- gFITCinfo$diagKnn - (gFITCinfo$PRt)^2 %*%
		rep(1, gFITCinfo$m)

	# We compute A^-1 m0 efficiently using the woodbury formula and cholesky representations

	L0 <- cholInverse(diag(gFITCinfo$m) + t(gFITCinfo$PRt) %*%
		(matrix(gFITCinfo$D^-1, gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt))

	gFITCinfo$AInvM0 <- as.vector(gFITCinfo$D^-1 * m0 - 
		(matrix(gFITCinfo$D^-1, gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt) %*% t(L0) %*% (L0 %*% 
		(t(gFITCinfo$PRt) %*%  (gFITCinfo$D^-1 * m0))))

	gFITCinfo$m0 <- m0

	gFITCinfo$RtRPAInvM0 <- t(gFITCinfo$R) %*% (t(gFITCinfo$PRt) %*% gFITCinfo$AInvM0)

	# We are done

	gFITCinfo
}

##
# Function that makes predictions using the generalized FITC approximation.
#
# @param	gFITCinfo	List with the problem information
#				(see initializegFITCinfo).
# @param	tauTilde 	n-dimensional vector with the inverse
#				variances of the Gaussian approximation
#				to the likelihood factor.
# @param	muTilde 	n-dimensional vector with the inverse
#				variances times the mean of the Gaussian
#				approximation to the likelihood factor.
# @param	Xtest		r x d matrix with the test data points.
#				
# @return	ret		List with the mean and variances of the
#				predictions at the test points.
#

predictFITC <- function(gFITCinfo, tauTilde, muTilde, Xtest) {

	# We compute the FITC prediction

	Xtest <- Xtest * matrix(gFITCinfo$l, nrow(Xtest),
		ncol(Xtest), byrow = T)
	pStar <- computeKnm(Xtest, gFITCinfo$XbarScaled, gFITCinfo$sigma)

	dStar <- computeDiagKnn(Xtest, gFITCinfo$sigma, gFITCinfo$sigma0) -
		(pStar %*%t(gFITCinfo$R))^2 %*% rep(1, gFITCinfo$m)

	Dnew <- gFITCinfo$D / (1 + gFITCinfo$D * tauTilde)
	Pnew <- matrix(1 / (1 + gFITCinfo$D * tauTilde), gFITCinfo$n,
		gFITCinfo$m) * gFITCinfo$P

	Rnew <- backsolve(rot180(t(chol(rot180(diag(gFITCinfo$m) +
		t(gFITCinfo$PRt) %*% (matrix(tauTilde / (1 + gFITCinfo$D *
		tauTilde), gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt))))),
		gFITCinfo$R)

	gammaNew <- t(Rnew) %*% (Rnew %*% (t(Pnew) %*% (muTilde + gFITCinfo$AInvM0)))

	# We obtain the new marginals

	mPrediction <- pStar %*% gammaNew 

	# We add the contribution of the prior mean

	mPrediction <- mPrediction + gFITCinfo$m0 - pStar %*% gFITCinfo$RtRPAInvM0

	vPrediction <- dStar + (pStar %*% t(Rnew))^2 %*% rep(1, gFITCinfo$m)

	list(m = mPrediction, v = vPrediction)
}


##
# Function that evaluates the normalization constant of the product
# of the FITC prior and a multivariate Gaussian density with diagonal
# correlation matrix.
#
# @param	gFITCinfo	List with the problem information
#				(see initializegFITCinfo).
# @param	tauTilde 	n-dimensional vector with the inverse
#				variances of the Gaussian approximation
#				to the likelihood factor.
# @param	muTilde 	n-dimensional vector with the inverse
#				variances times the mean of the Gaussian
#				approximation to the likelihood factor.
#
# @return	logZ		Logarithm of the normalization constant.
#

getFITCevidence <- function(gFITCinfo, tauTilde, muTilde) {

	Dnew <- gFITCinfo$D / (1 + gFITCinfo$D * tauTilde)
	Pnew <- matrix(1 / (1 + gFITCinfo$D * tauTilde), gFITCinfo$n,
		gFITCinfo$m) * gFITCinfo$P

	Rnew <- backsolve(rot180(t(chol(rot180(diag(gFITCinfo$m) +
		t(gFITCinfo$PRt) %*% (matrix(tauTilde / (1 + gFITCinfo$D *
		tauTilde), gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt))))),
		gFITCinfo$R)

	upsilon <- gFITCinfo$AInvM0 + muTilde
	aNew <- Dnew * upsilon
	gammaNew <- t(Rnew) %*% (Rnew %*% (t(Pnew) %*% upsilon))

	# We obtain the new marginal means

	mNew <- as.double(aNew) + as.double(Pnew %*% gammaNew)

	# (Rnew appears as a consequence of the matrix determinant lema)

	logZ <- -gFITCinfo$n / 2 * log(2 * pi) + sum(log(diag(Rnew))) -
		sum(log(diag(gFITCinfo$R))) - 0.5 * sum(log(1 + tauTilde *
		gFITCinfo$D)) + 0.5 * sum((muTilde + gFITCinfo$AInvM0) *
		mNew) - 0.5 * sum(muTilde^2 / tauTilde) - 0.5 * sum(gFITCinfo$m0 * gFITCinfo$AInvM0)

	logZ
}

##
# Function that computes the marginal means and variances of the product
# of the FITC prior and a multivariate Gaussian density with diagonal
# correlation matrix.
#
# @param	gFITCinfo	List with the problem information
#				(see initializegFITCinfo).
# @param	tauTilde 	n-dimensional vector with the inverse
#				variances of the Gaussian approximation
#				to the likelihood factor.
# @param	muTilde 	n-dimensional vector with the inverse
#				variances times the mean of the Gaussian
#				approximation to the likelihood factor.
#
# @return	ret		A list with the marginal means and variances.
#

computeTitledDistribution <- function(gFITCinfo, tauTilde, muTilde) {

	Dnew <- gFITCinfo$D / (1 + gFITCinfo$D * tauTilde)
	Pnew <- matrix(1 / (1 + gFITCinfo$D * tauTilde), gFITCinfo$n,
		gFITCinfo$m) * gFITCinfo$P

	Rnew <- backsolve(rot180(t(chol(rot180(diag(gFITCinfo$m) +
		t(gFITCinfo$PRt) %*% (matrix(tauTilde / (1 + gFITCinfo$D *
		tauTilde), gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt))))),
		gFITCinfo$R)

	upsilon <- gFITCinfo$AInvM0 + muTilde
	aNew <- Dnew * upsilon  
	gammaNew <- t(Rnew) %*% (Rnew %*% (t(Pnew) %*% upsilon))

	# We obtain the new marginal means

	vNew <- as.double(Dnew + (Pnew %*% t(Rnew))^2 %*% rep(1, gFITCinfo$m))
	mNew <- as.double(aNew) + as.double(Pnew %*% gammaNew)

	list(mNew = mNew, vNew = vNew)
}

##
# Function which rotates a matrix 180 degreees.
#

rot180 <- function(M) {

	matrix(rev(as.double(M)), nrow(M), ncol(M))
}

##
# Function which computes the cholesky decomposition of the inverse
# of a particular matrix.
#
# @param	M	m x m positive definite matrix.
#
# @return	L	m x m upper triangular matrix such that
#			M^-1 = L %*% t(L)
#

cholInverse <- function(M) {

	rot180(forwardsolve(t(chol(rot180(M))), diag(nrow(M))))
}

##
# Function which computes the kernel matrix for the pseudo-inputs.
#
# @param	Xbar	m x d matrix with the m pseudo inputs.
# @param	sigma	scalar with the amplitude of the GP.
# @param	sigma0	scalar with the noise level in the GP.
#
# @return	Kmm	m x m covariance matrix for the pseudo-inputs.
#

computeKmm <- function(Xbar, sigma, sigma0) {

	m <- nrow(Xbar)
	Q <- matrix(apply(Xbar^2, 1, sum), m, m)
	distance <- Q + t(Q) - 2 * Xbar %*% t(Xbar)
	K <- exp(-0.5 * distance)
	sigma * K + diag(sigma0, m) + diag(1e-4, m)
}

##
# Function which computes the kernel matrix between the pseudo-inputs and
# the data points.
#
# @param	X	n x d matrix with the n data points.
# @param	Xbar	m x d matrix with the m pseudo inputs.
# @param	sigma	scalar with the amplitude of the GP.
#
# @return	Knm	n x m covariance matrix between the data points and
#			the pseudo-inputs.
#

computeKnm <- function(X, Xbar, sigma) {

	n <- nrow(X)
	m <- nrow(Xbar)
	Q <- matrix(apply(X^2, 1, sum), n, m)
	Qbar <- matrix(apply(Xbar^2, 1, sum), n, m, byrow = T)
	distance <- Qbar + Q - 2 * X %*% t(Xbar)
	sigma * exp(-0.5 * distance)
}

##
# Function which computes the diagonal of the kernel matrix for the data
# points.
#
# @param	X 		n x d matrix with the n data points.
# @param	sigma		scalar with the amplitude of the GP.
# @param	sigma0		scalar with the noise level in the GP.
#
# @return	diagKnn		n-dimensional vector with the diagonal of the
#				kernel matrix for the data points.
#

computeDiagKnn <- function(X, sigma, sigma0) {

	n <- nrow(X)
	rep(sigma, n) + 1e-4 + sigma0
}

##
# Function which computes the derivative of the evidence
# with respect to the kernel hyper-parameters and the pseudo inputs.
#
# @param	gFITCinfo	List with the problem information
#				(see initializegFITCinfo).
# @param	tauTilde 	n-dimensional vector with the inverse
#				variances of the Gaussian approximation
#				to the likelihood factor.
# @param	muTilde 	n-dimensional vector with the inverse
#				variances times the mean of the Gaussian
#				approximation to the likelihood factor.
#
# @return	ret		A list with the derivatives of the evidence.
#

computeDerivativesEvidence <- function(gFITCinfo, tauTilde, muTilde) {

	# We compute the auxiliary matrices and vectors

	mHat <- muTilde / tauTilde - gFITCinfo$m0
	A <- (tauTilde^-1 + gFITCinfo$D)^-1 
	B <- cholInverse(diag(gFITCinfo$m) + t(gFITCinfo$PRt) %*%
		(matrix(A, gFITCinfo$n, gFITCinfo$m) * gFITCinfo$PRt))
	C <- matrix(A, gFITCinfo$n, gFITCinfo$m) * (gFITCinfo$PRt %*% t(B))
	e <- (tauTilde^-1 + gFITCinfo$D)^-1 * mHat - C %*% (t(C) %*% mHat)
	PRtR <- gFITCinfo$PRt %*% gFITCinfo$R
	M1 <- (t(PRtR) %*% C) %*% t(C) -
		t(PRtR) * matrix(A, gFITCinfo$m, gFITCinfo$n, byrow = T)
	M2 <- 0.5 * (((t(PRtR) %*% C) %*% t(C) - t(PRtR) *
		matrix(A, gFITCinfo$m, gFITCinfo$n, byrow = T)) %*% PRtR)
	v1 <- 1 / sqrt(2) * (C^2 %*% rep(1, gFITCinfo$m) - A)
	v2 <- rep(1 / sqrt(2), gFITCinfo$n)
	v3 <- t(gFITCinfo$PRt %*% gFITCinfo$R) %*% e
	e2 <- 1 / sqrt(2) * e
	M1prima <- M1 - 2 * t(PRtR) * matrix(v1 * v2 + e2^2,
		gFITCinfo$m, gFITCinfo$n, byrow = T)
	M2prima <- (t(PRtR) * matrix(v1 * v2 + e2^2, gFITCinfo$m,
		gFITCinfo$n, byrow = T)) %*% PRtR - M2

	# We evaluate the derivative of P, Kmm and diag(Knn) with respect to
	# sigma

	dPdSigma <- gFITCinfo$P / gFITCinfo$sigma * (gFITCinfo$sigma - 1e-3)
	dKmmdSigma <- (gFITCinfo$Kmm - diag(gFITCinfo$sigma0, gFITCinfo$m) -
		diag(1e-4, gFITCinfo$m)) / gFITCinfo$sigma *
		(gFITCinfo$sigma - 1e-3)
	dDiagKnndSigma <- (gFITCinfo$diagKnn - gFITCinfo$sigma0 - 1e-4) /
		gFITCinfo$sigma * (gFITCinfo$sigma - 1e-3)

	# We compute the derivative of the log evience with respect to sigma

	term1 <- sum((v1 * v2 + e2^2 ) * dDiagKnndSigma)
	term2 <- sum(M1prima * t(dPdSigma))
	term3 <- sum(M2prima * dKmmdSigma)
	term4 <- as.double(t(e) %*% dPdSigma %*% v3)
	term5 <- as.double(-0.5 * t(v3) %*% dKmmdSigma %*% v3)
	dLogZdSigma <- term1 + term2 + term3 + term4 + term5

	# We evaluate the derivative of P, Kmm and diag(Knn) with respect to
	# sigma0

	dPdSigma0 <- matrix(0, gFITCinfo$n, gFITCinfo$m)
	dKmmdSigma0 <- diag(gFITCinfo$sigma0, gFITCinfo$m)
	dDiagKnndSigma0 <- rep(gFITCinfo$sigma0, gFITCinfo$n)

	# We compute the derivative of the log evience with respect to sigma0

	term1 <- sum((v1 * v2 + e2^2 ) * dDiagKnndSigma0)
	term2 <- sum(M1prima * t(dPdSigma0))
	term3 <- sum(M2prima * dKmmdSigma0)
	term4 <- as.double(t(e) %*% dPdSigma0 %*% v3)
	term5 <- as.double(-0.5 * t(v3) %*% dKmmdSigma0 %*% v3)
	dLogZdSigma0 <- term1 + term2 + term3 + term4 + term5

	# We compute the derivatives with respect to the lengthscales

	dLogZdl <- rep(0, length(gFITCinfo$l))
	for (i in 1 : length(gFITCinfo$l)) {

		# We evaluate the derivative of P, Kmm and diag(Knn) with
		# respect to sigma0

		Q <- matrix(gFITCinfo$XScaled[ , i ]^2, gFITCinfo$n,
			gFITCinfo$m)
		Qbar <- matrix(gFITCinfo$XbarScaled[ , i ]^2, gFITCinfo$n,
			gFITCinfo$m, byrow = T)
		distance <- Qbar + Q - 2 * gFITCinfo$XScaled[ , i ] %*%
			t(gFITCinfo$XbarScaled[ , i ])
		dPdl <- -gFITCinfo$P * distance

		Q <- matrix(gFITCinfo$XbarScaled[ , i ]^2, gFITCinfo$m,
			gFITCinfo$m)
		Qbar <- matrix(gFITCinfo$XbarScaled[ , i ]^2, gFITCinfo$m,
			gFITCinfo$m, byrow = T)
		distance <- Qbar + Q - 2 * gFITCinfo$XbarScaled[ , i ] %*%
			t(gFITCinfo$XbarScaled[ , i ])
		dKmmdl <- -(gFITCinfo$Kmm -
			diag(gFITCinfo$sigma0, gFITCinfo$m) -
			diag(1e-4, gFITCinfo$m)) * distance

		dDiagKnndl <- rep(0, gFITCinfo$n)

		# We compute the derivative of the log evience with respect to
		# sigma0

		term1 <- sum((v1 * v2 + e2^2 ) * dDiagKnndl)
		term2 <- sum(M1prima * t(dPdl))
		term3 <- sum(M2prima * dKmmdl)
		term4 <- as.double(t(e) %*% dPdl %*% v3)
		term5 <- as.double(-0.5 * t(v3) %*% dKmmdl %*% v3)
		dLogZdl[ i ] <- (term1 + term2 + term3 + term4 + term5) * 0.5 *
			gFITCinfo$l[ i ]^(-0.5)
	}

	# We compute the derivatives with respect to the pseudo inputs

	dLogZdXbar <- matrix(0, gFITCinfo$m, gFITCinfo$d)
	for (i in 1 : length(gFITCinfo$l)) {

		# We evaluate the derivative of P, Kmm and diag(Knn) with
		# respect to the pseudo-inputs

		distance <- matrix(gFITCinfo$XScaled[ , i ], gFITCinfo$n,
			gFITCinfo$m) - matrix(gFITCinfo$XbarScaled[ , i ],
			gFITCinfo$n, gFITCinfo$m, byrow = T)
		dPdXbar <- gFITCinfo$P * distance * gFITCinfo$l[ i ]

		distance <- (matrix(gFITCinfo$XbarScaled[ , i ], gFITCinfo$m,
			gFITCinfo$m) - matrix(gFITCinfo$XbarScaled[ , i ],
			gFITCinfo$m, gFITCinfo$m, byrow = T))
		dKmmdXbar <- (gFITCinfo$Kmm - diag(gFITCinfo$sigma0,
			gFITCinfo$m) - diag(1e-4, gFITCinfo$m)) * distance *
			gFITCinfo$l[ i ]

		dDiagKnndXbar <- rep(0, gFITCinfo$n)

		# We compute the derivative of the log evience with respect to
		# the pseudo-inputs.

		term1 <- rep(sum((v1 * v2 + e2^2 ) * dDiagKnndXbar), gFITCinfo$m)
		term2 <- (M1prima * t(dPdXbar)) %*% rep(1, gFITCinfo$n)
		term3 <- -2 * (M2prima * dKmmdXbar) %*% rep(1, gFITCinfo$m)
		term4 <- as.double(t(dPdXbar) %*% e * v3)
		term5 <- -t(dKmmdXbar) %*% v3 * v3
		dLogZdXbar[ , i ] <- term1 + term2 + term3 + term4 + term5
	}

	# We compute the gradient with respect to the mean of the prior Gaussian Process

	dLogZdm0 <- sum(e)

	# We return a list with the gradient information

	list(dLogZdSigma = dLogZdSigma, dLogZdSigma0 = dLogZdSigma0,
		dLogZdl = dLogZdl, dLogZdXbar = dLogZdXbar, dLogZdm0 = dLogZdm0)
}


