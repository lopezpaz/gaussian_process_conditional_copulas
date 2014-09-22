###############################################################################
#
# GAUSSIAN PROCESS CONDITIONAL REGULAR VINE COPULAS 
# see vine_experiment() for a starting point
#
###############################################################################

library(copula)
source("epIDC.R")

###############################################################################
#
# FIT_COPULA 
#
# This is a wrapper that fits a (conditional) copula to two vectors of rank
# statistics (and a matrix of conditioning rank statistics). Three types of
# copulas are supported:
#
#   - method = "gc": Gaussian Copula, ignores conditioning covariates.
#   - method = "ep": Gaussian Processes with Expectation Propagation.
#
###############################################################################

safe <- function(x) { x[x > (1-1e-4)] <- (1-1e-4); x[x < 1e-4] <- 1e-4; x }

hGAC <- function(j, k, ktau) {
  rho <- sin(ktau*pi/2)
  as.matrix(pnorm((qnorm(j)-rho*qnorm(k))/sqrt(1-rho**2)))
}

fit_copula <- function(x_j, x_k, ktau, x_d = NULL,
              method = "gc", always=NULL, pseudo = 20) {
  x_d     <- cbind(x_d, always)
  tau_mle <- 2*asin(fitCopula(normalCopula(0.5),
  cbind(x_j,x_k))@copula@parameters)/pi

  if((method == "gc") || (length(x_d) == 0)) {
    return(list(pdf = function(j,k,D) dGAC(safe(j),safe(k),tau_mle),
                h   = function(j,k,D) hGAC(safe(j),safe(k),tau_mle),
                tau = tau_mle))
  } else if(method == "ep") {
    obj  <- epIDCexternal(cbind(x_j,x_k),x_d,dGAC,pseudo,tau_mle)

    pdf <- function(j,k,d)
      apply(dGAC(matrix(j,length(j),5000),
            matrix(k,length(k),5000), predictIDC(obj,d)),1,mean)
    h   <- function(j,k,d)
      apply(hGAC(matrix(j,length(j),5000),
            matrix(k,length(k),5000), predictIDC(obj,d)),1,mean)
    taus <- function(d) predictIDC(obj,d)

    return(list(pdf = pdf,
                h   = h,
                taus=taus))
  } 
}

###############################################################################
#
# HELPER FUNCTIONS 
#
# Primm's algorithm, Set Symmetric Difference, Safe numbers for copulas
#
###############################################################################

symdiff <- function(s1, s2) union(setdiff(s1,s2), setdiff(s2,s1))

prim <- function(G) {
  V     <- 1:nrow(G) 
  V_new <- 1 
  E_new <- NULL

  best_idx <- 1

  while(setequal(V, V_new) == FALSE) {
    best <- -Inf

    for(i in intersect(V, V_new)) {
      for(j in setdiff(V, V_new)) {
        if((G[i,j] > best) && (G[i,j] != 2)) {
          best <- G[i,j]
          best_idx <- c(i,j)
        }
      }
    }

    V_new <- c(V_new, best_idx[2])
    E_new <- c(E_new, list(best_idx))
  }

  E_new
}

###############################################################################
#
# VINE_FIT 
#
# This function fits a Regular Vine Distribution formed by Gaussian Copulas 
# of one of three kinds:
#   - method = "gc": Conditioned Copulas are asumed constant respect
#                            to conditioned covariates
#   - method = "ep": Conditioned Copulas parameters are functions estimated in
#                    terms of Gaussian Processes and Expectation Propagation
#
# The fitting technique is sequential: We start off by inferring a Maximum
# Spanning Tree given a complete graph in which each node represents a
# a variable. Each edge of this first tree will correspond to a copula. Each
# edge of subsequent trees will correspond to a conditioned copula.
#
###############################################################################

vine_fit <- function(d, ntrees = ncol(d) - 1, method = "ep",
                     always = NULL, pseudo = 20)
{
  # Initialize vector of copulas and parameters as empty
  cop_f <- cop_j <- cop_k <- cop_c <- cop_d <- cop_r <- NULL
  
  vine  <- list(cop_f = cop_f, cop_j = cop_j,
                cop_k = cop_k, cop_d = cop_d)

  # Get Kendall Tau's for all pairs of variables (TODO: Use fast Kendall)
  W <- cor(d, method="kendall")

  # Get Maximum Spanning Tree given the absolute value of
  # Kendall's Tau
  T <- prim(abs(W))

  # Fill Copulas and Parameters for first tree 
  for(tedge in T) {
    j <- tedge[1]
    k <- tedge[2]
    # cop_j: First argument of the copula
    vine$cop_j <- c(vine$cop_j, list(j))
    # cop_k: Second argument of the copula
    vine$cop_k <- c(vine$cop_k, list(k))
    # cop_c: Constraing set of the copula
    vine$cop_c <- c(vine$cop_c, list(c(j,k)))
    # cop_d: Conditioned set of the copula
    vine$cop_d <- c(vine$cop_d, list(NULL))
    # cop_r: Father copulas
    vine$cop_r <- c(vine$cop_r, list(NULL))
    # cop_f: Functions of the copula (pdf and h)
    vine$cop_f <- c(vine$cop_f, list(fit_copula(d[,j], d[,k], 
                    W[j,k], method=method, always=always, pseudo=pseudo)))

  }

  # This loop will build the rest of the trees, that contain
  # Conditional Copulas
  if(ntrees > 1) {
    for(ntree in 2:ntrees) {
      new_jk <- new_c <- new_d <- new_U <- new_V <- new_R <- NULL
      # Kendall Matrix and Edge Relationship matrix 
      W <- matrix(2, ncol(d)-ntree+1, ncol(d)-ntree+1)
      R <- matrix(2, ncol(d)-ntree+1, ncol(d)-ntree+1)
      # Given previous tree, we need to find pairs of edges
      # with a node in common:
      seen_edges <- 1
      for(e1 in 1:(length(T)-1)) {
        for(e2 in (e1+1):length(T)) {
          if(length(intersect(T[[e1]],T[[e2]])) == 1) {
            # We found two edges with a node in common. This pair now produces
            # a single new edge in the new graph.
            ce1 <- e1+(0.5*(ntree-2)*(2*ncol(d)-ntree+1))
            ce2 <- e2+(0.5*(ntree-2)*(2*ncol(d)-ntree+1))
            # Compute new conditioned, conditioning and constraint set
            # For the new formed copula.
            this_jk <- symdiff  (vine$cop_c[[ce1]], vine$cop_c[[ce2]])
            this_c  <- union    (vine$cop_c[[ce1]], vine$cop_c[[ce2]])
            this_d  <- intersect(vine$cop_c[[ce1]], vine$cop_c[[ce2]])
            # Get conditioned vector U
            U <- condition_data(vine, d, this_jk[1], this_d, always)
            # Get conditioned vector V
            V <- condition_data(vine, d, this_jk[2], this_d, always)
            # Get associated Kendall's tau correlation
            W[e1,e2] <- W[e2,e1] <- cor(U,V,method="kendall")
            R[e1,e2] <- R[e2,e1] <- seen_edges
            seen_edges <- seen_edges + 1;
            # Save data and parameters to avoid recalculation
            new_jk <- c(new_jk, list(this_jk))
            new_c  <- c(new_c , list(this_c))
            new_d  <- c(new_d , list(this_d))
            new_U  <- c(new_U, list(U))
            new_V  <- c(new_V, list(V))
            new_R  <- c(new_R, list(c(ce1,ce2)))
          }
        }
      }
      
      # Get next maximum spanning tree
      T <- prim(abs(W))

      for(e in 1:length(T)) {
        edge_idx <- R[T[[e]][1],T[[e]][2]]
        # Insert computed parameters in the vine model
        vine$cop_j <- c(vine$cop_j, list(new_jk[[edge_idx]][1]))
        vine$cop_k <- c(vine$cop_k, list(new_jk[[edge_idx]][2]))
        vine$cop_c <- c(vine$cop_c, list(new_c [[edge_idx]]))
        vine$cop_d <- c(vine$cop_d, list(new_d [[edge_idx]]))
        vine$cop_r <- c(vine$cop_r, list(new_R [[edge_idx]]))
        vine$cop_f <- c(vine$cop_f, list(fit_copula(
            as.matrix(new_U[[edge_idx]]),
            as.matrix(new_V[[edge_idx]]),
            W[T[[e]][1],T[[e]][2]],
            as.matrix(d[,new_d[[edge_idx]]]),
            method=method,always=always, pseudo=pseudo)))
       }
    }
  }
  vine
}

###############################################################################
#
# CONDITION_DATA 
#
# This function applies the recursive rule that allow us to compute the
# conditional CDF's that are passed as arguments to the conditioned copulas.
#
# TODO: Right now back-tracking is not implemented, i.e, conditioning on 10
# levels will require computing the conditioning UP TO level 9, instead of
# storing the conditioning data and performing only the last level.
#
###############################################################################

condition_data <- function(vine, d, xi, di, always = NULL) {
  if(length(di) == 0) {
    return(d[,xi])
  } else {
    for(i in 1:length(vine$cop_f)) {
      yi <- setdiff(di, vine$cop_d[[i]])
      if(length(yi) == 1) {
         if(((vine$cop_j[[i]] == xi) && (vine$cop_k[[i]] == yi)) ||
            ((vine$cop_j[[i]] == yi) && (vine$cop_k[[i]] == xi))) {
            di <- di[-which(di==yi)]
            return(vine$cop_f[[i]]$h(condition_data(vine, d, xi, di, always),
                                condition_data(vine, d, yi, di, always),
                                cbind(as.matrix(d[,vine$cop_d[[i]]]), always)))
         }
      }
    }
  }
}

###############################################################################
#
# EVALUATE_VINE 
#
# This function loops over the array of copulas contained in the Regular
# Vine Distribution and accumulates the total log-likelihood of the vine at
# x.
#
###############################################################################

vine_eval <- function(vine, x, always = NULL, eval.levels = NULL)
{
  res <- 0

  for(i in 1:length(vine$cop_f))
  {
    if(is.null(eval.levels) == FALSE) {
      if(length(vine$cop_d[[i]]) == eval.levels) {
        return(res);
      }
    }

    # Apply recursive rule to get j|D
    j <- condition_data(vine, x, vine$cop_j[[i]], vine$cop_d[[i]], always)
    # Apply recursive rule to get k|D
    k <- condition_data(vine, x, vine$cop_k[[i]], vine$cop_d[[i]], always)
    # Conditioning covariates
    d <- x[,vine$cop_d[[i]]]

    # Evaluate bivariate copula
    res <- res + log(vine$cop_f[[i]]$pdf(as.matrix(j), 
                                         as.matrix(k),
                                         cbind(as.matrix(d),always)))
  }
  res
}

###############################################################################
#
# VINE_EXPERIMENT
#
# 1. Reads data from file fname
# 2. Maps the data to the unit hypercube
# 3. Splits data in train/test partitions, with percentage of train 30%
# 4. Fits Expectation-Propagation Vine, Local Likelihood Vine and Gaussian Vine
# 5. Evaluates the log-likelihood of the three models
# 6. Saves result in file "vine.res"
#
###############################################################################

vine_experiment <- function(fname, npoints = 300, levels = 4)
{
  # Read data
  d <- as.matrix(read.table(fname))

  # Apply feature-wise Probability Integral Transforms
  for(i in 1:ncol(d)) d[,i] <- safe(ecdf(d[,i])(d[,i]))

  d         <- d[sample(1:nrow(d), npoints),replace=TRUE]
  train_idx <- sample (1:nrow(d), nrow(d)*0.5)
  test_idx  <- setdiff(1:nrow(d), train_idx)
  d_train   <- as.matrix(d[train_idx,])
  d_test    <- as.matrix(d[test_idx ,])
  vine_ep   <- vine_fit(d_train, method="ep",ntrees=levels)
  vine_gc   <- vine_fit(d_train, method="gc",ntrees=levels)
  er_ep     <- mean(vine_eval(vine_ep,d_test,eval.levels=levels))
  er_gc     <- mean(vine_eval(vine_gc,d_test,eval.levels=levels))
  print(paste(fname, er_gc, er_ep))
}

vine_experiment("weather.txt")

