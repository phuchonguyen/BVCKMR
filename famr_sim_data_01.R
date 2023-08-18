## In this scenario:
# 1. The factors that explains most of the variation in the exposures predicts the outcomes.
# 2. Exposures have linear effects and linear interactions with time
# 3. The outcomes are independent, random intercepts are independent across outcomes

library(tidyverse)

q = 2 # Random sloMe and random interceMt
## latent factor dimension
k <- 2
## mixture dimensions
M <- 5 * k # each factor loads on 5 predictors mainly
## covariates dimension
p <- 2
## outcome dimension
Q <- 1
## number of repeated measurement
T <- 10
N = n*T
## id of subjects, stack observations
id <- rep(1:n, each=T) # NOTICE: the order is changed compared to input to the other models
## observation time of each observation
## regular on a grid
time <- rep(1:T, n)  # NOTICE: the order is changed compared to input to the other models
## noise variance for mixtures
sigma_x <- rep(0.5, M) #runif(M, 0.2, 0.5)

# Only sample predictors once
if (i == 1) {
  ## sample effective latent factors
  eta <- matrix(rnorm(n * k), n, k)
  eta_pred <- matrix(rnorm(n * k), n, k)
  ## define a loading matrix for the effective factors
  Theta <- matrix(0, p, k)
  for (kk in 0:1) {
    Theta[kk * (p / 2) + (1:(p / 2)), (kk + 1)] <- rnorm(p/2, 0, 1)
  }
  ## generate mixture n x M
  Z <-
    eta %*% t(Theta) + MASS::mvrnorm(n, mu = rep(0, M), Sigma = diag(sigma_x))
  Z_pred <-
    eta_pred %*% t(Theta) + MASS::mvrnorm(n, mu = rep(0, M), Sigma = diag(sigma_x))
  Z <- apply(Z, 2, scale)
  Z_pred <- apply(Z_pred, 2, scale)
  ## generate covariates n x p
  X <- cbind(rbinom(N, 1, prob = 0.6), rnorm(N))
  X_pred <- cbind(rbinom(N, 1, prob = 0.6), rnorm(N))
  ## generate linear regression coefs
  B <- matrix(runif(q * k, 0, 3.5) * sample(c(-1, 1), q * k, replace = TRUE), q, k)
  B_t <- matrix(runif(q * k, 0, 1.6) * sample(c(-1, 1), q * k, replace = TRUE), q, k)
  ## convert factors' effects to effects in Z
  Ax <-
    solve(t(Theta) %*% diag(1 / sigma_x ^ 2) %*% Theta + diag(1, k, k))
  Ax <- Ax %*% t(Theta) %*% diag(1 / sigma_x ^ 2)
  Bx <- B %*% Ax
  Bx_t <- B_t %*% Ax
  ## effects of covariates
  B_z <- matrix(runif(q * p_covariates, 0.3, 0.5), q, p_covariates)
  ## generate random intercept
  xi_sigma <- 1
  ## generate observation error
  E_sigma <- 0.5
}

## generate outcomes
Y <- eta[id, ] %*% B[i,] +
  (eta[id, ] * tcrossprod(scale(time), rep(1, k))) %*% B_t[i,] +
  X[id, ] %*% B_z[i,] + 
  rnorm(n, 0, sqrt(xi_sigma))[id] +
  rnorm(N, 0, sqrt(E_sigma))

## for evaluation
trueh_pred <- eta_pred[id, ] %*% B[i,] +
  (eta_pred[id, ] * tcrossprod(scale(time), rep(1, k))) %*% B_t[i,]
Y_pred_oracle <- trueh_pred + X_pred[id, ] %*% B_z[i,]
Y_pred <- Y_pred_oracle +
  rnorm(n, 0, sqrt(xi_sigma))[id] +
  rnorm(N, 0, sqrt(E_sigma))

# Make U matrix
age = unique(scale(time))
U				  = matrix(data=rep(0,T*2*n), nrow=T)
U[,1]			= rep(1, T)
U[,2]			= age
counter 		= 1
while(counter < n) {
  U_i					= matrix(data=rep(0,T*2*n), nrow=T)
  U_i[,counter*2+1]	= rep(1, T)
  U_i[,counter*2+2]	= age
  U 					= rbind(U, U_i)
  counter 				= counter+1
}
# Make W matrix
W				= matrix(data=rep(0,T*2*n), nrow=T)
W[,1]			= rep(1, T)
W[,n+1]			= age
counter 		= 2
while(counter <= n) {
  W_i				= matrix(data=rep(0,T*2*n), nrow=T)
  W_i[,counter]	= rep(1, T)
  W_i[,n+counter]	= age
  W 				= rbind(W, W_i)
  counter 		= counter+1
}

WTW = t(W) %*% W
XTX	= t(X) %*% X
UTU = t(U) %*% U