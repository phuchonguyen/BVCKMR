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
## sample effective latent factors
eta <- matrix(rnorm(n * k), n, k)
eta_pred <- matrix(rnorm(n * k), n, k)
## define a loading matrix for the effective factors
Theta <- matrix(rnorm(M * k, 0, 0.1), M, k)
for (kk in 0:1) {
  Theta[kk * (M / 2) + (1:(M / 2)), (kk + 1)] <- runif(M / 2, 1.2, 2.2)
}
# for (kk in 0:1) {
#   Theta[kk*(M/2) + (1:5), (kk+3)] <- runif(5, 0.2, 0.5)
# }
## generate mixture n x M
Z <-
  eta %*% t(Theta) + MASS::mvrnorm(n, mu = rep(0, M), Sigma = diag(sigma_x))
Z_pred <-
  eta_pred %*% t(Theta) + MASS::mvrnorm(n, mu = rep(0, M), Sigma = diag(sigma_x))
## generate covariates n x p
X <- cbind(rbinom(n * T, 1, prob = 0.6), rnorm(n * T))
X_pred <- cbind(rbinom(n * T, 1, prob = 0.6), rnorm(n * T))
## generate linear regression coefs
B <- matrix(0, Q, k)
B[, 1:2] <- runif(Q * 2, 2, 4) * sample(c(-1, 1), Q * 2, replace = TRUE)
B_t <- matrix(0, Q, k)
B_t[, 1:2] <- runif(Q, 0, 2) * sample(c(-1, 1), Q * 2, replace = TRUE)
## convert factors' effects to effects in Z
Ax <-
  solve(t(Theta) %*% diag(1 / sigma_x ^ 2) %*% Theta + diag(1, k, k))
Ax <- Ax %*% t(Theta) %*% diag(1 / sigma_x ^ 2)
Bx <- B %*% Ax
Bx_t <- B_t %*% Ax
## effects of covariates
B_z <- matrix(runif(Q * p, 0.1, 0.2), Q, p)
## generate random intercept
xi <- MASS::mvrnorm(n, rep(0, Q), diag(1, Q, Q))
## generate observation error
E <- MASS::mvrnorm(n * T, rep(0, Q), diag(0.5, Q, Q))
## generate outcomes
Y <- eta[id, ] %*% t(B) +
  (eta[id, ] * tcrossprod(scale(time), rep(1, k))) %*% t(B_t) +
  X[id, ] %*% t(B_z) + xi[id,] + E

## for evaluation
trueh_pred <- eta_pred[id, ] %*% t(B) +
  (eta_pred[id, ] * tcrossprod(scale(time), rep(1, k))) %*% t(B_t)
Y_pred_oracle <- trueh_pred + X_pred[id, ] %*% t(B_z)
Y_pred <- Y_pred_oracle +
  MASS::mvrnorm(n, rep(0, Q), diag(1, Q, Q))[id,] +
  MASS::mvrnorm(n * T, rep(0, Q), diag(0.5, Q, Q))

# Make U matrix
age = 1:T
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