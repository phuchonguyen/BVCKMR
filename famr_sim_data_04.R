########################
# Data simulation code #
########################

library(MASS); library(stats)

#n 				= 100 # Number of subjects
q 				= 2 # Random slope and random intercept
p 				= 2 # Number of confounders
M				  = 10 # Number of metals
k         = 9 # Rank of covariance of metals
# age 			= c(-5:5)
T				  = 10 # Number of time points
age       = 1:T
N				  = n*T

#############################
# Simulating Z as multivariate normal matrix
Theta = matrix(rnorm(M*k), M, k)
Z     = matrix(mvrnorm(n = n, mu=c(rep(0,M)), 
                       Sigma = tcrossprod(Theta, Theta)), 
               byrow=FALSE, nrow=n)
Z_pred     = matrix(mvrnorm(n = n, mu=c(rep(0,M)), 
                       Sigma = tcrossprod(Theta, Theta)), 
               byrow=FALSE, nrow=n)

Z 			= apply(Z, 2, scale)
Z_pred 			= apply(Z_pred, 2, scale)

X				= cbind( 
  matrix( rep(rnorm(n, mean=0, sd=1), each=T), byrow=FALSE ), 
  matrix( rep(scale(matrix(sample(1:2, n, prob=c(0.4, 0.6), replace = TRUE))), each=T) )
)
X_pred				= cbind( 
  matrix( rep(rnorm(n, mean=0, sd=1), each=T), byrow=FALSE ), 
  matrix( rep(scale(matrix(sample(1:2, n, prob=c(0.4, 0.6), replace = TRUE))), each=T) )
)

beta		= runif(p, 0.5, 1)

b				= c( t(mvrnorm(n = n, mu=rep(0,q), Sigma = matrix(data = c(0.25,0.005,0.005,0.05), nrow=q))))
b_pred	= c( t(mvrnorm(n = n, mu=rep(0,q), Sigma = matrix(data = c(0.25,0.005,0.005,0.05), nrow=q))))
h1.coef = runif(5, 0.5, 1) * sample(c(-1, 1), size = 5, replace = TRUE)
h2.coef = runif(1, 0.25, 0.5)
res.sd = 1

# Make U matrix

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

Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + 
  X %*% beta + 
  rep(cbind(Z[,1]^2, - Z[,2]^2, Z[,1]*Z[,2], Z[,3], Z[,4]) %*% h1.coef, each=T) + 
  matrix(rep(h2.coef*(Z[,1]^2 - Z[,2]^2), each=T)) * rep(age, n) +
  U %*% b

# For validation
trueh_pred = rep(cbind(Z_pred[,1]^2, - Z_pred[,2]^2, Z_pred[,2]*Z_pred[,3], Z_pred[,2], Z_pred[,3]) %*% h1.coef, each=T) + 
  matrix(rep(h2.coef*(Z_pred[,1]^2 - Z_pred[,2]^2), each=T)) * rep(age, n)
Y_pred_oracle <- trueh_pred + X_pred %*% beta
Y_pred = matrix(rnorm(n=N, mean=0, sd=res.sd)) + 
  Y_pred_oracle +
  U %*% b_pred

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
