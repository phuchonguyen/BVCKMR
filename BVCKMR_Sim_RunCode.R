##############################
# BVCKMR Simulation Run Code #
##############################

source("BVCKMR_Sim_Pred.R")
RANDOM_SEED = 1
set.seed(RANDOM_SEED)

# Set tasks
doMCMC 			= TRUE
doginv			= FALSE
docofactor	= TRUE
cofactor		= 1e-5	# randomly generated kernel matrix is not always invertible.

if (doMCMC)
{
	n			    = 100
	N         = 5 # Number of outcomes
	res = list()
	models = list()
	for (i in 1:N) {
	  start_time <- Sys.time()
	  # source("BVCKMR_Sim_Data.R")
	  source("famr_sim_data_04.R")
	  num.reps 	= 10000
	  sel 		  = 5000:num.reps
	  #cat(c("Job started at:",date()),fill=TRUE)
	  source("BVCKMR_Sim_Initialize.R")
	  source("BVCKMR_MCMC.R")
	  #cat(c("Job finished at:",date()),fill=TRUE)
	  MCMC = list(Sigsq = sigsq.post, Lam1 = lambda1.post, H = h.post, Tau = tausq.post, Beta = beta.post, DINV = D.inv.post, b = b.post, Z = Z, Y = Y, X = X, U = U, W = W, sel = sel, cofactor = cofactor, q = q, M = M, n = n)  
	  models[[i]] = MCMC
	  save(models, file=paste0("output/sim-04-MCMC", RANDOM_SEED, ".RData"))
	  
	  # Analysis
	  
	  Summary.h 	= matrix(NA, nrow=q, ncol=4)
	  
	  # true.h1 = 0.5*(Z[,1]^2 + Z[,2]^2 - Z[,4]^2 - Z[,5]^2 + 0.5*Z[,4]*Z[,5] + Z[,4] + Z[,5])
	  # true.h2 = 0.5*(Z[,1]^2 - Z[,2]^2 - Z[,1]*Z[,2] + Z[,3]^2 + Z[,4] - Z[,5])
	  
	  # Posterior mean of estimation functions at new Z
	  hpred.fit = newh.postmean(Znew = Z_pred, sel = sel)
	  hpred.h1 = hpred.fit$postmean[(n+1):(2*n)]
	  hpred.h2 = hpred.fit$postmean[1:n]
	  hpred = rep(hpred.h1, each=T) + 
	    matrix(rep(hpred.h2, each=T)) * rep(age, n)
	  ypred = X_pred %*% colMeans(MCMC$Beta) + hpred
	    
	  # Get fitted functions
	  # hhat.1	= apply(MCMC$H[sel, 1:n], 2, mean)
	  # hhat.2	= apply(MCMC$H[sel, (n+1):(2*n)], 2, mean)
	  
	  # Get predictive MSE
	  # model.1 = lm(hpred.h1~true.h1)
	  # model.2 = lm(hpred.h2~true.h2)
	  model.h = lm(hpred ~ trueh_pred)
	  
	  # Summary.h = list(summary(model.1)$coef[1,1], 
	  #                   summary(model.1)$coef[2,1], 
	  #                   summary(model.1)$r.sq, 
	  #                   sqrt( sum( (hpred.h1 - true.h1)^2) / n))
	  source("BVCKMR_Sim_RelativeImportance.R")
	  
	  res[[i]] = list("pmse" = sum( (ypred - Y_pred)^2) / (n*T),
	                  "lm_intercept" = summary(model.h)$coef[1,1], 
	                  "lm_slope" = summary(model.h)$coef[2,1], 
	                  "lm_r2" = summary(model.h)$r.sq,
	                  "fnorm_rank_importance" = norm(predh.rank-h.rank, type="F")
	                  )
	  save(res, file=paste0("output/sim-04-res", RANDOM_SEED, ".RData"))
	  
	  end_time <- Sys.time()
	  cat('\n Model took: ', end_time - start_time)
	}
}

beepr::beep()
