##############################
# BVCKMR Simulation Run Code #
# Scenario 3:
# Correlated exposure with 2 latent factors (same as Scenario 1 and 2)
# Quadratic outcome functions. 
# The same variables are important for all outcomes, but the effect sizes are 
# different.
# Independent random effects.
# Same as Scenario 4 except k=2 instead of k=9
##############################

source("BVCKMR_Sim_Pred.R")
source("install_package_temporarily.R")
RANDOM_SEED = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(RANDOM_SEED)

# Set tasks
doMCMC 			= TRUE
doginv			= FALSE
docofactor	= TRUE
cofactor		= 1e-5	# randomly generated kernel matrix is not always invertible.

if (doMCMC)
{
	n			    = 100
	Q         = 5 # Number of outcomes
	res = list()
	models = list()
	for (i in 1:Q) {
	  start_time <- Sys.time()
	  # source("BVCKMR_Sim_Data.R")
	  source("famr_sim_data_03.R")
	  num.reps 	= 10000
	  sel 		  = 5000:num.reps
	  #cat(c("Job started at:",date()),fill=TRUE)
	  source("BVCKMR_Sim_Initialize.R")
	  source("BVCKMR_MCMC.R")
	  #cat(c("Job finished at:",date()),fill=TRUE)
	  MCMC = list(Sigsq = sigsq.post, Lam1 = lambda1.post, H = h.post, Tau = tausq.post, Beta = beta.post, DINV = D.inv.post, b = b.post, Z = Z, Y = Y, X = X, U = U, W = W, sel = sel, cofactor = cofactor, q = q, M = M, n = n)  
	  # models[[i]] = MCMC
	  # save(models, file=paste0("sim-gp-output/scenario-004/MCMC-", RANDOM_SEED, ".RData"))
	  
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
	  # model.h = lm(hpred ~ trueh_pred)
	  
	  # Summary.h = list(summary(model.1)$coef[1,1], 
	  #                   summary(model.1)$coef[2,1], 
	  #                   summary(model.1)$r.sq, 
	  #                   sqrt( sum( (hpred.h1 - true.h1)^2) / n))
	  source("BVCKMR_Sim_RelativeImportance.R")
	  qs = c(0.25, 0.75, 0.50)
	  true.mat = rep(NA, M)
	  for (j in 1:M) {
	    # Get true.mat
	    cross.sec 		= rbind(apply(Z, 2, median), apply(Z, 2, median))
	    cross.sec[,j] 	= c(quantile(Z[,j], qs[2]), quantile(Z[,j], qs[1]))
	    h1.tmp = cbind(cross.sec[,1]*cross.sec[,2], 
	                   0.5*cross.sec[,7]^2, cross.sec[,8]) %*% h1.coef
	    h2.tmp = cbind(cross.sec[,1]^2, -cross.sec[,2]^2) %*% h2.coef # The rank of effect at all time points are the same
	    true.mat[j] = (h1.tmp[1] + h2.tmp[1]) - (h1.tmp[2] + h2.tmp[2])
	  }
	  # Rank metals based on mean importance
	  predh.rank = rank(abs(mat[1,] + mat[2,]), ties.method = "max")
	  h.rank = rank(abs(true.mat), ties.method = "max")
	  
	  res[[i]] = list("pmse" = sum( (ypred - Y_pred)^2) / (n*T),
	                  "pmse_mean" = sum((mean(Y) - Y_pred)^2) / (n*T),
	                  "pmse_true" = sum( (Y_pred_oracle - Y_pred)^2) / (n*T),
	                  "fnorm2_rank_importance" = norm(predh.rank-h.rank, type="F")^2,
	                  "spearman_rank_importance" = cor(predh.rank, h.rank, method="spearman")
	                  )
	  save(res, file=paste0("sim-bvckmr-output/scenario-003/res-", RANDOM_SEED, ".RData"))
	  
	  end_time <- Sys.time()
	  cat('\n Model took: ', end_time - start_time)
	}
}
