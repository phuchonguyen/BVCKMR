##############################
# BVCKMR Simulation Run Code #
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
	n			    = 200
	QQ        = 5 # Number of outcomes
	res = list()
	models = list()
	for (i in 1:QQ) {
	  start_time <- Sys.time()
	  # source("BVCKMR_Sim_Data.R")
	  source("famr_sim_data_01.R")
	  num.reps 	= 20000
	  sel 		  = 10000:num.reps
	  #cat(c("Job started at:",date()),fill=TRUE)
	  source("BVCKMR_Sim_Initialize.R")
	  source("BVCKMR_MCMC.R")
	  #cat(c("Job finished at:",date()),fill=TRUE)
	  MCMC = list(Sigsq = sigsq.post, Lam1 = lambda1.post, H = h.post, Tau = tausq.post, Beta = beta.post, DINV = D.inv.post, b = b.post, Z = Z, Y = Y, X = X, U = U, W = W, sel = sel, cofactor = cofactor, q = q, M = M, n = n)  
	  # models[[i]] = MCMC
	  # save(models, file=paste0("sim-gp-output/scenario-001/MCMC-", RANDOM_SEED, ".RData"))
	  
	  # Analysis
	  # Summary.h 	= matrix(NA, nrow=q, ncol=4)
	  # true.h1 = 0.5*(Z[,1]^2 + Z[,2]^2 - Z[,4]^2 - Z[,5]^2 + 0.5*Z[,4]*Z[,5] + Z[,4] + Z[,5])
	  # true.h2 = 0.5*(Z[,1]^2 - Z[,2]^2 - Z[,1]*Z[,2] + Z[,3]^2 + Z[,4] - Z[,5])
	  
	  # Posterior mean of estimation functions at new Z
	  hpred.fit = newh.postmean(Znew = Z_pred, sel = sel)
	  hpred.h1 = hpred.fit$postmean[(n+1):(2*n)]
	  hpred.h2 = hpred.fit$postmean[1:n]
	  hpred = rep(hpred.h1, each=T) + 
	    matrix(rep(hpred.h2, each=T)) * rep(age, n)
	  ypred = Z_pred[id, ] %*% colMeans(MCMC$Beta) + hpred
	    
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
	  # Rank metals based on mean importance
	  predh.rank = rep(0, M)
	  h.rank = rep(0, M)
	  for (mm in 1:M) {
	    predh.rank[mm] = sum(abs(mat[1,mm] + mat[2,mm]*age))
	    h.rank[mm] = sum(abs(Bx[i, mm] + Bx_t[i, mm]*age))
	  }
	  predh.rank = rank(predh.rank, ties.method = "max")
	  h.rank = rank(h.rank, ties.method = "max")
	  
	  res[[i]] = list("pmse" = sum( (ypred - Y_pred)^2) / (n*T),
	                  "pmse_mean" = sum((mean(Y) - Y_pred)^2) / (n*T),
	                  "pmse_true" = sum( (Y_pred_oracle - Y_pred)^2) / (n*T),
	                  "fnorm2_rank_importance" = norm(predh.rank-h.rank, type="F")^2,
	                  "spearman_rank_importance" = cor(predh.rank, h.rank, method="spearman")
	                  )
	  save(res, file=paste0("sim-bvckmr-output/scenario-001/res-", RANDOM_SEED, ".RData"))
	  
	  end_time <- Sys.time()
	  cat('\n Model took: ', end_time - start_time)
	}
}
