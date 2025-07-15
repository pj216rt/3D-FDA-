library(tidyverse)
library(mvtnorm)
library(spam)  #not entirely sure what the spam call in the sampler does
library(splines)
library(MASS)
library(MCMCpack) #for the inverse wishart distribution
require(gridExtra)
library(bayesplot)
library(grid)
library(future.apply)

source("data_gen.R")

#Goldsmith's modifiyed Gibbs sampler
golds.gibbs.three.d <- function(Y, fixef.form, id, data, 
                                Kt, N.iter = 1000, N.burn = 200, alpha = .1){
  
  set.seed(1)
  
  ## fixed and random effect design matrices
  W.des = model.matrix( fixef.form, data = data)
  Z.des = model.matrix( ~ 0 + as.factor(id) + (-1):as.factor(id))
  W.des = as.spam(W.des)
  Z.des = as.spam(Z.des)
  
  I = dim(Z.des)[2]
  D = dim(Y)[2]/3 #divide by three here
  Ji = as.numeric(apply(Z.des, 2, sum))
  IJ = sum(Ji)
  p = dim(W.des)[2]
  
  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)
  Gamma = kronecker(diag(1,3), Theta)
  
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2
  # not doing anything to the penalty matrix here
  #need to wait until we declare the variance terms first
  
  SUBJ = factor(apply(Z.des %*% 1:dim(Z.des)[2], 1, sum))
  
  ## find first observation
  firstobs = rep(NA, length(unique(SUBJ)))
  for(i in 1:length(unique(SUBJ))){
    firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
  }
  Wi = W.des[firstobs,]
  
  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  IIP = kronecker(kronecker(diag(1, I, I), P.mat), diag(1,3)) #kronecker this by a 2*2 matrix
  WIk = kronecker(Wi, diag(1, 3*Kt, 3*Kt)) #double the knots
  tWIW = t(WIk) %*% IIP %*% WIk
  tWI = t(WIk) %*% IIP
  
  # initial estimation and hyperparameter choice
  vec.bz = solve(kronecker(t(Z.des)%*% Z.des, t(Gamma) %*% Gamma)) %*% t(kronecker(Z.des, Gamma)) %*% Y.vec
  bz = matrix(vec.bz, nrow = 3*Kt, ncol = I) #double the knots
  
  w.temp = kronecker(t(Wi), diag(1, Kt, Kt))
  vec.bw = solve(tWIW) %*% tWI %*% vec.bz
  bw = matrix(vec.bw, nrow = 3*Kt, ncol = p) #need to double the knots
  
  Yhat = as.matrix(Z.des %*% t(bz) %*% t(Gamma))
  varhat = var(as.vector(Y - Yhat))
  
  Psi = diag(varhat*IJ, 3*D, 3*D)
  v = IJ
  inv.sig = solve(Psi/v)
  
  Az = I*Kt / 2  #set Az equal to Az1
  Az1 = I*Kt / 2 
  Az2 = I*Kt / 2 
  Bz = sum(diag((t(bz[1:Kt, ]) - Wi %*% t(bw[1:Kt, ])) %*% P.mat %*% t(t(bz[1:Kt, ]) - Wi %*% t(bw[1:Kt, ]))))
  Bz1 = sum(diag((t(bz[(Kt+1):(2*Kt), ]) - Wi %*% t(bw[(Kt+1):(2*Kt), ])) %*% P.mat %*% t(t(bz[(Kt+1):(2*Kt), ]) - Wi %*% t(bw[(Kt+1):(2*Kt), ]))))
  Bz2 = sum(diag((t(bz[(2*Kt+1):(3*Kt), ]) - Wi %*% t(bw[(2*Kt+1):(3*Kt), ])) %*% P.mat %*% t(t(bz[(2*Kt+1):(3*Kt), ]) - Wi %*% t(bw[(2*Kt+1):(3*Kt), ]))))
  
  Aw = Kt / 2  #set Aw equal to Aw1
  Aw1 = Kt / 2
  Aw2 = Kt / 2
  Bw = sapply(1:p, function(u) max(1, sum(diag( t(bw[1:Kt,u]) %*% P.mat %*% (bw[1:Kt,u])))))
  Bw1 = sapply(1:p, function(u) max(1, sum(diag( t(bw[(Kt+1):(2*Kt),u]) %*% P.mat %*% (bw[(Kt+1):(2*Kt),u])))))
  Bw2 = sapply(1:p, function(u) max(1, sum(diag( t(bw[(2*Kt+1):(3*Kt),u]) %*% P.mat %*% (bw[(2*Kt+1):(3*Kt),u])))))
  
  
  ## matrices to store within-iteration estimates
  BW = array(NA, c(3*Kt, p, N.iter)) #change t0 2*Kt
  BW[,,1] = bw
  BZ = array(NA, c(3*Kt, I, N.iter)) #change to 2*Kt
  BZ[,,1] = bz
  INV.SIG = array(NA, c(3*D, 3*D, N.iter)) #change to 2D
  INV.SIG[,,1] = inv.sig
  LAMBDA.BW = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW[1,] = lambda.bw = Aw/Bw
  LAMBDA.BW.1 = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW.1[1,] = lambda.bw.1 = Aw1/Bw1
  LAMBDA.BW.2 = matrix(NA, nrow = N.iter, ncol = p)
  LAMBDA.BW.2[1,] = lambda.bw.2 = Aw2/Bw2
  LAMBDA.BZ = rep(NA, N.iter)
  LAMBDA.BZ[1] = lambda.ranef = Az/Bz
  LAMBDA.BZ.1 = rep(NA, N.iter)
  LAMBDA.BZ.1[1] = lambda.ranef.1 = Az1/Bz1
  LAMBDA.BZ.2 = rep(NA, N.iter)
  LAMBDA.BZ.2[1] = lambda.ranef.2 = Az2/Bz2
  
  y.post = array(NA, dim = c(IJ, 3*D, (N.iter - N.burn))) #change dimensions to 2*D
  
  cat("Beginning Sampler \n")
  pb <- txtProgressBar(min = 0, max = N.iter, initial = 0, style = 3)
  for(i in 1:N.iter){
    setTxtProgressBar(pb,i)
    # if(i %% 100 == 0){ #adding a print statement to tell us where we are
    #   print(i)
    # }
    
    #stick these outside of the subj loop?
    combined <- diag(c(lambda.ranef, lambda.ranef.1, lambda.ranef.2))
    two.P = kronecker(P.mat, combined)  #constructing P \otimes two variances
    
    #sigma_w_k.pre.kron <- c(LAMBDA.BW[i, ], LAMBDA.BW.1[i, ])
    sigma_w_k.pre.kron <- c(lambda.bw, lambda.bw.1, lambda.bw.2)
    #print(sigma_w_k.pre.kron)
    sigma_w_k.pre.kron.1 <- diag(sigma_w_k.pre.kron, nrow = length(sigma_w_k.pre.kron))
    sigma_w_k.post.kron <- kronecker(sigma_w_k.pre.kron.1, P.mat)
    #print(dim(sigma_w_k.post.kron))
    
    ###############################################################
    ## update b-spline parameters for subject random effects
    ###############################################################
    for(subj in 1:length(unique(SUBJ))){
      
      t.designmat.Z = t(kronecker(rep(1, Ji[subj]), Gamma))   #change to Gamma
      
      
      #print(dim(two.P))
      sigma = solve(t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% t(t.designmat.Z) +
                      two.P)
      mu = sigma %*% (t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% (as.vector(t(Y[which(SUBJ == unique(SUBJ)[subj]),]))) +
                        (two.P) %*% bw %*% t(Wi[subj,]))
      
      bz[,subj] = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = 3*Kt, ncol = 1)
    }
    ranef.cur = Z.des %*% t(bz) %*% t(Gamma)
    
    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################
    #error in evaluating the argument 'a' in selecting a method for function 
    #'solve': NAs in argument 5 and 'NAOK = FALSE' (dotCall64)
    
    #maybe here?
    #IIP.1 <- kronecker(kronecker(diag(1, I, I), P.mat), combined)
    IIP.1 <- kronecker(diag(1, I, I), two.P)
    tWIW.1 <- t(WIk) %*% IIP.1 %*% WIk
    
    tWI.1 = t(WIk) %*% IIP.1
    
    sigma = solve(tWIW.1 + sigma_w_k.post.kron)
    #print(determinant(sigma))
    mu = sigma %*% (tWI.1 %*% as.vector(bz))
    bw = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = 3*Kt, ncol = p)
    # 
    beta.cur = t(bw) %*% t(Gamma)
    
    ###############################################################
    ## update inverse covariance matrix
    ###############################################################
    
    resid.cur = Y - ranef.cur
    inv.sig = solve(riwish(v + IJ, Psi + t(resid.cur) %*% resid.cur))
    
    ###############################################################
    ## update variance components
    ###############################################################
    
    ## lambda for beta's
    for(term in 1:p){
      a.post = Aw + Kt/2
      a.post.1 = Aw1 + Kt/2
      a.post.2 = Aw2 + Kt/2
      b.post = Bw[term] + 1/2 * bw[(1:Kt),term] %*% P.mat %*% bw[(1:Kt),term]
      b.post.1 = Bw1[term] + 1/2 * bw[(Kt+1):(2*Kt),term] %*% P.mat %*% bw[(Kt+1):(2*Kt),term]
      b.post.2 = Bw1[term] + 1/2 * bw[(2*Kt+1):(3*Kt),term] %*% P.mat %*% bw[(2*Kt+1):(3*Kt),term]
      lambda.bw[term] = rgamma(1, a.post, b.post)
      lambda.bw.1[term] = rgamma(1, a.post.1, b.post.1)
      lambda.bw.2[term] = rgamma(1, a.post.2, b.post.2)
    }
    
    ## lambda for random effects
    a.post = Az + I*Kt/2
    a.post.1 = Az1 + I*Kt/2
    a.post.2 = Az2 + I*Kt/2
    b.post = Bz + .5 * sum(sapply(1:I, function(u) (t(bz[(1:Kt),u]) - Wi[u,] %*% t(bw[(1:Kt), ])) %*% P.mat %*% t(t(bz[(1:Kt),u]) - Wi[u,] %*% t(bw[(1:Kt), ])) ))
    b.post.1 = Bz1 + .5 * sum(sapply(1:I, function(u) (t(bz[(Kt+1):(2*Kt),u]) - Wi[u,] %*% t(bw[(Kt+1):(2*Kt), ])) %*% P.mat %*% t(t(bz[(Kt+1):(2*Kt),u]) - Wi[u,] %*% t(bw[(Kt+1):(2*Kt), ])) ))
    b.post.2 = Bz1 + .5 * sum(sapply(1:I, function(u) (t(bz[(2*Kt+1):(3*Kt),u]) - Wi[u,] %*% t(bw[(2*Kt+1):(3*Kt), ])) %*% P.mat %*% t(t(bz[(2*Kt+1):(3*Kt),u]) - Wi[u,] %*% t(bw[(2*Kt+1):(3*Kt), ])) ))
    lambda.ranef = rgamma(1, a.post, b.post)
    lambda.ranef.1 = rgamma(1, a.post.1, b.post.1)
    lambda.ranef.2 = rgamma(1, a.post.2, b.post.2)
    
    ###############################################################
    ## save this iteration's parameters
    ###############################################################
    
    BW[,,i] = as.matrix(bw)
    BZ[,,i] = as.matrix(bz)
    
    INV.SIG[,,i] = inv.sig
    LAMBDA.BW[i,] = lambda.bw
    LAMBDA.BW.1[i,] = lambda.bw.1
    LAMBDA.BW.2[i,] = lambda.bw.2
    LAMBDA.BZ[i] = lambda.ranef
    LAMBDA.BZ.1[i] = lambda.ranef.1
    LAMBDA.BZ.2[i] = lambda.ranef.2
    
    #print(lambda.bw)
    
    if(i > N.burn){
      y.post[,,i - N.burn] = ranef.cur
    }
  }
  close(pb)
  
  ###############################################################
  ## compute posteriors for this dataset
  ###############################################################
  
  ## main effects
  beta.post = array(NA, dim = c(p, 3*D, (N.iter - N.burn)))
  for(n in 1:(N.iter - N.burn)){
    beta.post[,,n] = t(BW[,, n + N.burn]) %*% t(Gamma)
  }
  
  beta.pm = apply(beta.post, c(1,2), mean)
  beta.LB = apply(beta.post, c(1,2), quantile, c(.025))
  beta.UB = apply(beta.post, c(1,2), quantile, c(.975))
  
  ## random effects
  b.pm = matrix(NA, nrow = I, ncol = 3*D)
  for(i in 1:I){
    b.post = matrix(NA, nrow = (N.iter - N.burn), ncol = 3*D)
    for(n in 1:(N.iter - N.burn)){
      b.post[n,] = BZ[,i, n + N.burn] %*% t(Gamma)
    }
    b.pm[i,] = apply(b.post, 2, mean)
  }
  
  ## covariance matrix
  sig.pm = solve(apply(INV.SIG, c(1,2), mean))
  
  ## export fitted values
  ranef.pm = Z.des %*% b.pm
  Yhat = apply(y.post, c(1,2), mean)
  
  ####
  ###mcmc diagnostics
  #b_w matrix
  # b.w.mcmc <- data.frame(matrix(BW, nrow=dim(BW)[3], byrow=TRUE))
  # b.w.mcmc <- mcmc(tail(b.w.mcmc, n = (N.iter - N.burn)))
  # print(effectiveSize(b.w.mcmc)[which.min(effectiveSize(b.w.mcmc))])
  # 
  # #b_z matrix
  # b.z.mcmc <- data.frame(matrix(BZ, nrow=dim(BZ)[3], byrow=TRUE))
  # b.z.mcmc <- mcmc(tail(b.z.mcmc, n = (N.iter - N.burn)))
  # print(effectiveSize(b.z.mcmc)[which.min(effectiveSize(b.z.mcmc))])
  # 
  # #lambda bz terms
  # lambda.bzs <- mcmc(tail(cbind(LAMBDA.BZ, LAMBDA.BZ.1, LAMBDA.BZ.2), n = (N.iter - N.burn)))
  # print(effectiveSize(lambda.bzs))
  # 
  # #lambda bw terms
  # lambda.bws <- mcmc(tail(cbind(LAMBDA.BW, LAMBDA.BW.1, LAMBDA.BW.2), n = (N.iter - N.burn)))
  # print(effectiveSize(lambda.bws))
  
  
  
  #returning what we want to return
  ret = list(beta.pm, beta.LB, beta.UB, ranef.pm, sig.pm, Yhat)
  names(ret) = c("beta.pm", "beta.LB", "beta.UB", "ranef.pm", "sig.pm", "Yhat")
  
  return(ret)
  
}

#multichain_golds.gibbs
multichain.golds.gibbs <- function(n.chains = 4, Y, fixef.form, id, data, Kt,
                                   N.iter=1000, N.burn=200, alpha=0.1){
  #plan for parallel execution
  plan(multisession, workers = n.chains)
  
  chains_results <- future_lapply(1:n.chains, function(chain_id) {
    set.seed(chain_id)
    result <- golds.gibbs.three.d(Y, fixef.form, id, data, Kt, N.iter, N.burn, alpha)
    cat("Finished chain", chain_id, "\n")
    return(result)
  }, future.seed = TRUE) #independent random seeds across chains
  
  #reset to sequential execution
  plan(sequential)
  
  #get average of output across chains
  matrix_names <- names(chains_results[[1]])
  
  averaged_matrices <- lapply(matrix_names, function(name) {
    matrices <- lapply(chains_results, function(lst) lst[[name]])  #Extract matrices
    Reduce("+", matrices) / length(matrices)  #Element-wise average
  })
  names(averaged_matrices) <- matrix_names
  
  return(averaged_matrices)
}


gold.test <- multichain.golds.gibbs(Y=generated.data$obs, fixef.form = ~. -1 -subj,
                                 id = generated.data$raw.data$subj,
                                 data = generated.data$raw.data, Kt=5, N.iter = 1000,
                                 N.burn = 200, alpha = 0.01)

#compute coverage
coverage <- rowMeans(gold.test$beta.LB <= beta.coefficients & beta.coefficients <= gold.test$beta.UB)

#IMSE
IMSE <- rowMeans((gold.test$beta.pm - beta.coefficients)^2)
