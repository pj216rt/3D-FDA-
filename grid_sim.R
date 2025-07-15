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

#functions
#generate smooth beta functions
beta.gen.fun.tri <- function(grid.size, knots, parameters, deg){
  x <- seq(0, grid.size, length.out = grid.size)
  
  #b spline basis matrix
  bs_basis <- bs(x, degree = deg)  
  
  beta_list <- vector("list", 3)
  
  #going to loop over each dimension
  for (dim in 1:3) {
    random_coeffs_matrix <- matrix(rnorm(ncol(bs_basis)*parameters), 
                                   nrow = ncol(bs_basis), ncol = parameters)
    
    #compute the beta matrix for dimension dim
    beta_matrix <- bs_basis %*% random_coeffs_matrix
    
    #store this in the list
    beta_list[[dim]] <- t(beta_matrix)
  }
  
  #combine all three matrices column wise
  coefs <- do.call(cbind, beta_list)
  
  return(coefs)
}

#generating 3 dimenaional data, assuming balanced design
gen.3d.data <- function(num.subj = 50, num.visits = 5,
                        beta.vals, var.x = 1, var.z=1){
  #make sure you supply beta vals
  if (missing(beta.vals)) {
    message("Beta values not supplied.")
  }
  
  grid.size.3d <- ncol(beta.vals)
  n.params <-nrow(beta.vals)
  
  #total observations
  tot.obs <- num.subj*num.visits
  
  #generate the X design matrix
  X <- matrix(rnorm(n=(num.subj*n.params), sd = sqrt(var.x)), 
              nrow = num.subj, ncol = n.params) %>%
    as.data.frame() %>%
    mutate(subj = row_number())
  X$subj = as.factor(X$subj)
  X1 <- X[rep(seq_len(nrow(X)), each = num.visits), ]
  X.des <- model.matrix(~. -1 -subj, data = X1)
  
  #generate fixed effects
  fixef = as.matrix(X.des) %*% as.matrix(beta.vals)
  
  #random effects
  Z.des = model.matrix( ~ 0 + subj + (-1):subj, data = X1)
  subj.ranef <- matrix(rnorm(n=(num.subj*grid.size.3d), sd=sqrt(var.z)), 
                       nrow = num.subj, ncol = grid.size.3d)
  ranef <- Z.des %*% subj.ranef
  
  #level 1 residuals
  eps <- rnorm(n = tot.obs)
  
  Yij.true <- fixef + ranef
  Yij.obs <- fixef + ranef + eps
  
  #return "observed" data and other stuff in a list
  
  output.list <- list()
  output.list$obs <- Yij.obs
  output.list$true <- Yij.true
  output.list$x_design <- X.des
  output.list$z_design <- Z.des
  output.list$raw.data <- X1
  return(output.list)
  
}
#function to compile data for STAN
compile.stan.dat <- function(gener.dat, knots, subjects, 
                             J, alpha, true_betas){
  #get the grid size for one dimension
  D <- ncol(gener.dat$obs)/3
  Kt = knots
  I = subjects
  J=5
  p = ncol(gener.dat$raw.data)-1 #remove the subject identifier
  IJ=I*J
  
  #spline function matrix from Goldsmith, then kronecker multiply this by I_3
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)
  Gamma = kronecker(diag(1,3), Theta)
  
  #constructing the penalty matrix
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2
  P.mat.inv <- solve(P.mat) #get the inverse of the penalty matrix
  
  SUBJ = factor(apply(gener.dat$z_design %*% 1:dim(gener.dat$z_design)[2], 1, sum))
  
  ## find first observation for each subject.  Only need baseline
  firstobs = rep(NA, length(unique(SUBJ)))
  for(i in 1:length(unique(SUBJ))){
    firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
  }
  Wi = gener.dat$x_design[firstobs,]
  
  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(gener.dat$obs))
  IIP = kronecker(kronecker(diag(1, I, I), P.mat), diag(1,3)) #kronecker this by a 3*3 matrix
  WIk = kronecker(Wi, diag(1, 3*Kt, 3*Kt))
  tWIW = t(WIk) %*% IIP %*% WIk
  tWI = t(WIk) %*% IIP
  
  # initial estimation and hyperparameter choice
  vec.bz = solve(kronecker(t(gener.dat$z_design)%*% gener.dat$z_design, t(Gamma) %*% Gamma)) %*% t(kronecker(gener.dat$z_design, Gamma)) %*% Y.vec
  bz = matrix(vec.bz, nrow = 3*Kt, ncol = I)
  
  w.temp = kronecker(t(Wi), diag(1, Kt, Kt))
  vec.bw = solve(tWIW) %*% tWI %*% vec.bz
  bw = matrix(vec.bw, nrow = 3*Kt, ncol = p)
  
  Yhat = as.matrix(gener.dat$z_design %*% t(bz) %*% t(Gamma))
  varhat = var(as.vector(gener.dat$obs - Yhat))
  
  Psi = diag(varhat*IJ, 3*D, 3*D)
  v = IJ
  inv.sig = solve(Psi/v)
  
  #set initial values
  Az = I*Kt / 2  #set Az equal to Az1 equal to Az2
  Az1 = I*Kt / 2 
  Az2 = I*Kt / 2
  Bz = sum(diag((t(bz[1:Kt, ]) - Wi %*% t(bw[1:Kt, ])) %*% P.mat %*% t(t(bz[1:Kt, ]) - Wi %*% t(bw[1:Kt, ]))))
  Bz1 = sum(diag((t(bz[(Kt+1):(2*Kt), ]) - Wi %*% t(bw[(Kt+1):(2*Kt), ])) %*% P.mat %*% t(t(bz[(Kt+1):(2*Kt), ]) - Wi %*% t(bw[(Kt+1):(2*Kt), ]))))
  Bz2 = sum(diag((t(bz[((2*Kt)+1):(3*Kt), ]) - Wi %*% t(bw[((2*Kt)+1):(3*Kt), ])) %*% P.mat %*% t(t(bz[((2*Kt)+1):(3*Kt), ]) - Wi %*% t(bw[((2*Kt)+1):(3*Kt), ]))))
  
  Aw = Kt / 2  #set Aw equal to Aw1
  Aw1 = Kt / 2
  Aw2 = Kt / 2
  Bw = sapply(1:p, function(u) max(1, sum(diag( t(bw[1:Kt,u]) %*% P.mat %*% (bw[1:Kt,u])))))
  Bw1 = sapply(1:p, function(u) max(1, sum(diag( t(bw[(Kt+1):(2*Kt),u]) %*% P.mat %*% (bw[(Kt+1):(2*Kt),u])))))
  Bw2 = sapply(1:p, function(u) max(1, sum(diag( t(bw[((2*Kt)+1):(3*Kt),u]) %*% P.mat %*% (bw[((2*Kt)+1):(3*Kt),u])))))
  
  #compile this stuff into a list to feed into STAN
  stan.dat <- list(I = I, J = J, IJ = IJ, D = 3*D, p = ncol(gener.dat$x_design), 
                   Kt = Kt, Y = gener.dat$obs, X = Wi,
                   Z = gener.dat$z_design, THETA = Gamma, PenMatInv = P.mat.inv,
                   Y_true=gener.dat$true, beta_true=true_betas)
  
  return(stan.dat)
}

#function to extract beta values
#want a function to extract the values of various values
fixed.effects.analysis <- function(stan_object, ci.level=0.95){
  #specifcying the alpha for the given CI level
  alpha <- (1 - ci.level) / 2
  BETA.results <- stan_object$summary(
    variables = "BETA",
    mean,
    ~quantile(.x, probs = c(alpha, 1 - alpha))
  )
  
  #rename columns to clarify
  colnames(BETA.results) <- c("variable", "mean", "CI_lower", "CI_upper")
  BETA.results <- BETA.results %>%
    mutate(
      parameter = as.integer(gsub("BETA\\[([0-9]+),([0-9]+)\\]", "\\1", variable)),  
      grid_point = as.integer(gsub("BETA\\[([0-9]+),([0-9]+)\\]", "\\2", variable)) 
    )
  
  n_predictors <- n_distinct(BETA.results$parameter)
  n_dims <- 3 #three dimensional
  n_per_dim <- n_distinct(BETA.results$grid_point)/n_dims
  print(n_per_dim)
  
  beta_summary <- BETA.results %>%
    mutate(
      predictor = parameter,
      dim = (grid_point - 1) %/% n_per_dim + 1,
      position = (grid_point - 1) %% n_per_dim + 1) %>%
    dplyr::select(predictor, dim, position, mean, CI_lower, CI_upper)
  
  return(beta_summary)
}

#coverage calculation:
compute.coverage <- function(true_values, summary.dat, grid.size){
  #pivot true values
  true_values_long <- true_values %>%
    as.data.frame() %>%
    mutate(predictor = 1:nrow(true_values)) %>% 
    pivot_longer(cols = starts_with("V"), names_to = "position", 
                 values_to = "true_value") %>%
    mutate(position = as.integer(gsub("V", "", position)),  
           dim = (position - 1) %/% grid.size + 1) %>%
    group_by(dim, predictor, position)
  
  coverage_data <- summary.dat %>%
    inner_join(true_values_long, by = c("predictor", "position", "dim")) %>%
    mutate(
      #Check if the true value lies within the CI bounds
      coverage = ifelse(true_value >= CI_lower & true_value <= CI_upper, 1, 0)
    )
  
  #Compute the coverage rate for each predictor
  coverage_summary <- coverage_data %>%
    group_by(predictor) %>%
    summarise(coverage_rate = mean(coverage))
  
  return(coverage_summary)
}

#MSE calculation
compute.mse.y <- function(stan_object){
  #compute mean and 95% CI
  MSE_Y_summary <- stan_object$summary(
    variables = "MSE_Y",
    mean,
    ~quantile(.x, probs = c(0.025, 0.975))
  )
  colnames(MSE_Y_summary) <- c("variable", "mean", "CI_lower", "CI_upper")
  
  return(MSE_Y_summary)
}

#compute integrated mean squared error
#function to compute integrated mean squared error
compute.IMSE <- function(summary.dat, true_values, grid.size) {
  #pivot true values
  true_values_long <- true_values %>%
    as.data.frame() %>%
    mutate(predictor = 1:nrow(true_values)) %>% 
    pivot_longer(cols = starts_with("V"), names_to = "position", 
                 values_to = "true_value") %>%
    mutate(position = as.integer(gsub("V", "", position)),  
           dim = (position - 1) %/% grid.size + 1) %>%
    group_by(dim, predictor, position) %>% group_by(dim, position, predictor)
  
  # #join the data to compute squared difference, then pivot wider so that we
  # #have a p*3D matrix of the squared differences
  MSE.dat <- summary.dat %>%
    inner_join(true_values_long, by = c("predictor", "position", "dim")) %>% 
    mutate(squared_diff = (mean - true_value)^2)
  
  #Calculate the IMSe for each predictor
  IMSE_vals <- MSE.dat %>% group_by(predictor) %>%
    summarise(IMSE = mean(squared_diff))
  
  
  return(IMSE_vals)
}

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
# multichain.golds.gibbs <- function(n.chains = 4, Y, fixef.form, id, data, Kt,
#                                    N.iter=1000, N.burn=200, alpha=0.1){
#   #plan for parallel execution
#   plan(multisession, workers = n.chains)
#   
#   chains_results <- future_lapply(1:n.chains, function(chain_id) {
#     set.seed(chain_id)
#     result <- golds.gibbs.three.d(Y, fixef.form, id, data, Kt, N.iter, N.burn, alpha)
#     cat("Finished chain", chain_id, "\n")
#     return(result)
#   }, future.seed = TRUE) #independent random seeds across chains
#   
#   #reset to sequential execution
#   plan(sequential)
#   
#   #get average of output across chains
#   matrix_names <- names(chains_results[[1]])
#   
#   averaged_matrices <- lapply(matrix_names, function(name) {
#     matrices <- lapply(chains_results, function(lst) lst[[name]])  #Extract matrices
#     Reduce("+", matrices) / length(matrices)  #Element-wise average
#   })
#   names(averaged_matrices) <- matrix_names
#   
#   return(averaged_matrices)
# }


#plotting example of simulated data
test <- beta.gen.fun.tri(grid.size = 25, knots = 5, parameters = 5, deg = 3)

cols_per_dim <- ncol(test) /3
dim_labels <- c("X", "Y", "Z")

df <- as.data.frame(test) %>%
  mutate(param = row_number()) %>%
  pivot_longer(
    cols = -param,
    names_to = "colname",
    values_to = "value"
  ) %>%
  mutate(
    col_index = as.integer(gsub("V", "", colname)),
    dim = dim_labels[ceiling(col_index / cols_per_dim)],
    time = (col_index - 1) %% cols_per_dim + 1
  )

ggplot(df, aes(x = time, y = value, group = param)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~dim) +
  labs(x = "Time", y = "Coefficient", title = "Functional Coefficients by Dimension")


#generating the data
test1 <- gen.3d.data(beta.vals = test)

df <- as.data.frame(test1$obs) %>%
  mutate(subject = rep(1:50, each = 5))



#long format
df_long <- df %>%
  pivot_longer(
    cols = -subject,
    names_to = "time_var",
    values_to = "value"
  ) %>%
  mutate(
    col_index = as.numeric(gsub("V", "", time_var)),
    time = (col_index - 1) %% 25 + 1,
    dimension = c("X", "Y", "Z")[ceiling(col_index / 25)]
  ) %>%
  group_by(subject) %>%
  mutate(replication = rep(1:5, each = 75)) %>%
  ungroup()

#randomly select 3 subjects
highlight_subjects <- sample(unique(df_long$subject), 3)
df_long <- df_long %>%
  mutate(
    highlight = ifelse(subject %in% highlight_subjects, "highlight", "dimmed"),
    subject_label = ifelse(highlight == "highlight", as.character(subject), NA),
    color_val = ifelse(highlight == "highlight", as.character(subject), "dimmed"),
    alpha_val = ifelse(highlight == "highlight", 1, 0.1)
  )

#plotting with vidris, for colorblind friendly visualization
ggplot(df_long, aes(
  x = time,
  y = value,
  group = interaction(subject, replication),
  color = color_val,
  alpha = alpha_val
)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~dimension, nrow = 1) +
  scale_color_manual(
    values = c("dimmed" = "gray80", setNames(viridisLite::viridis(3), as.character(highlight_subjects))),
    breaks = as.character(highlight_subjects),
    name = "Subject"
  ) +
  scale_alpha_identity() +
  labs(
    x = "Time",
    y = "Observed Position",
    title = "Simulated Functional Data by Dimension (3 Subjects Highlighted)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5)
  )





#grid conditions
grid_sizes <- c(25)
knots_gen_values <- c(5)
knots_sample_values <- c(5, 10)
num_predictors <- c(3, 10)
num_subjects <- c(50)
var_fixed_effects <- c(0.1, 1)
var_random_effects <- c(0.1, 1)

#create grid
param_grid <- expand.grid(
  D = grid_sizes,
  Kt_gen = knots_gen_values,
  Kt_sample = knots_sample_values,
  p = num_predictors,
  I = num_subjects,
  var_fixed = var_fixed_effects,
  var_random = var_random_effects
)

#Simulation
results_list <- list()
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  #generate splines
  generated.beta <- beta.gen.fun.tri(grid.size = params$D, knots = params$Kt_gen, 
                                     deg = 3, parameters = params$p)
  
  #generate data
  generated.data <- gen.3d.data(num.subj = params$I, num.visits = 5, 
                                beta.vals = generated.beta, var.x = params$var_fixed,
                                var.z = params$var_random)
  
  #compile data for STAN
  compiled.stan.data <- compile.stan.dat(gener.dat = generated.data, 
                                           knots = params$Kt_sample, subjects = params$I,
                                           alpha=0.01, true_betas = generated.beta)
  
  #run STAN 
  fit  <- cmdstanr::cmdstan_model("stan_files/corr_rand_effects_model.stan")
  
  #track start time
  start_time <- Sys.time()
  
  results <- fit$sample(
    data = compiled.stan.data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 10 # print update every 100 iters
  )
  
  end_time <- Sys.time()
  stan_runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  #summarize the beta values
  summarized_beta_values <- fixed.effects.analysis(stan_object = results)
  
  #compute coverage, mse_y, and mse_beta
  coverage <- compute.coverage(true_values = as.data.frame(generated.beta), 
                           summary.dat = summarized_beta_values, grid.size = params$D)
  
  #IMSE
  IM_squared <- compute.IMSE(summary.dat = summarized_beta_values, 
                             true_values = generated.beta, grid.size = params$D)
  
  
  ###GOLDSMITH MODEL.  Going to be four chains run in parallel
  start_time <- Sys.time()
  golds.model <- golds.gibbs.three.d(Y=generated.data$obs, fixef.form = ~. -1 -subj,
                         id = generated.data$raw.data$subj,
                         data = generated.data$raw.data, Kt=params$Kt_sample, N.iter = 8000,
                         N.burn = 4000, alpha = 0.01)
  end_time <- Sys.time()
  gibbs_runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  #coverage
  gibbs.cover <- rowMeans(golds.model$beta.LB <= generated.beta & generated.beta <= golds.model$beta.UB)
  
  #gibbs IMSE
  gibbs.IMSE <- rowMeans((golds.model$beta.pm - generated.beta)^2)
  
  #store results
  results_list[[i]] <- list(
    params = params,
    # stan.fixed_effects = summarized_beta_values,
    stan.coverage.calc = coverage,
    stan.integrated_MSE = IM_squared,
    stan.run.time = stan_runtime,
    gibbs.coverage <- gibbs.cover,
    gibbs.integrated_MSE <- gibbs.IMSE,
    gibbs.run.time <- gibbs_runtime
  )
  
  
  #print indicator showing how far along grid run is
  print(paste("Completed iteration", i, "of", nrow(param_grid)))
}




##EXTRACTION STUFF###
#extracting the results of the list
params <- lapply(results_list, function(x) x[[1]])
params.df <- as.data.frame(do.call(rbind, params))

stan.run.times <- sapply(results_list, function(x) x[[4]])
params.df$stan_times <- stan.run.times

gibbs.run.times <- sapply(results_list, function(x) x[[7]])
params.df$gibbs_times <- gibbs.run.times

stan.coverage <- lapply(results_list, function(x) x[[2]])
stan.mean.cov <- sapply(stan.coverage, function(tbl) mean(tbl[[2]]))
params.df$mean.stan.coverages <- stan.mean.cov

gibbs.coverage <- lapply(results_list, function(x) x[[5]])
mean.gibbs.cover <- sapply(gibbs.coverage, mean)
params.df$mean.gibbs.coverages <- mean.gibbs.cover

stan.imse <- lapply(results_list, function(x) x[[3]])
stan.mean.cov <- sapply(stan.imse, function(tbl) mean(tbl[[2]]))
params.df$imse.stan <- stan.mean.cov

gibbs.imse <- lapply(results_list, function(x) x[[6]])
mean.gibbs.imse <- sapply(gibbs.imse, mean)
params.df$imse.gibbs <- mean.gibbs.imse

#save to csv
write.csv(params.df, "partial_results.csv")


str(params.df)
