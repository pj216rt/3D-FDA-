#loading in pervious files
#source("beta_generation.R")
source("data_gen.R")

#need to get the data that the STAN function needs into a list
#specify the number of knots
knots <- 5
subjects <- 50
replications <- 5
alph <- 0.01

compile.stan.dat <- function(gener.dat, knots, subjects, J, alpha, true_betas){
  #get the grid size for one dimension
  D <- ncol(gener.dat$obs)/3
  Kt = knots
  I = subjects
  J=replications
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

compiled.stan.data <- compile.stan.dat(gener.dat = generated.data, knots = knots, 
                                       subjects = subjects,
                         alpha = alph, true_betas = beta.coefficients)

#compiling the stan code for use
fit1  <- cmdstanr::cmdstan_model("stan_files/corr_rand_effects_model.stan")
results <- fit1$sample(
  data = compiled.stan.data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 10 # print update every 100 iters
)
