//STAN model with correlated random effects
//function to compute kronecker product
functions {
  matrix kronecker_prod(matrix A, matrix B) {
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
    for (j in 1:n) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start = (j - 1) * q + 1;
      col_end = (j - 1) * q + q;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
}

  real matMean(matrix A){
    real sums;
    sums = 0;
    for (i in 1:rows(A)){
      for (j in 1:cols(A)){
        sums = sums + A[i,j];
      }
    }
    return(sums/(cols(A)*rows(A)));
  }
}

data {
	int<lower=0> I;                  // number of subjects
	int<lower=0> J;                  // number of visits per subject
	int<lower=0> IJ;                 // total number of subjects
	int<lower=0> D;                  // grid length
	int<lower=0> p;                  // number of parameters

	int Kt;                          // number of spline basis functions
	
	//hyperparameters
	//vector<lower=0> [3] Az;
	//vector<lower=0> [3] Bz;
	
	// real<lower=0> Aw;
	// real<lower=0> Aw1;
	// real<lower=0> Aw2;
	// vector<lower=0> [p] Bw;
	// vector<lower=0> [p] Bw1;
	// vector<lower=0> [p] Bw2; 

	matrix [IJ, D] Y;                // outcome matrix  
	matrix [I, p] X;                //fixed effects design matrix
  matrix [IJ, I] Z;                //random effects design matrix
	
	matrix[D,3*Kt] THETA;                 // B-spline evaluation matrix
	cov_matrix[Kt] PenMatInv;        //inverse of the penalty matrix
	//matrix[3,3] Ident2;             //2*2 identity matrix
	
	matrix[IJ, D] Y_true;
	matrix[p, D] beta_true;
}

parameters {
  //matrices
	matrix[3*Kt, I] B_Z;               // matrix of fixed effect spline coefficients
	matrix[3*Kt, p] B_W;               //matrix of random effect spline coefficients
	
	//covariances
	cov_matrix [3] Omega_z; //covariance matrix for the B_Z component
	//array holding the cov matrices for each 
	array [p] cov_matrix [3] Omega_W;
	real<lower=0> lev1_sigma;
}

transformed parameters {
  //matrices
  cov_matrix[3*Kt] B_Z_cov_mat;
  array[p] cov_matrix[3*Kt] B_W_cov_mat;
  
  //constructing the covariance structures
  B_Z_cov_mat = kronecker_prod(Omega_z, PenMatInv);
  for(i in 1:p){
    B_W_cov_mat[i] = kronecker_prod(Omega_W[i], PenMatInv);
  }
  
  //beta matrix
  matrix[p, D] BETA;
  BETA = (THETA*B_W)';
}

model {
  //priors
	for (k in 1:p){
	  //v > p+1
	  Omega_W[k] ~ inv_wishart(4, diag_matrix(rep_vector(1, 3)));
	}

	for(k in 1:p){
	  B_W[,k] ~ multi_normal(rep_vector(0, 3*Kt), B_W_cov_mat[k]);
	}

	Omega_z ~ inv_wishart(4, diag_matrix(rep_vector(1, 3)));

	for(i in 1:I){
	  B_Z[, i] ~ multi_normal(B_W*X[i]', B_Z_cov_mat);
  }
	
	//level 1 variance
  lev1_sigma ~ inv_gamma(0.01, 0.01);
  
  for(n in 1:IJ){
	  Y[n] ~ normal(Z[n]*B_Z'*THETA', lev1_sigma);
	}
}

generated quantities {
  vector [p] MSE_BETA;
  for (i in 1:p){
    MSE_BETA[i] = dot_self((beta_true[i, ]- BETA[i,]));
  }
  
  matrix [IJ, D] y_hat;
  y_hat = Z*B_Z'*THETA';
  
  real MSE_Y;
  MSE_Y = matMean(crossprod(Y_true-y_hat));
}
