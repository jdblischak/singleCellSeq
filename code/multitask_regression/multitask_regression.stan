data {
  int<lower=0> N; 
  int<lower=0> K; 
  int<lower=0> P;
  vector[N] y[K];
  matrix[N,P] x[K];
}

parameters {
  vector[P] mu; 
  cov_matrix[P] Sigma; 
  real<lower=0> noiseVar; 
}

model {
  for (k in 1:K) {
    y[k] ~ multi_normal( x[k] * mu, x[k] * Sigma * x[k]' + diag_matrix(rep_vector(noiseVar, N)));
  }
}
