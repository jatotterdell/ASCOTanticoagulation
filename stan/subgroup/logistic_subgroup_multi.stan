// Simple logistic regression model for Bernoulli outcome

data {
  int N;
  int K; // number of coefficients
  int Q; // number of interacting covariates
  int P; // number of terms per covariate
  matrix[N, K] X; // fixed effects
  matrix[N, P*Q] W; // interaction subgroup effects
  array[N] int y;
  vector[K] beta_sd;
}

parameters {
  vector[K] beta_raw;
  vector[Q, K] gamma_raw;
  cholesky_factor_corr[P] L;
  vector<lower=0>[P] omega;
}

transformed parameters {
  vector[K] beta = beta_sd .* beta_raw;
  matrix[Q, K] gamma = L * gamma_raw;
}

model {
  // lp
  vector[N] eta = X * beta + W * to_vector(gamma);
  // prior
  target += lkj_corr_cholesky_lpdf(L | 1);
  for q in 1:Q {
    target += gamma[q]
  }
  target += std_normal_lpdf(beta_raw) +
            std_normal_lpdf(gamma_raw) +
            lkj_corr_cholesky_lpdf(L | 1) +
            student_t_lpdf(omega | 3, 0, 1) - student_t_lccdf(0 | 3, 0, 1);
  // likelihood
  target += bernoulli_logit_lpmf(y | eta);
}
