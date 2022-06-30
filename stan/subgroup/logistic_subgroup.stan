// Simple logistic regression model for Bernoulli outcome

data {
  int N;
  int K;
  int Q;
  matrix[N, K] X; // fixed effects
  matrix[N, Q] W; // interaction subgroup effects
  array[N] int y;
  vector[K] beta_sd;
}

parameters {
  vector[K] beta_raw;
  vector[Q] gamma_raw;
  real<lower=0> omega;
}

transformed parameters {
  vector[K] beta = beta_sd .* beta_raw;
  vector[Q] gamma = omega * gamma_raw;
}

model {
  // lp
  vector[N] eta = X * beta + W * gamma;
  // prior
  target += std_normal_lpdf(beta_raw) +
            std_normal_lpdf(gamma_raw) +
            student_t_lpdf(omega | 3, 0, 1) - student_t_lccdf(0 | 3, 0, 1);
  // likelihood
  target += bernoulli_logit_lpmf(y | eta);
}
