data {
  int N;          // number of participants
  int K;          // number of design parameters
  int M_region;   // number of regions
  int M_site;     // total number of sites across all regions
  int M_epoch;    // number of epochs
  array[N] int y; // outcome
  matrix[N, K] X; // design matrix, including treatment design, covariates, and eligibilities
  vector[K] beta_sd; // prior for design coefficient parameters
  array[N] int<lower=1> epoch; // epoch indicator
}

parameters {
  vector[K] beta_raw;                 // design coefficients
  vector[M_epoch-1] epsilon_epoch;    // epoch coefficients, first term = 0.
  real<lower=0> tau_epoch;            // cohort sd
}

transformed parameters {
  vector[K] beta = beta_sd .* beta_raw;
  vector[M_epoch] gamma_epoch;
  // define random-walk(1) prior
  gamma_epoch[1] = 0.0;
  gamma_epoch[2:M_epoch] = tau_epoch * cumulative_sum(epsilon_epoch);
}

model {
  // lp
  vector[N] eta = X * beta + gamma_epoch[epoch];
  // prior
  target += std_normal_lpdf(beta_raw)
          + std_normal_lpdf(epsilon_epoch)
          + student_t_lpdf(tau_epoch | 3, 0, 1) - student_t_lccdf(0 | 3, 0, 1);
  // likehoood
  target += bernoulli_logit_lpmf(y | eta);
}

generated quantities {
  array[N] int y_ppc;
  vector[N] eta =
    X * beta + gamma_epoch[epoch];
  for(n in 1:N) {
    y_ppc[n] = bernoulli_logit_rng(eta[n]);
  }
}
