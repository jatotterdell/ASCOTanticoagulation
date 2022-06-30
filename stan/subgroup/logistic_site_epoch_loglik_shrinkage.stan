data {
  int N;          // number of participants
  int K;          // number of design parameters
  int Kint;
  int M_region;   // number of regions
  int M_site;     // total number of sites across all regions
  int M_epoch;    // number of epochs
  array[N] int y; // outcome
  matrix[N, K] X; // design matrix, including treatment design, covariates, and eligibilities
  matrix[N, Kint] Xint; // interaction design terms
  vector[K] beta_sd; // prior for design coefficient parameters

  array[M_site] int<lower=1> region_by_site; // region indicator for each site
  array[N] int<lower=1> site;                // site indicator
  array[N] int<lower=1> epoch;               // epoch indicator
}

parameters {
  vector[K] beta_raw;                 // design coefficients
  vector[Kint] beta_int_raw;             // interaction terms
  vector[M_site] epsilon_site;        // site coefficients
  vector[M_epoch-1] epsilon_epoch;    // epoch coefficients, first term = 0.
  vector<lower=0>[M_region] tau_site; // region-specific site sd
  real<lower=0> tau_epoch;            // cohort sd
  real<lower=0> tau_int;
}

transformed parameters {
  vector[K] beta = beta_sd .* beta_raw;
  vector[Kint] beta_int = tau_int * beta_int_raw;
  vector[M_site] gamma_site;
  vector[M_epoch] gamma_epoch;
  vector[N] eta;
  // use country specific site variation
  for(m in 1:M_site) {
    gamma_site[m] = tau_site[region_by_site[m]] * epsilon_site[m];
  }
  // define random-walk(1) prior
  gamma_epoch[1] = 0.0;
  gamma_epoch[2:M_epoch] = tau_epoch * cumulative_sum(epsilon_epoch);
  eta = X * beta + Xint * beta_int + gamma_site[site] + gamma_epoch[epoch];
}

model {
  // prior
  target += std_normal_lpdf(beta_raw)
          + std_normal_lpdf(beta_int_raw)
          + std_normal_lpdf(epsilon_epoch)
          + std_normal_lpdf(epsilon_site)
          + student_t_lpdf(tau_int | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5)
          + student_t_lpdf(to_vector(tau_site) | 3, 0, 1) - M_region*student_t_lccdf(0 | 3, 0, 1)
          + student_t_lpdf(tau_epoch | 3, 0, 1) - student_t_lccdf(0 | 3, 0, 1);
  // likehoood
  target += bernoulli_logit_lpmf(y | eta);
}

generated quantities {
  array[N] int y_ppc;
  vector[N] log_lik;
  for(n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | eta[n]);
    y_ppc[n] = bernoulli_logit_rng(eta[n]);
  }
}
