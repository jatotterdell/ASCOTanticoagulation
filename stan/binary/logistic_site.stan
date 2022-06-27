data {
  int N;          // number of participants
  int K;          // number of design parameters
  int M_region;   // number of regions
  int M_site;     // total number of sites across all regions
  array[N] int y; // outcome
  matrix[N, K] X; // design matrix, including intercept, treatment design, covariates, eligibilities, region
  vector[K] beta_sd; // prior for design coefficient parameters
  array[M_site] int<lower=1> region_by_site; // region indicator for each site
  array[N] int<lower=1> site;                // site indicator
}

parameters {
  vector[K] beta_raw;                 // design coefficients
  vector[M_site] epsilon_site;        // site coefficients
  vector<lower=0>[M_region] tau_site; // region-specific site sd
}

transformed parameters {
  vector[K] beta = beta_sd .* beta_raw;
  vector[M_site] gamma_site;
  // use country specific site variation
  for(m in 1:M_site) {
    gamma_site[m] = tau_site[region_by_site[m]] * epsilon_site[m];
  }
}

model {
  // lp
  vector[N] eta = X * beta + gamma_site[site];
  // prior
  target += std_normal_lpdf(beta_raw)
          + std_normal_lpdf(epsilon_site)
          + student_t_lpdf(to_vector(tau_site) | 3, 0, 1) - M_region*student_t_lccdf(0 | 3, 0, 1);
  // likehoood
  target += bernoulli_logit_lpmf(y | eta);
}

generated quantities {
  array[N] int y_ppc;
  vector[N] eta = X * beta + gamma_site[site];
  for(n in 1:N) {
    y_ppc[n] = bernoulli_logit_rng(eta[n]);
  }
}
