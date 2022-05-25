data {
  int N;          // number of participants
  int K;          // number of design parameters
  int M_region;   // number of regions
  int M_site;     // total number of sites across all regions
  array[N] int y; // outcome
  matrix[N, K] X; // design matrix, including intercept, treatment design, covariates, and eligibilities, excluding region and site
  vector[K] beta_sd; // prior for design coefficient parameters

  array[N] int<lower=1> region;              // region indicator for individual
  array[M_site] int<lower=1> region_by_site; // region indicator for each site
  array[N] int<lower=1> site;                // site indicator
}

parameters {
  vector[K] beta_raw;                 // design coefficients
  vector[M_region-1] beta_region_raw; // region effect, first term = 0.
  vector[M_site] epsilon_site;        // site coefficients
  vector<lower=0>[M_region] tau_site; // region-specific site sd
}

transformed parameters {
  vector[K] beta = beta_sd .* beta_raw;
  vector[M_region] beta_region;
  vector[M_site] gamma_site;
  // use country specific site variation
  for(m in 1:M_site) {
    gamma_site[m] = tau_site[region_by_site[m]] * epsilon_site[m];
  }
  beta_region[1] = 0.0;
  beta_region[2:M_region] = beta_region_raw;
}

model {
  // lp
  vector[N] eta =
    X * beta + beta_region[region] + gamma_site[site];
  // prior
  target += std_normal_lpdf(beta_raw)
          + std_normal_lpdf(beta_region)
          + std_normal_lpdf(epsilon_site)
          + student_t_lpdf(to_vector(tau_site) | 3, 0, 1) - M_region*student_t_lccdf(0 | 3, 0, 1);
  // likehoood
  target += bernoulli_logit_lpmf(y | eta);
}
