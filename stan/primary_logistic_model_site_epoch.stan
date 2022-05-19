data {
  int N;          // number of participants
  int K;          // number of treatment parameters
  int M_region;   // number of regions
  int M_site;     // total number of sites across all regions
  int M_epoch;    // number of epochs
  array[N] int y; // outcome
  matrix[N, K] X; // design matrix, including treatment design, covariates, and eligibilities
  vector[K] beta_sd; // prior for design coefficient parameters

  array[N] int<lower=1> region;              // region indicator for individual
  array[M_site] int<lower=1> region_by_site; // region indicator for each site
  array[N] int<lower=1> site;                // site indicator
  array[N] int<lower=1> epoch;               // epochepoch indicator
}

parameters {
  real beta_int_raw;                  // overall intercept
  vector[K] beta_raw;                 // design coefficients
  vector[M_region-1] beta_region_raw; // region effect, first term = 0.
  vector[M_site] epsilon_site;        // site coefficients
  vector[M_epoch-1] epsilon_epoch;    // epoch coefficients, first term = 0.
  vector<lower=0>[M_region] tau_site; // region-specific site sd
  real<lower=0> tau_epoch;            // cohort sd
}

transformed parameters {
  real beta_int = 2.5*beta_int_raw;
  vector[K] beta = beta_sd .* beta_raw;
  vector[M_region] beta_region;
  vector[M_site] gamma_site;
  vector[M_epoch] gamma_epoch;
  // use country specific site variation
  for(m in 1:M_site) {
    gamma_site[m] = tau_site[region_by_site[m]] * epsilon_site[m];
  }

  beta_region[1] = 0.0;
  beta_region[2:M_region] = beta_region_raw;
  gamma_epoch[1] = 0.0;
  gamma_epoch[2:M_epoch] = tau_epoch * cumulative_sum(epsilon_epoch);
}

model {
  // lp
  vector[N] eta =
    beta_int + X * beta + beta_region[region] +
    gamma_site[site] + gamma_epoch[epoch];
  // prior
  target += std_normal_lpdf(beta_raw)
          + std_normal_lpdf(beta_region)
          + std_normal_lpdf(beta_int_raw)
          + std_normal_lpdf(epsilon_epoch)
          + std_normal_lpdf(epsilon_site)
          + student_t_lpdf(to_vector(tau_site) | 3, 0, 1) - M_region*student_t_lccdf(0 | 3, 0, 1)
          + student_t_lpdf(tau_epoch | 3, 0, 1) - student_t_lccdf(0 | 3, 0, 1);
  // likehoood
  target += bernoulli_logit_lpmf(y | eta);
}
