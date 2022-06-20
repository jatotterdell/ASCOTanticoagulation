//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ normal(mu, sigma);
}

data {
  int N;          // number of participants
  int K;          // number of design parameters
  int J;          // number of outcome levels
  int M_region;   // number of regions
  int M_site;     // total number of sites across all regions
  int M_epoch;    // number of epochs
  array[N] int y; // outcome level
  matrix[N, K] X; // design matrix, including treatment design, covariates, and eligibilities
  vector[K] beta_sd; // prior for design coefficient parameters
  vector<lower=0> [J] p_par; // dirichlet prior hyper-parameters

  array[N] int<lower=1> region;              // region indicator for individual
  array[N] int<lower=1> site;                // site indicator
  array[N] int<lower=1> epoch;               // epoch indicator
}

parameters {
  simplex[J] p;
  vector[K] beta_raw;                 // design coefficients
  vector[M_region-1] beta_region_raw; // region effect, first term = 0.
  vector[M_site] xi_site_raw;        // site coefficients
  vector[M_epoch-1] epsilon_epoch;    // epoch coefficients, first term = 0.
  vector<lower=0>[M_region] sigma_site; // region-specific site sd
  real<lower=0> sigma_epoch;            // cohort sd
}

transformed parameters {
  ordered[J-1] alpha = logit(cumulative_sum(p[1:(J-1)]));
  vector[K] beta = beta_sd .* beta_raw;
  vector[M_region] beta_region;
  vector[M_site]       xi_site;
  vector[M_epoch]  tau_epoch;

  beta_region[1] = 0.0;
  beta_region[2:M_region] = beta_region_raw;

  for(n in 1:N) {
    xi_site[site[n]] = xi_site_raw[site[n]] * sigma_site[region[n]];
  }

  // define random-walk(1) prior
  tau_epoch[1] = 0.0;
  tau_epoch[2:M_epoch] = sigma_epoch * cumulative_sum(epsilon_epoch);
}

model {
  // lp
  vector[N] eta =
    X * beta + beta_region[region] + xi_site[site] + tau_epoch[epoch];
  // prior
  // - beta ~ Normal(0, beta_sd)
  // - p(eta = 0) ~ Dirichlet(p_par)
  target += dirichlet_lpdf(p | p_par)
          + std_normal_lpdf(beta_raw)
          + std_normal_lpdf(beta_region)
          + std_normal_lpdf(xi_site_raw)
          + std_normal_lpdf(epsilon_epoch)
          + student_t_lpdf(to_vector(sigma_site) | 3, 0, 1) - M_region*student_t_lccdf(0 | 3, 0, 1)
          + student_t_lpdf(sigma_epoch | 3, 0, 1) - student_t_lccdf(0 | 3, 0, 1);
  // likehoood
  target += ordered_logistic_lpmf(y | eta, alpha);
}

generated quantities {
  array[N] int y_ppc;
  vector[N] eta =
    X * beta + beta_region[region] + xi_site[site] + tau_epoch[epoch];
  for(n in 1:N) {
    y_ppc[n] = ordered_logistic_rng(eta[n], alpha);
  }
}