data {
  int N;          // number of participants
  int K;          // number of design parameters
  int J;          // number of outcome levels
  int M_region;   // number of regions
  array[N] int y; // outcome level
  matrix[N, K] X; // design matrix, including treatment design, covariates, and eligibilities
  vector[K] beta_sd; // prior for design coefficient parameters
  vector<lower=0> [J] p_par; // dirichlet prior hyper-parameters
  array[N] int<lower=1> region; // region indicator for individual

}

parameters {
  simplex[J] p;
  vector[K] beta_raw;                 // design coefficients
  vector[M_region-1] beta_region_raw; // region effect, first term = 0.
}

transformed parameters {
  ordered[J-1] alpha = logit(cumulative_sum(p[1:(J-1)]));
  vector[K] beta = beta_sd .* beta_raw;
  vector[M_region] beta_region;
  beta_region[1] = 0.0;
  beta_region[2:M_region] = beta_region_raw;
}

model {
  // lp
  vector[N] eta = X * beta + beta_region[region];
  // prior
  // - beta ~ Normal(0, beta_sd)
  // - p(eta = 0) ~ Dirichlet(p_par)
  target += dirichlet_lpdf(p | p_par)
          + std_normal_lpdf(beta_raw)
          + std_normal_lpdf(beta_region);
  // likehoood
  target += ordered_logistic_lpmf(y | eta, alpha);
}

generated quantities {
  array[N] int y_ppc;
  vector[N] eta =
    X * beta + beta_region[region];
  for(n in 1:N) {
    y_ppc[n] = ordered_logistic_rng(eta[n], alpha);
  }
}
