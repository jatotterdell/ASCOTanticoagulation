data {
  int N;          // number of participants
  int K;          // number of design parameters
  int J;          // number of outcome levels
  array[N] int y; // outcome level
  matrix[N, K] X; // design matrix, including treatment design, covariates, and eligibilities
  vector[K] beta_sd; // prior for design coefficient parameters
  vector<lower=0> [J] p_par; // dirichlet prior hyper-parameters
}

parameters {
  simplex[J] p;
  vector[K] beta_raw;                 // design coefficients
}

transformed parameters {
  ordered[J-1] alpha = logit(cumulative_sum(p[1:(J-1)]));
  vector[K] beta = beta_sd .* beta_raw;
}

model {
  // lp
  vector[N] eta =
    X * beta;
  // prior
  // - beta ~ Normal(0, beta_sd)
  // - p(eta = 0) ~ Dirichlet(p_par)
  target += dirichlet_lpdf(p | p_par)
          + std_normal_lpdf(beta_raw)
    ;
  // likehoood
  target += ordered_logistic_lpmf(y | eta, alpha);
}

generated quantities {
  array[N] int y_ppc;
  vector[N] eta = X * beta;
  for(n in 1:N) {
    y_ppc[n] = ordered_logistic_rng(eta[n], alpha);
  }
}
