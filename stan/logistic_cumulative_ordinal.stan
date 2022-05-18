// cumulative logistic model for ordinal outcome

data {
  int N; // number of observations
  int K; // number of design parameters
  int J; // number of outcome levels
  matrix[N, K] X; // design matrix (no intercept)
  array[N] int y; // outcome level
  vector[K] beta_sd; // prior standard deviation for beta parameters
  vector<lower=0> [J] p_par; // dirichlet prior hyper-parameters
}

parameters {
  vector[K] beta_raw;
  simplex[J] p;
}

transformed parameters {
  vector[K] beta = beta_sd .* beta_raw;
  ordered[J-1] alpha = logit(cumulative_sum(p[1:(J-1)]));
}

model {
  // lp
  vector[N] eta = X * beta;
  // prior
  // - beta ~ Normal(0, beta_sd)
  // - p(eta = 0) ~ Dirichlet(p_par)
  target += std_normal_lpdf(beta_raw) +
            dirichlet_lpdf(p | p_par);
  // likelihood
  target += ordered_logistic_lpmf(y | eta, alpha);
}
