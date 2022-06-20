// Simple logistic regression model for Bernoulli outcome

data
{
  int N;             // number of participants
  int K;             // number of design parameters
  int M_region;      // number of regions
  int M_site;        // total number of sites across all regions
  int y[N];          // outcome
  matrix[N, K] X;    // design matrix, including intercept, treatment design, covariates, and eligibilities, excluding region and site
  vector[K] beta_sd; // prior for design coefficient parameters

  int<lower=1> region[N];              // region indicator for individual
  int<lower=1> region_by_site[M_site]; // region indicator for each site
  int<lower=1> site[N];                // site indicator
}

parameters
{
  vector[K] beta_raw;
}

transformed parameters
{
  vector[K] beta = beta_sd .* beta_raw;
}

model
{
  // lp
  vector[N] eta = X * beta;
  // prior
  target += std_normal_lpdf(beta_raw);
  // likelihood
  target += bernoulli_logit_lpmf(y | eta);
}
