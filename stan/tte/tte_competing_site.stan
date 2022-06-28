// time to event model - competing risk
// prototype discrete cause-specific time-to-event multinomial logit model
// to account for death in evaluating time to recovery

functions {
  vector mlogit (vector x) {
    return exp(x - log1p(sum(exp(x))));
  }
  row_vector mlogit (row_vector x) {
    return exp(x - log1p(sum(exp(x))));
  }
}

data {
  int N; // number of observations
  int R; // number of event types
  int K; // number of design parameters
  int T; // number of time-points
  int M_site;
  int M_region;
  array[N, R + 1] int<lower=0,upper=1> y; // multinomial outcome
  matrix[N, K] X;
  array[N] int time;
  array[M_site] int<lower=1> region_by_site; // region indicator for each site
  array[N] int<lower=1> site;                // site indicator
  vector[K] beta_sd; // prior for design coefficient parameters
}

parameters {
  matrix[T, R] alpha;
  matrix[K, R] beta_raw;
  matrix[M_site, R] epsilon_site;        // site coefficients
  matrix<lower=0>[M_region, R] tau_site; // region-specific site sd
}

transformed parameters {
  matrix[K, R] beta;
  matrix[M_site, R] gamma_site;
  for(r in 1:R) {
    beta[, r] = beta_sd .* beta_raw[, r];
    // use country specific site variation
    for(m in 1:M_site) {
      gamma_site[m, r] = tau_site[region_by_site[m], r] * epsilon_site[m, r];
    }
  }
}

model {
  matrix[N, R] eta;
  matrix[N, R + 1] lambda;
  for (n in 1:N) {
    for (r in 1:R) {
      eta[n, r] = alpha[time[n], r] + X[n] * beta[,r] +  gamma_site[site[n], r];
    }
    lambda[n, 2:(R+1)] = mlogit(eta[n]);
    lambda[n, 1] = 1 - sum(lambda[n, 2:(R+1)]);
    y[n] ~ multinomial(to_vector(lambda[n]));
  }
  to_vector(alpha) ~ normal(0, 10);
  to_vector(beta) ~ normal(0, 1);
  to_vector(epsilon_site) ~ normal(0, 1);
  to_vector(tau_site) ~ student_t(3, 0, 1);
}
