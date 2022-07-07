// time to event model - competing risk
// discrete cause-specific time-to-event multinomial logit model
// to account for death in evaluating time to recovery

functions {
  // multinomial logit transformation on vector x
  vector mlogit (vector x) {
    return exp(x - log1p(sum(exp(x))));
  }
  // multinomial logit transformation on row_vector x
  row_vector mlogit (row_vector x) {
    return exp(x - log1p(sum(exp(x))));
  }
}

data {
  int N; // number of observations
  int R; // number of event types
  int K; // number of design parameters
  int T; // number of time-points
  int M_site; // number of sites
  int M_region; // number of regions
  int M_epoch; // number of epochs
  array[N, R + 1] int<lower=0,upper=1> y; // multinomial outcome
  matrix[N, K] X; // Design
  array[N] int time; // Study day
  array[M_site] int<lower=1> region_by_site; // region indicator for each site
  array[N] int<lower=1> site;                // site indicator
  array[N] int epoch; // epoch indicator
  vector[K] beta_sd; // prior for design coefficient parameters
}

parameters {
  matrix[T, R] alpha_raw;
  matrix[K, R] beta_raw;
  vector<lower=0>[R] tau_alpha;
  matrix[M_epoch-1, R] epsilon_epoch;
  vector<lower=0>[R] tau_epoch;
  matrix[M_site, R] epsilon_site;        // site coefficients
  matrix<lower=0>[M_region, R] tau_site; // region-specific site sd
}

transformed parameters {
  matrix[T, R] alpha;
  matrix[K, R] beta;
  matrix[M_epoch, R] gamma_epoch;
  matrix[M_site, R] gamma_site;
  for(r in 1:R) {
    beta[, r] = beta_sd .* beta_raw[, r];
    alpha[1,r] = 10*alpha_raw[1,r];
    for(i in 2:T) {
      alpha[i,r] = alpha[i-1,r] + tau_alpha[r] * alpha_raw[i, r];
    }
    for(m in 1:M_site) {
      gamma_site[m, r] = tau_site[region_by_site[m], r] * epsilon_site[m, r];
    }
    gamma_epoch[1, r] = 0.0;
    gamma_epoch[2:M_epoch, r] = tau_epoch[r] * cumulative_sum(epsilon_epoch[, r]);
  }
}

model {
  matrix[N, R] eta;
  matrix[N, R + 1] lambda;
  for (n in 1:N) {
    for (r in 1:R) {
      eta[n, r] = alpha[time[n], r] + X[n] * beta[,r] +  gamma_site[site[n], r] + gamma_epoch[epoch[n], r];
    }
    lambda[n, 2:(R+1)] = mlogit(eta[n]);
    lambda[n, 1] = 1 - sum(lambda[n, 2:(R+1)]);
    y[n] ~ multinomial(to_vector(lambda[n]));
  }
  to_vector(alpha_raw) ~ normal(0, 1);
  to_vector(beta_raw) ~ normal(0, 1);
  to_vector(epsilon_site) ~ normal(0, 1);
  to_vector(tau_site) ~ student_t(3, 0, 1);
  to_vector(epsilon_epoch) ~ normal(0, 1);
  tau_epoch ~ student_t(3, 0, 1);
}
