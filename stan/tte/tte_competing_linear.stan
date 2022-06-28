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
  array[N, R + 1] int<lower=0,upper=1> y; // multinomial outcome
  matrix[N, K] X;
  array[N] int time;
  vector[K] beta_sd; // prior for design coefficient parameters
}

parameters {
  vector[R] alpha;
  matrix[K, R] beta_raw;
}

transformed parameters {
  matrix[K, R] beta;
  for(r in 1:R) {
    beta[, r] = beta_sd .* beta_raw[, r];
  }
}

model {
  matrix[N, R] eta;
  matrix[N, R + 1] lambda;
  for (n in 1:N) {
    for (r in 1:R) {
      eta[n, r] = alpha[r]*time[n] + X[n] * beta[,r];
    }
    lambda[n, 2:(R+1)] = mlogit(eta[n]);
    lambda[n, 1] = 1 - sum(lambda[n, 2:(R+1)]);
    y[n] ~ multinomial(to_vector(lambda[n]));
  }
  alpha ~ normal(0, 10);
  to_vector(beta_raw) ~ normal(0, 1);
}
