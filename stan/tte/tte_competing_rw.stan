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
  matrix[T, R] alpha_raw;
  matrix[K, R] beta_raw;
  vector<lower=0>[R] tau_alpha;
}

transformed parameters {
  matrix[T, R] alpha;
  matrix[K, R] beta;
  // define random-walk(1) prior
  for(r in 1:R) {
    beta[, r] = beta_sd .* beta_raw[, r];
    alpha[1,r] = 10*alpha_raw[1,r];
    for(i in 2:T) {
      alpha[i,r] = alpha[i-1,r] + tau_alpha[r] * alpha_raw[i, r];
    }
  }
}

model {
  matrix[N, R] eta;
  matrix[N, R + 1] lambda;
  for (n in 1:N) {
    for (r in 1:R) {
      eta[n, r] = alpha[time[n], r] + X[n] * beta[,r];
    }
    lambda[n, 2:(R+1)] = mlogit(eta[n]); // Transform to hazard
    lambda[n, 1] = 1 - sum(lambda[n, 2:(R+1)]); // Unrecovered
    y[n] ~ multinomial(to_vector(lambda[n]));
  }
  to_vector(alpha_raw) ~ normal(0, 1); // alpha[r] | tau[r] ~ RW(1, tau[r])
  to_vector(beta_raw) ~ normal(0, 1);  // beta[k,r] ~ Normal(0, 1)
  tau_alpha ~ student_t(3, 0, 1);      // tau[r] ~ t(3, 0, 1)
}
