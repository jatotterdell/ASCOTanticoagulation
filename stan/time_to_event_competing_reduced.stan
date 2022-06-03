// time to event model - competing risk
// discrete cause-specific time-to-event multinomial logit model
// to account for death in modelling time to recovery

functions {

  // exp(x) / (1 + sum(exp(x)))
  vector mlogit (vector x) {
    return exp(x - log1p(sum(exp(x))));
  }

  // exp(x) / (1 + sum(exp(x)))
  row_vector mlogit (row_vector x) {
    return exp(x - log1p(sum(exp(x))));
  }

}

data {
  int N; // number of observations
  int R; // number of event types
  int K; // number of design parameters
  array[N, R + 1] int<lower=0,upper=1> y; // multinomial outcome
  matrix[N, K] X; // design matrix (assumed same across all outcomes)
}

parameters {
  matrix[K, R] beta; // beta coefficents for K terms in X by each of R events
}

model {
  matrix[N, R] eta;
  matrix[N, R + 1] lambda;
  for (r in 1:R) {
    eta[, r] = X * beta[, r]; // linear-predictor
  }
  for (n in 1:N) {
      lambda[n, 2:(R+1)] = mlogit(eta[n]); // cause-specific hazard
      lambda[n, 1] = 1 - sum(lambda[n, 2:(R+1)]); // conditional survival (no recovery or death)
      y[n] ~ multinomial(to_vector(lambda[n])); // multi-logit likelihood
  }
  to_vector(beta) ~ normal(0, 10);
}
