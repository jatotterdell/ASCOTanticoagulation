// Cumulative logistic model in which
// the first logit has a distinct parameter
// for each covariate, but all other logits
// have shared parameters.

functions {
  // make_cutpoints
  // - p: outcome level probabilities
  vector make_cutpoints(vector p, real scale) {
    int C = rows(p) - 1;
    vector[C] cutpoints;
    real running_sum = 0;
    for(c in 1:C) {
      running_sum += p[c];
      cutpoints[c] = logit(running_sum);
    }
    return scale * cutpoints;
  }

  // Pr(y == k) for k=1,...,K
  // - c: cut-points for outcome levels
  // - eta: linear predictor for each pattern
  matrix Pr(vector c, vector eta, vector eta1) {
    int N = num_elements(eta);
    int K = num_elements(c) + 1;
    matrix[N, K] out;
    // for stability, work on log-scale
    for (n in 1:N) {
      out[n, 1] = log1m_exp(-log1p_exp(-eta1[n] + c[1]));
      out[n, 2] = log_diff_exp(-log1p_exp(-eta1[n] + c[1]), -log1p_exp(-eta[n] + c[2]));
      out[n, K] = -log1p_exp(-eta[n] + c[K-1]);
      for (k in 3:(K - 1)) {
        out[n, k] = log_diff_exp(-log1p_exp(-eta[n] + c[k-1]), -log1p_exp(-eta[n] + c[k]));
      }
    }
    return exp(out);
  }

  // log-likelihood (multinomial)
  // - p: level probabilities for each pattern
  // - y: observed count for each level for each pattern
  vector log_lik(matrix p, matrix y) {
    int N = rows(y);
    int K = cols(y);
    vector[N] out;
    for(n in 1:N) {
      out[n] = 0.0;
      for(k in 1:K) {
        out[n] += y[n, k] * log(p[n, k]);
      }
    }
    return out;
  }
}

data {
  int N; // number of records
  int J; // number of response levels
  int K; // number of covariates
  matrix[N, J] y; // response record x level
  matrix[N, K] X; // design matrix
  vector[J] p_par; // prior for Dirichlet on cuts
  vector[K] beta_sd; // prior SD for beta coefficients
}

parameters {
  simplex[J] p0;      // outcome level probabilities for reference level
  vector[K] beta_raw; // covariate coefficients
  vector[K] beta1_raw; // covariate coefficients
}

transformed parameters {
  vector[J-1] alpha;    // outcome level cuts for reference pattern
  matrix[N, J] p;   // matrix of level probabilities for covariate pattern and outcome level
  vector[N] loglik; // store covariate pattern loglik contribution
  vector[K] beta = beta_sd .* beta_raw;
  vector[K] beta1 = beta_sd .* beta1_raw;
  vector[N] eta = X * beta;
  vector[N] eta1 = X * beta1;
  alpha = make_cutpoints(p0, 1);
  p = Pr(alpha, eta, eta1);
  loglik = log_lik(p, y);
}

model {
  // Prior model
  target += normal_lpdf(beta_raw | 0, 1) +
            normal_lpdf(beta1_raw | 0, 1);
  target += dirichlet_lpdf(p0 | p_par);
  // Observational model
  target += loglik;
}
