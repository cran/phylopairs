data {
  int<lower=0> N; // Number of observations (lineage pairs)
  int<lower=0> K; // Number of predictors, including the intercept
  vector[N] Y; // Response variable
  matrix[N, K] X; // Design matrix with intercept column and observations for each predictor variable in subsequent columns
  matrix[N, N] Cp; // The lineage-pair covariance matrix (used only in models 2 and 3)
  real coef_mean; // mean for the normal prior on the coefficients
  real<lower=0> coef_sd; // sd for the normal prior on coefficients
  real sig2_mean; // mean of the lognormal prior on sig2_scale (used only in models 2 and 3)
  real sig2_sd; // sd of the lognormal prior on sig2_scale (used only in models 2 and 3)
  real sig2_loc; // location parameter of cauchy prior on sigma_resid (used only in models 1 and 3)
  real sig2_sc; // scale parameter of cauchy prior on sigma_resid (used only in models 1 and 3)
  int<lower=1, upper=3> model_type; // 1 for OLS, 2 for phylogenetic GLS, 3 for mixed model
}

parameters {
  vector[K] Coef; // coefficient vector to be estimated
  real<lower=0> sigma_resid[model_type == 2 ? 0 : 1]; // standard deviation of residual errors (used only in models 1 and 3)
  real<lower=0> sig2_scale[model_type == 2 || model_type == 3 ? 1 : 0]; // scaling parameter for phylo effects (used in models 2 and 3)
  vector[model_type == 3 ? N : 0] pair_effects; // random effects (intercepts) for species pairs (used only in model 3)
}

transformed parameters {
  vector[N] mu;
  
  if (model_type == 3) {
    mu = X * Coef + pair_effects;
  } else {
    mu = X * Coef;
  }
}

model {
  // Priors
  Coef ~ normal(coef_mean, coef_sd);

  if (model_type == 1) {
    sigma_resid[1] ~ cauchy(sig2_loc, sig2_sc);

    // Likelihood for response
    Y ~ normal(mu, sigma_resid[1]);

  } else if (model_type == 2) {
    sig2_scale[1] ~ lognormal(sig2_mean, sig2_sd);

    // Likelihood for response
    Y ~ multi_normal_cholesky(mu, sqrt(sig2_scale[1]) * cholesky_decompose(Cp));

  } else if (model_type == 3) {
    sigma_resid[1] ~ cauchy(sig2_loc, sig2_sc);
    sig2_scale[1] ~ lognormal(sig2_mean, sig2_sd);
    pair_effects ~ multi_normal(rep_vector(0, N), sig2_scale[1] * Cp);

    // Likelihood for response
    Y ~ normal(mu, sigma_resid[1]);
  }
}

generated quantities {
  vector[N] y_pred; // Predicted values of Y
  vector[N] loglik; // Log likelihood of each observation
  
  if (model_type == 1) {
    for (n in 1:N) {
      y_pred[n] = normal_rng(mu[n], sigma_resid[1]);
      loglik[n] = normal_lpdf(Y[n] | mu[n], sigma_resid[1]);
    }
  } else if (model_type == 2) {
    for (n in 1:N) {
      y_pred[n] = normal_rng(mu[n], sqrt(sig2_scale[1] * Cp[n, n]));
      loglik[n] = normal_lpdf(Y[n] | mu[n], sqrt(sig2_scale[1] * Cp[n, n]));
    }
  } else if (model_type == 3) {
    for (n in 1:N) {
      y_pred[n] = normal_rng(mu[n], sigma_resid[1]);
      loglik[n] = normal_lpdf(Y[n] | mu[n], sigma_resid[1]);
    }
  }
}
