data {
  int<lower=0> N; // Number of observations (lineage pairs)
  int<lower=0> K; // Number of predictors, including the intercept
  vector<lower=0,upper=1>[N] Y; // Response variable, bounded between 0 and 1
  matrix[N, K] X; // Design matrix with intercept column and observations for each predictor variable in subsequent columns
  matrix[N, N] Cp; // The lineage-pair covariance matrix
  real coef_mean; // mean for the normal prior on the coefficients
  real<lower=0> coef_sd; // sd for the normal prior on coefficients
  real<lower=0> phi_shape; // shape parameter for the gamma prior on phi
  real<lower=0> phi_rate; // rate parameter for the gamma prior on phi
  real sig2_mean; // mean of the lognormal prior on sig2_scale
  real sig2_sd; // sd of the lognormal prior on sig2_scale
  int<lower=1, upper=4> link_choice; // 1 = logit, 2 = probit, 3 = cloglog, 4 = loglog
  int<lower=1, upper=2> model_type; // 1 = beta.ols, 2 = beta.mm
}

parameters {
  vector[K] Coef; // coefficient vector to be estimated
  real<lower=0> phi; // precision parameter for beta distribution
  real<lower=0> sig2_scale[model_type == 2 ? 1 : 0]; // scaling parameter for variance in species-specific effects
  vector[model_type == 2 ? N : 0] pair_effects; // random effects (intercepts) for species pairs
}

transformed parameters {
  vector[N] eta;
  
  if (model_type == 2) {
    eta = X * Coef + pair_effects;
  } else {
    eta = X * Coef;
  }
}

model {
  // Priors
  Coef ~ normal(coef_mean, coef_sd);
  phi ~ gamma(phi_shape, phi_rate);

  if (model_type == 2) {
    sig2_scale[1] ~ lognormal(sig2_mean, sig2_sd);
    pair_effects ~ multi_normal(rep_vector(0, N), sig2_scale[1] * Cp);
  }

  // Link function
  vector[N] mu;
  if (link_choice == 1) {
    mu = inv_logit(eta); // logit link
  } else if (link_choice == 2) {
    for (n in 1:N) {
      mu[n] = normal_cdf(eta[n], 0, 1); // probit link
    }
  } else if (link_choice == 3) {
    mu = 1 - exp(-exp(eta)); // cloglog link
  } else if (link_choice == 4) {
    mu = exp(-exp(-eta)); // loglog link
  }

  // Likelihood for response
  Y ~ beta(mu * phi, (1 - mu) * phi);
}

generated quantities {
  vector[N] loglik;
  vector[N] mu; // Mean of the Beta distribution

  // Calculate the predicted values based on the link function
  if (link_choice == 1) {
    mu = inv_logit(eta); // logit link
  } else if (link_choice == 2) {
    for (n in 1:N) {
      mu[n] = normal_cdf(eta[n], 0, 1); // probit link
    }
  } else if (link_choice == 3) {
    mu = 1 - exp(-exp(eta)); // cloglog link
  } else if (link_choice == 4) {
    mu = exp(-exp(-eta)); // loglog link
  }

  // Calculate log likelihood in a vectorized manner
  for (n in 1:N) {
    loglik[n] = beta_lpdf(Y[n] | mu[n] * phi, (1 - mu[n]) * phi);
  }
}
