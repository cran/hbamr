data {
  int<lower = 1> N;                       // n of individuals
  int<lower = 1> J;                       // n of items
  int<lower = 1> G;                       // n of groups
  int<lower = 1> N_obs;                   // n of observations
  array[N_obs] int<lower = 1> ii;         // index i in matrix
  array[N_obs] int<lower = 1> jj;         // index j in matrix
  array[N] int<lower = 1> gg;             // group-index of i
  int<lower = 1> B;                       // length of scale -1 / 2
  int<lower = 1, upper = J> L;            // left pole
  int<lower = 1, upper = J> R;            // right pole
  vector<lower = -B, upper = B>[N_obs] Y; // reported stimuli positions
  vector<lower = -B, upper = B>[N] V;     // reported self-placements
  int<lower = 0, upper = 1> CV;           // indicator of cross-validation
  vector<lower = 0, upper = 1>[N_obs] holdout; // holdout for cross-validation
  real<lower = 0> sigma_alpha;            // sd of prior on alpha
  real<lower = 0> sigma_beta;             // sd of prior on log(beta)
  real<lower = 0> sigma_mu_alpha;         // sd of prior on mu_alpha
  real<lower = 0> sigma_mu_beta;          // sd of prior on mu_beta
  int<lower = 0, upper = 1> MCMC;         // indicator of fitting method
}

transformed data {
  real nu = 6;                            // concentration of etas
  real tau = B / 4.0;                     // scale of prior on errors
  real eta_scale = tau * J;
  real psi = 6;                           // implies 13% prior prob. of flipping
  vector<lower = 0, upper = 1>[N_obs] not_holdout = 1 - holdout;
  real mean_mu_simplexes = 1.0 / G;       // for later scaling of simplexes
  real sd_mu_simplexes = sqrt(mean_mu_simplexes * (1 - mean_mu_simplexes) / (50 * G + 1));
}

parameters {
  matrix[N, 2] alpha_raw;                 // shift parameter, split, raw
  matrix[N, 2] beta_raw;                  // stretch parameter, split, raw
  ordered[2] theta_lr;                    // left and right pole
  array[J] real theta_raw;                // remaining stimuli
  simplex[G] mu_alpha_raw;                // group-level mean of alpha, raw
  simplex[G] mu_beta_raw;                 // group-level mean of log(beta), raw
  vector<lower = 0>[N] eta;               // mean ind. error variance x J^2
  simplex[J] rho;                         // stimuli-shares of variance
  vector[N] lambda_raw;                   // raw mixing proportion, flipping
}

transformed parameters {
  array[J] real theta;                    // latent stimuli position
  matrix[N, 2] alpha0;                    // shift parameter, split
  matrix[N, 2] beta0;                     // stretch parameter, split
  vector[N_obs] log_lik;                  // pointwise log-likelihood for Y
  vector<lower = 0, upper = 1>[N] lambda = inv_logit(psi + lambda_raw * 5); // prob. of non-flipping
  vector[G] mu_alpha = ((mu_alpha_raw - mean_mu_simplexes) / sd_mu_simplexes) * sigma_mu_alpha;
  vector[G] mu_beta = ((mu_beta_raw - mean_mu_simplexes) / sd_mu_simplexes) * sigma_mu_beta;
  theta = theta_raw;
  theta[L] = theta_lr[1];                 // safeguard to ensure identification
  theta[R] = theta_lr[2];

  for (i in 1:N) {
    alpha0[i, 1] = alpha_raw[i, 1] * sigma_alpha + mu_alpha[gg[i]]; // non-centered specifications
    alpha0[i, 2] = alpha_raw[i, 2] * sigma_alpha + mu_alpha[gg[i]];
    beta0[i, 1] = exp(beta_raw[i, 1] * sigma_beta + mu_beta[gg[i]]);
    beta0[i, 2] = -exp(beta_raw[i, 2] * sigma_beta); // group-level mean not added for flipped state
  }

  for (n in 1:N_obs) {
    log_lik[n] = log_mix( lambda[ii[n]],
      normal_lpdf(Y[n] | alpha0[ii[n], 1] + beta0[ii[n], 1] * theta[jj[n]],
        sqrt(eta[ii[n]]) * rho[jj[n]]),
      normal_lpdf(Y[n] | alpha0[ii[n], 2] + beta0[ii[n], 2] * theta[jj[n]],
        sqrt(eta[ii[n]]) * rho[jj[n]]) );
  }
}

model {
  theta_raw ~ normal(0, B / 2.0);
  theta_lr ~ normal(0, B / 2.0);
  alpha_raw[, 1] ~ normal(0, 1);
  alpha_raw[, 2] ~ normal(0, 1);
  beta_raw[, 1] ~ normal(0, 1);
  beta_raw[, 2] ~ normal(0, 1);
  mu_alpha_raw ~ dirichlet(rep_vector(50, G));
  mu_beta_raw ~ dirichlet(rep_vector(50, G));
  eta ~ scaled_inv_chi_square(nu, eta_scale);
  rho ~ dirichlet(rep_vector(50, J));
  lambda_raw ~ normal(0, 1);

  if (CV == 0)
    target += sum(log_lik);
  else
    target += sum(log_lik .* not_holdout);
}

generated quantities {
  vector[N] kappa;
  vector[N] chi;
  real<lower = 0> min_rho = min(rho);
  if (MCMC == 1)
    kappa = to_vector(bernoulli_rng(lambda));
  else
    kappa = to_vector(round(lambda)); // Rounding to MAP instead of sampling
  vector[N] alpha = (kappa .* alpha0[, 1]) + ((1 - kappa) .* alpha0[, 2]);
  vector[N] beta = (kappa .* beta0[, 1]) + ((1 - kappa) .* beta0[, 2]);
  if (MCMC == 1)
    chi = (V - to_vector(normal_rng(0, sqrt(eta) * min_rho)) - alpha) ./ beta;
  else
    chi = (V - alpha) ./ beta;
}
