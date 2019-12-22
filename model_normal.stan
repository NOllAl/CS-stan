functions {
#include "ode.stan";
}

data {
  int<lower = 0> N;           // number of measurement times
  real ts[N];                 // measurement times > 0
  real y_init[2];             // initial measured populations
  real<lower = 0> y[N, 2];    // measured populations
}

parameters {
  real<lower = 0> theta[4];   // { alpha, beta, gamma, delta }
  real<lower = 0> z_init[2];  // initial population
  real<lower = 0> sigma[2];   // measurement errors
}

transformed parameters {
  real z[N, 2]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-5, 1e-3, 5e2);
}

model {
  theta[{1, 3}] ~ normal(1, 0.5);
  theta[{2, 4}] ~ normal(0.05, 0.05);
  sigma ~ lognormal(-1, 1);
  z_init ~ normal(log(10), 1);
  for (k in 1:2) {
    y_init[k] ~ normal(z_init[k], sigma[k]);
    y[ , k] ~ normal(z[, k], sigma[k]);
  }
}

generated quantities {
  real y_init_rep[2];
  real y_rep[N, 2];
  for (k in 1:2) {
    y_init_rep[k] = normal_rng(z_init[k], sigma[k]);
    for (n in 1:N)
      y_rep[n, k] = normal_rng(z[n, k], sigma[k]);
  }
}
