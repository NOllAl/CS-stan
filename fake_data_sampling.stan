functions {
#include "ode.stan";
}

data {
  int<lower=0> N;           // number of measurement times
  real ts[N];               // measurement times > 0
  real z_init[2];           // initial measured populations
  real<lower=0> theta[4];   // { alpha, beta, gamma, delta}
  real<lower=0> sigma[2];   // measurement errors
}

transformed parameters {
  real z[N, 2]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-5, 1e-3, 5e2);
}

generated quantities {
  real<lower=0> y_init[2];
  real<lower=0> y[N, 2];
  for (k in 1:2) {
    y_init[k] = lognormal_rng(log(z_init[k]), sigma[k]);
    y[ , k] = lognormal_rng(log(z[ , k]), sigma[k]);
  }
}