functions {
#include "ode.stan";
}

data {
  int<lower = 0> N;           // number of measurement times
  int<lower = 0> N_cv;        // number of future cross-validation observations
  real ts[N + N_cv];                 // measurement times > 0
  real y_init[2];             // initial measured populations
  real<lower = 0> y[N + N_cv, 2];    // measured populations
}

transformed data {
  int<lower = 0> N_tot = N + N_cv;
}

parameters {
  real<lower = 0> theta[4];   // { alpha, beta, gamma, delta }
  real<lower = 0> z_init[2];  // initial population
  real<lower = 0> sigma[2];   // measurement errors
}

transformed parameters {
  real z[N_tot, 2];
  z = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-5, 1e-3, 5e2);
}

model {
  theta[{1, 3}] ~ normal(1, 0.5);
  theta[{2, 4}] ~ normal(0.05, 0.05);
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(10), 1);
  for (k in 1:2) {
    y_init[k] ~ lognormal(log(z_init[k]), sigma[k]);
    y[1:N , k] ~ lognormal(log(z[1:N, k]), sigma[k]);
  }
}

generated quantities {
  real log_likelihood[N_cv, 2];
 
  for (k in 1:2) {
    for (n in 1:N_cv) {
      log_likelihood[n, k] = lognormal_lpdf(y[N + n, k] | log(z[n + N, k]), sigma[k]);
    }
  }
  
}
