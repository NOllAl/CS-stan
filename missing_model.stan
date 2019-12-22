functions {
#include "ode.stan";
}

data {
  int<lower = 0> J;           // number of observations
  int<lower = 0> species[J];  // species associated with observation j
  int tt[J];                  // time associated with observation j
  real<lower = 0> yy[J];      // Acutal observations
  
  int<lower = 0> J_mis;       // number of missing observations
  int<lower = 0> species_mis[J_mis]; // species associated with missing obs j
  int tt_mis[J_mis];          // time associated with missing observation
  
  int<lower = 0> N;           // number of timesteps (could also be calculated as max(tt);
  real ts[N];                 // measurement times > 0
  real y_init[2];             // initial measured populations
}

parameters {
  real<lower = 0> theta[4];   // { alpha, beta, gamma, delta }
  real<lower = 0> z_init[2];  // initial population
  real<lower = 0> sigma[2];   // measurement errors
}

transformed parameters {
  real z[N, 2];
  
  z = integrate_ode_rk45(
    dz_dt, z_init, 0, ts, theta,
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
  }
  
  for (j in 1:J) {
    yy[j] ~ lognormal(log(z[tt[j], species[j]]), sigma[species[j]]);
  }
}

generated quantities {
  real y_mis[J_mis];
  
  for (j in 1:J_mis) {
    y_mis[J_mis] = lognormal_rng(log(z[tt_mis[j], species_mis[j]]), sigma[species_mis[j]]);
  }
}