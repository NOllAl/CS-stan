real[] dz_dt(real t,       // time
             real[] z,     // system state {prey, predator}
             real[] theta, // parameters
             real[] x_r,   // unused data
             int[] x_i){
    real u = z[1];
    real v = z[2];

    real alpha = theta[1];
    real beta = theta[2];
    real gamma = theta[3];
    real delta = theta[4];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;
    return { du_dt, dv_dt };
}