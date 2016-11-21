data {
  int N;                               # sample size
  int K;                               # number of items
  matrix[N, K] X;                      # items
  vector[N] y;                         # target variable
}

parameters {
  vector<lower=-10, upper=10>[N] f;    # latent variable
  real b0;                             # intercept in the structural model
  real beta;                           # coefficient for f in the structural model

  vector<lower=0, upper=1>[K] lambda;  # loadings
  vector[K] fints;                     # intercepts for the factor model


  vector<lower=0>[K] u;                # uniquenesses
  real<lower=0> sigma;                 # error scale for target
  // real<lower=0> tau;                # factor sd; fix to 1 for identification
}

model {
  vector[N] mu;

  // f ~ normal(0, tau);
  f ~ normal(0, 1);
  lambda ~ normal(0, 1);

  u ~ student_t(3, 0 , 5);
  sigma ~ student_t(3, 0 , 5);
  // tau ~ student_t(3, 0 , 5);

  b0 ~ normal(0, 10);
  beta ~ normal(0, 10);

  for (k in 1:K){
    X[,k] ~ normal(fints[k] + lambda[k]*f, u[k]);
  }

  mu = b0 + beta*f;
  y ~ normal(mu, sigma);
}