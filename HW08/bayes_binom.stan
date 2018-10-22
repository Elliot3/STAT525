data {
  // Data:
  int<lower=1> n; // number of observations
  int<lower=0> y; // number of "yes" answers
  
  // Input hyperparameters as data:
  real<lower=0> alpha;
  real<lower=0> beta;
}

parameters {
  
  // Probabilty of yes:
  real theta;
  
}

model {
  
  // Likelihood:
  y ~ binomial(n, theta);
  
  // Priors: 
  theta ~ beta(alpha, beta);

}

