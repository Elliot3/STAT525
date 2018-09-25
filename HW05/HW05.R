## No scientific notation

options(scipen = 999)

## Load in packages

library(fGarch)



########## Problem 1 ##########



## Input the known parameters and data

y <- c(15,9,15,8,8,20,25,31,13,12,11,15,11,13,13,11,11,12,6,14)
n <- length(y)
S <- 10000
alpha <- 5
beta <- alpha * 15



##### Part i #####



## Sample from the Gamma

samp_1i <- rgamma(n = S, shape = n + alpha, rate = beta + sum(y))



##### Part ii #####



## Construct the replicate data sets

y_rep <- matrix(0, nrow = S, ncol = n)

for(s in 1:S) {
    
    y_rep[s,] <- rexp(n = n, rate = samp_1i[s])
    
}



##### Part iii #####



## Compute the required test statistics for each replicate data set

y_rep_min <- apply(y_rep, 1, min)
y_rep_q25 <- apply(y_rep, 1, quantile, probs = 0.25)
y_rep_med <- apply(y_rep, 1, median)
y_rep_sd <- apply(y_rep, 1, sd)

par(mfrow = c(2, 2))

hist(x = y_rep_min)
abline(v = min(y), col = 'blue')

hist(x = y_rep_q25)
abline(v = quantile(y, 0.25), col = 'blue')

hist(x = y_rep_med)
abline(v = median(y), col = 'blue')

hist(x = y_rep_sd)
abline(v = sd(y), col = 'blue')



##### Part iv #####



## Compute the posterior predictive p-values for each statistic

pval_min <- mean(y_rep_min > min(y))
pval_q25 <- mean(y_rep_q25 > quantile(y, 0.25))
pval_med <- mean(y_rep_med > median(y))
pval_sd <- mean(y_rep_sd > sd(y))

## Probability of waiting at most 10 minutes

wait_10 <- mean(y_rep < 10)

## Output the results

cat("Minimum p-value: ", round(pval_min, 4))
cat("25% Quantile p-value: ", round(pval_q25, 4))
cat("Median p-value: ", round(pval_med, 4))
cat("Standard Deviation p-value: ", round(pval_sd, 4))
cat("Probability 10 Minute Wait: ", round(wait_10, 4))



########## Problem 2 ##########



## Input the known parameters and data

z <- log(y)
n <- length(z)
S <- 10000
alpha <- 0.01
beta <- 0.01

## Draw the sample from the posterior

s_squared <- sd(z)^2
post_tau <- rgamma(n = S, 
                  shape = alpha + (n - 1)/2,
                  rate = beta + (n - 1) * (s_squared / 2))

post_sigma <- 1/sqrt(post_tau)

Q_mu <- n * post_tau
ell_mu <- n * post_tau * mean(z)
post_mu <- rnorm(n = S, mean = (Q_mu^-1) * ell_mu, sd = sqrt(Q_mu^-1))



##### Part i #####



## Construct the replicate data set

z_rep <- matrix(0, nrow = S, ncol = n)

for(s in 1:S) {
    
    z_rep[s,] <- rnorm(n = n, mean = post_mu[s], sd = post_sigma[s])
    
}

## Compute the required test statistics for each replicate data set

z_rep_min <- apply(z_rep, 1, min)
z_rep_q25 <- apply(z_rep, 1, quantile, probs = 0.25)
z_rep_med <- apply(z_rep, 1, median)
z_rep_sd <- apply(z_rep, 1, sd)

par(mfrow = c(2, 2))

hist(x = z_rep_min)
abline(v = min(z), col = 'blue')

hist(x = z_rep_q25)
abline(v = quantile(z, 0.25), col = 'blue')

hist(x = z_rep_med)
abline(v = median(z), col = 'blue')

hist(x = z_rep_sd)
abline(v = sd(z), col = 'blue')

## Compute the posterior predictive p-values for each statistic

pval_min <- mean(z_rep_min > min(z))
pval_q25 <- mean(z_rep_q25 > quantile(z, 0.25))
pval_med <- mean(z_rep_med > median(z))
pval_sd <- mean(z_rep_sd > sd(z))

## Output the results

cat("Minimum p-value: ", round(pval_min, 4))
cat("25% Quantile p-value: ", round(pval_q25, 4))
cat("Median p-value: ", round(pval_med, 4))
cat("Standard Deviation p-value: ", round(pval_sd, 4))



##### Part ii #####



## Probability of waiting at most 10 minutes

wait_10 <- mean(z_rep < log(10))

## Output the results

cat("Probability 10 Minute Wait: ", round(wait_10, 4))



########## Problem 3 ##########



## Load in the data and known parameters

y <- c(1.960,
      1.823,
      0.968,
      2.753,
      0.779,
      3.003,
      2.463,
      2.197,
      2.806,
      3.137,
      -1.813,
      -4.242,
      -19.780,
      8.597,
      -3.134,
      7.486,
      12.621,
      -11.511,
      -12.746,
      17.244)
n <- length(y)
S <- 10000
alpha <- 0.01
beta <- 0.01



##### Part i #####



## Draw the sample from the posterior

s_squared <- sd(y)^2
post_tau <- rgamma(n = S, 
                   shape = alpha + (n - 1)/2,
                   rate = beta + (n - 1) * (s_squared / 2))

post_sigma <- 1/sqrt(post_tau)

Q_mu <- n * post_tau
ell_mu <- n * post_tau * mean(y)
post_mu <- rnorm(n = S, mean = (Q_mu^-1) * ell_mu, sd = sqrt(Q_mu^-1))



##### Part ii #####



## Construct the replicate data set

y_rep <- matrix(0, nrow = S, ncol = n)

for(s in 1:S) {
    
    y_rep[s,] <- rnorm(n = n, mean = post_mu[s], sd = post_sigma[s])
    
}



##### Part iii #####



## Compute the variance difference of the first 10 minus second 10 observation

y_1_var <- var(y[1:10])
y_2_var <- var(y[11:20])

var_diff <- y_1_var - y_2_var

## Compute the same variance for the replicate data

var_diff_vec <- numeric()

for (i in 1:S) {
    
    y_1_temp <- var(y_rep[i, 1:10])
    y_2_temp <- var(y_rep[i, 11:20])
    var_diff_vec[i] <- y_1_temp - y_2_temp
    
}

## Plot the histogram of the results

par(mfrow = c(1, 1))

hist(var_diff_vec, breaks = 50)
abline(v = var_diff, col = 'blue')



##### Part iv #####



## Posterior predictive p-value

post_pred <- mean(var_diff < var_diff_vec)

## Output the results

cat("Probability 10 Minute Wait: ", round(post_pred, 4))



########## Problem 4 ##########



##### Exercise 10.6 #####



### Part a ####



## Input the known parameters

S <- 100
mu <- 0
var <- 1
deg_free <- 3

## Sample theta from g

theta <- rstd(n = S, mean = mu, sd = sqrt(var), nu = deg_free)

## Calculate and plot the importance ratios

import_rat <- dnorm(x = theta, mean = mu, sd = sqrt(var)) /
    dstd(x = theta, mean = mu, sd = sqrt(var), nu = deg_free)

hist(import_rat, breaks = 20)



### Part b ###



## Compute the importance sampled EX and VarX

ex_imp <- ((1/S) * sum(theta * import_rat)) / ((1/S) * sum(import_rat))
var_imp <- ((1/S) * sum(theta^2 * import_rat)) / ((1/S) * sum(import_rat)) - ex_imp^2

## Compute the true EX and VarX

ex_true <- 0
var_true <- 1

## Output the results

cat("Importance Sampled Expected Value: ", round(ex_imp, 4))
cat("True Expected Value: ", round(ex_true, 4))
cat("Importance Sampled Variance: ", round(var_imp, 4))
cat("True Variance: ", round(var_true, 4))



### Part c ###



## Input the known parameters

S <- 10000
mu <- 0
var <- 1
deg_free <- 3

## Sample theta from g

theta <- rstd(n = S, mean = mu, sd = sqrt(var), nu = deg_free)

## Calculate and plot the importance ratios

import_rat <- dnorm(x = theta, mean = mu, sd = sqrt(var)) /
    dstd(x = theta, mean = mu, sd = sqrt(var), nu = deg_free)

hist(import_rat, breaks = 20)

## Compute the importance sampled EX and VarX

ex_imp <- ((1/S) * sum(theta * import_rat)) / ((1/S) * sum(import_rat))
var_imp <- ((1/S) * sum(theta^2 * import_rat)) / ((1/S) * sum(import_rat)) - ex_imp^2

## Compute the true EX and VarX

ex_true <- 0
var_true <- 1

## Output the results

cat("Importance Sampled Expected Value: ", round(ex_imp, 4))
cat("True Expected Value: ", round(ex_true, 4))
cat("Importance Sampled Variance: ", round(var_imp, 4))
cat("True Variance: ", round(var_true, 4))



### Part d ###



## Compute the effective sample size

S_eff <- 1 / (sum((import_rat / sum(sample(import_rat)))^2))

## Output the results

cat("Effective Sample Size: ", round(S_eff, 4))



##### Exercise 10.7 #####



### Part a ####



## Input the known parameters

S <- 100
mu <- 0
var <- 1
deg_free <- 3

## Sample theta from g

theta <- rnorm(n = S, mean = mu, sd = sqrt(var))

## Calculate and plot the importance ratios

import_rat <- dstd(x = theta, mean = mu, sd = sqrt(var), nu = deg_free) /
    dnorm(x = theta, mean = mu, sd = sqrt(var))

hist(import_rat, breaks = 20)



### Part b ###



## Compute the importance sampled EX and VarX

ex_imp <- ((1/S) * sum(theta * import_rat)) / ((1/S) * sum(import_rat))
var_imp <- ((1/S) * sum(theta^2 * import_rat)) / ((1/S) * sum(import_rat)) - ex_imp^2

## Compute the true EX and VarX

ex_true <- 0
var_true <- 1

## Output the results

cat("Importance Sampled Expected Value: ", round(ex_imp, 4))
cat("True Expected Value: ", round(ex_true, 4))
cat("Importance Sampled Variance: ", round(var_imp, 4))
cat("True Variance: ", round(var_true, 4))



### Part c ###



## Input the known parameters

S <- 10000
mu <- 0
var <- 1
deg_free <- 3

## Sample theta from g

theta <- rnorm(n = S, mean = mu, sd = sqrt(var))

## Calculate and plot the importance ratios

import_rat <- dstd(x = theta, mean = mu, sd = sqrt(var), nu = deg_free) /
    dnorm(x = theta, mean = mu, sd = sqrt(var))

hist(import_rat, breaks = 20)

## Compute the importance sampled EX and VarX

ex_imp <- ((1/S) * sum(theta * import_rat)) / ((1/S) * sum(import_rat))
var_imp <- ((1/S) * sum(theta^2 * import_rat)) / ((1/S) * sum(import_rat)) - ex_imp^2

## Compute the true EX and VarX

ex_true <- 0
var_true <- 1

## Output the results

cat("Importance Sampled Expected Value: ", round(ex_imp, 4))
cat("True Expected Value: ", round(ex_true, 4))
cat("Importance Sampled Variance: ", round(var_imp, 4))
cat("True Variance: ", round(var_true, 4))



### Part d ###



## Compute the effective sample size

S_eff <- 1 / (sum((import_rat / sum(sample(import_rat)))^2))

## Output the results

cat("Effective Sample Size: ", round(S_eff, 4))