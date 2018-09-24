## No scientific notation

options(scipen = 999)



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
    
    y_rep[s,] = rexp(n = n, rate = samp_1i[s])
    
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
    
    z_rep[s,] = rnorm(n = n, mean = post_mu[s], sd = post_sigma[s])
    
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

wait_10 <- mean(z_rep < 10)

## Output the results

cat("Probability 10 Minute Wait: ", round(wait_10, 4))





