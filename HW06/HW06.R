########## Workspace Prep ##########



## Load in the necessaey packages

library(coda)




########## Problem 1 ##########



##### Part i #####



##### Part ii #####



### Unchanged code from Lecture



year = (1851:1962)
y = c(4, 5, 4, 1, 0, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6, 3, 3, 5, 4, 5, 3,
      1, 4, 4, 1, 5, 5, 3, 4, 2, 5, 2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1,
      1, 1, 3, 0, 0, 1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1,
      0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2, 3, 3, 1, 1, 2,
      1, 1, 1, 1, 2, 4, 2, 0, 0, 0, 1, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 1, 0, 1)

n = length(y)

# Hyperparameters:
alpha = 0.5; beta = 0.01

# Pick some initial values:
k = ceiling(n/2) # midpoint is CP
lambda_1 = mean(y[1:k])
lambda_2 = mean(y[(k+1):n])

# Sample from the posterior distribution:
S = 10^4

# Storage:
post_lambda_1 = array(0, c(S, 1))
post_lambda_2 = array(0, c(S, 1))
post_k = array(0, c(S, 1));

for(s in 1:S){
    # Sample lambda_1:
    lambda_1 = rgamma(n = 1,
                      shape = alpha + sum(y[1:k]),
                      rate = beta + k)
    # Sample lambda_2:
    # NOTE: assuming k < n throughout!
    lambda_2 = rgamma(n = 1,
                      shape = alpha + sum(y[(k+1):n]),
                      rate = beta + (n-k))
    # Sample k:
    log_g = cumsum(y)*log(lambda_1/lambda_2) + (1:n)*(lambda_2 - lambda_1)
    #pm<-exp(log_g)/sum(exp(log_g)); k<-min((1:n)[runif(1)<cumsum(pm)])
    k = sample(1:n, 1, prob = exp(log_g))
    
    # Store:
    post_lambda_1[s] = lambda_1
    post_lambda_2[s] = lambda_2
    post_k[s] = k
}



### Edited code from lecture, factors in beta prior



# Hyperparameters:
alpha = 0.5; beta = 0.01

## Set the beta params

b_alpha <- 0.1
b_beta <- 0.1

# Pick some initial values:
k = ceiling(n/2) # midpoint is CP
lambda_1 = mean(y[1:k])
lambda_2 = mean(y[(k+1):n])

# Sample from the posterior distribution:
S = 10^4

# Storage:
post_lambda_1_new = array(0, c(S, 1))
post_lambda_2_new = array(0, c(S, 1))
post_k_new = array(0, c(S, 1))
post_beta_new = array(0, c(S, 1))

for(s in 1:S){
    
    # Sample lambda_1:
    lambda_1 = rgamma(n = 1,
                      shape = alpha + sum(y[1:k]),
                      rate = beta + k)
    # Sample lambda_2:
    # NOTE: assuming k < n throughout!
    lambda_2 = rgamma(n = 1,
                      shape = alpha + sum(y[(k+1):n]),
                      rate = beta + (n-k))
    # Sample k:
    log_g = cumsum(y)*log(lambda_1/lambda_2) + (1:n)*(lambda_2 - lambda_1)
    #pm<-exp(log_g)/sum(exp(log_g)); k<-min((1:n)[runif(1)<cumsum(pm)])
    k = sample(1:n, 1, prob = exp(log_g))
    
    # Sample beta:
    
    beta <- rgamma(n = 1,
                   shape = (2 * alpha) + b_alpha,
                   rate = lambda_1 + lambda_2 + b_beta)
    
    # Store:
    post_lambda_1_new[s] = lambda_1
    post_lambda_2_new[s] = lambda_2
    post_k_new[s] = k
    post_beta_new[s] = beta
}



##### Part iii #####



## Compute the HPD intervals for beta = 0.01

ci_lambda_1 <- as.vector(HPDinterval(as.mcmc(post_lambda_1), prob = 0.95)[1, 1:2])
ci_lambda_2 <- as.vector(HPDinterval(as.mcmc(post_lambda_2), prob = 0.95)[1, 1:2])
ci_k = as.vector(HPDinterval(as.mcmc(post_k + 1850), prob = 0.95)[1, 1:2])

cat("HPD Interval Lambda_1 with Beta = 0.01: ", "(", ci_lambda_1[1], ",", ci_lambda_1[2], ")")
cat("HPD Interval Lambda_2 with Beta = 0.01: ", "(", ci_lambda_2[1], ",", ci_lambda_2[2], ")")
cat("HPD Interval k + 1850 with Beta = 0.01: ", "(", ci_k[1], ",", ci_k[2], ")")

## Compute the HPD intervals for beta hyperprior

ci_lambda_1_new <- as.vector(HPDinterval(as.mcmc(post_lambda_1_new), prob = 0.95)[1, 1:2])
ci_lambda_2_new <- as.vector(HPDinterval(as.mcmc(post_lambda_2_new), prob = 0.95)[1, 1:2])
ci_k_new = as.vector(HPDinterval(as.mcmc(post_k_new + 1850), prob = 0.95)[1, 1:2])

cat("HPD Interval Lambda_1 with Beta Hyperprior: ", "(", ci_lambda_1_new[1], ",", ci_lambda_1_new[2], ")")
cat("HPD Interval Lambda_2 with Beta Hyperprior: ", "(", ci_lambda_2_new[1], ",", ci_lambda_2_new[2], ")")
cat("HPD Interval k + 1850 with Beta Hyperprior: ", "(", ci_k_new[1], ",", ci_k_new[2], ")")



##### Part iv #####



par(mfrow = c(1, 1))

hist(post_beta_new, freq = FALSE)
x <- rgamma(n = S, shape = b_alpha, rate = b_beta)
lines(sort(x), y = dgamma(sort(x), shape = b_alpha, rate = b_beta), lwd = 2)
abline(v = 0.01, lwd = 5)



########## Problem 2 ##########



##### Part i #####



##### Part ii #####



nu <- 4
S <- 10000

post_gamma <- numeric()

for (s in 1:S) {
    
    post_gamma[s] <- rgamma(n = 1, shape = (nu / 2), rate = (nu / 2))
    
}

post_norm <- dnorm(x = post_gamma, mean = 0, sd = sqrt(1 / post_gamma))

t_dens <- dt(x = post_gamma, df = nu)

par(mfrow = c(1, 2))

hist(post_norm)
hist(t_dens)



##### Part iii #####



########## Problem 3 ##########



##### Part i #####



##### Part ii #####



##### Part iii #####



##### Part iv #####



##### Part v #####













