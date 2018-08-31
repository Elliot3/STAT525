
### Load in the necessary packages

library(ggplot2)
library(coda)



##### Problem 2 #####



## Set the parameter values

y <- 0
n <- 20
alpha <- 2
beta <- 20

## Set the theta values

theta <- seq(from = 0, to = 1, length.out = 1000)



### Part i ###



## Plot the posterior distribution

post_beta_dens <- dbeta(x = theta, shape1 = y + alpha, shape2 = (n - y) + beta)

ggplot() +
    geom_line(aes(x = theta, y = post_beta_dens), col = "blue") +
    geom_line(aes(x = theta, y = dbeta(x = theta, shape1 = alpha, shape2 = beta)), col = "green") +
    geom_vline(aes(xintercept = (y / n)), col = "red") +
    labs(x = "Theta", y = "Density", title = "Beta Density",
         subtitle = "Blue - Posterior Distribution   /   Green - Prior Distribution   /   Red - MLE")



### Part ii ###



## Draw a sample from my posterior distribution

post_beta_rand <- rbeta(n = 10000, shape1 = y + alpha, shape2 = (n - y) + beta)

## Compute the HPD interval

bayes_int <- as.vector(HPDinterval(as.mcmc(post_beta_rand), prob = 0.95)[1, 1:2])

## Calculate the MLE

y_hat <- (y / n)

## Compute the frequentist interval

freq_int <- c((y_hat) - (1.96 * (sqrt(y_hat * (1 - y_hat)) / 2)),
              (y_hat) + (1.96 * (sqrt(y_hat * (1 - y_hat)) / 2)))



### Part iii ###



## Compute the log odds of our posterior beta distrbution

odds_samp <- log(post_beta_rand / (1 - post_beta_rand))

## Compute the HPS interval

bayes_int_odds <- as.vector(HPDinterval(as.mcmc(odds_samp), prob = 0.95)[1, 1:2])



##### Problem 3 #####



### Part i ###



## Get the poisson parameter

lambda <- (217 + 66) / 155

## Generate poisson data with the parameter

pois_data <- rpois(1000, lambda = lambda)

## Sample mean and variance

samp_mean <- mean(pois_data)
samp_var <- var(pois_data)

## Compute the gamma parameters

find_gamma_params <- function(mu, sigma_sq) {
    
    beta <- mu / sigma_sq
    alpha <- beta * mu
    return(c(alpha, beta))
    
}

gamma_params <- find_gamma_params(samp_mean, samp_var)



### Part ii ###












