
## Load in necessary packages

library(ggplot2)



##### Problem 3 #####



## Set mean and variance

mu <- 0.6
var <- 0.09



##### Part a



## Estimate the beta parameters

est_beta_params <- function(mu, var) {
    
    alpha <- (((1 - mu) / var) - (1 / mu)) * mu^2
    beta <- alpha * ((1 / mu) - 1)
    return(round(c(alpha, beta), 2))
    
}

beta_params <- est_beta_params(mu, var)

## Establish the x-axis values

theta <- seq(0, 1, length.out = 100)

## Set the beta density

beta_density <- dbeta(x = theta,
                      shape1 = beta_params[1],
                      shape2 = beta_params[2])

beta_density <- beta_density[-length(beta_density)]
theta <- theta[-length(theta)]

## Build the density plot

ggplot() +
    geom_line(aes(x = theta, y = beta_density)) +
    labs(x = "Theta", y = "Density", title = "Beta Density")



##### Part b



## Set the sample size

n <- 1000

## Set the number of yes votes

y <- n * 0.65

## Set post alpha and beta

post_alpha <- y + beta_params[1]
post_beta <- (n - y) + beta_params[2]

## Posterior mean and variance

post_mean <- round(post_alpha / (post_alpha + post_beta), 4)
post_var <- round(sqrt((post_alpha * post_beta) /
                           ((post_alpha + post_beta)^2 * (post_alpha + post_beta + 1))), 4)

## Calculate the posterior density

post_density <- dbeta(x = theta,
                      shape1 = post_alpha,
                      shape2 = post_beta)

## Build the density plot

ggplot() +
    geom_line(aes(x = theta, y = post_density)) +
    labs(x = "Theta", y = "Density", title = "Beta Density") +
    xlim(c(0.5, 0.75))



##### Part c



## First sensitivity check



## Set mean and variance

mu <- 0.5
var <- (1/12)

## Set post alpha and beta

post_alpha <- y + 1
post_beta <- (n - y) + 1

## Establish the x-axis values

theta <- seq(0, 1, length.out = 100)

## Set the beta density

beta_density <- dbeta(x = theta,
                      shape1 = post_alpha,
                      shape2 = post_beta)

beta_density <- beta_density[-length(beta_density)]
theta <- theta[-length(theta)]

## Build the density plot

ggplot() +
    geom_line(aes(x = theta, y = beta_density)) +
    labs(x = "Theta", y = "Density", title = "Beta Density") +
    xlim(c(0.5, 0.75))



## Second sensitivity check



## Set mean and variance

mu <- 0.3
var <- 0.001

## Estimate the beta parameters

beta_params <- est_beta_params(mu, var)

## Set post alpha and beta

post_alpha <- y + beta_params[1]
post_beta <- (n - y) + beta_params[2]

## Establish the x-axis values

theta <- seq(0, 1, length.out = 100)

## Set the beta density

beta_density <- dbeta(x = theta,
                      shape1 = post_alpha,
                      shape2 = post_beta)

beta_density <- beta_density[-length(beta_density)]
theta <- theta[-length(theta)]

## Build the density plot

ggplot() +
    geom_line(aes(x = theta, y = beta_density)) +
    labs(x = "Theta", y = "Density", title = "Beta Density") +
    xlim(c(0.5, 0.75))



## Third sensitivity check



## Set mean and variance

mu <- 0.9
var <- 0.0005

## Estimate the beta parameters

beta_params <- est_beta_params(mu, var)

## Set post alpha and beta

post_alpha <- y + beta_params[1]
post_beta <- (n - y) + beta_params[2]

## Establish the x-axis values

theta <- seq(0, 1, length.out = 100)

## Set the beta density

beta_density <- dbeta(x = theta,
                      shape1 = post_alpha,
                      shape2 = post_beta)

beta_density <- beta_density[-length(beta_density)]
theta <- theta[-length(theta)]

## Build the density plot

ggplot() +
    geom_line(aes(x = theta, y = beta_density)) +
    labs(x = "Theta", y = "Density", title = "Beta Density") +
    xlim(c(0.5, 0.75))

