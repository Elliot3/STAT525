## Load in the necessary packages

library(gamair)
library(coda)

## Load in the data

data("hubble")

x <- hubble$x
y <- hubble$y

## Setup the parameters

alpha <- 0.01
beta <- 0.01
n <- 24
theta_hat <- sum(x * y) / sum(x^2)
theta_hat <- 1/theta_hat*3.09e19/(60^2*24*365)
s_squared <- (1 / (n - 1)) * sum(y - (x * theta_hat))^2

## Make 10,000 draws from [tau_y|y], tau^s

draws_1 <- rgamma(n = 10000, shape = alpha + (n - 1), rate = beta + (n - 1) * (s_squared / 2))

## Make 10,000 draws from [theta|y, tau_y = t^s]

var_y <- var(y)

draws_2 <- rnorm(n = 10000,
                 mean = ((1/var_y) * sum(x^2))^(-1) * ((1/var_y) * sum(x * y))^(-1),
                 sd = ((1/var_y) * sum(x^2))^(-1)
                 )

## Construct the HPD interval

HPDinterval(as.mcmc(draws_2), prob = 0.95)
