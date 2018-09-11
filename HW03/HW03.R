########## Problem 2 ##########





## Set options

options(scipen = 999)

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
s_squared <- (1 / (n - 1)) * sum(y - (x * theta_hat))^2



##### Part iii #####



## Make 10,000 draws from [tau_y|y], tau^s

draws_1 <- rgamma(n = 10000, shape = alpha + ((n - 1) / 2),
                  rate = beta + (n - 1) * (s_squared / 2))

## Make 10,000 draws from [theta|y, tau_y = t^s]

var_y <- 1 / draws_1

draws_2 <- rnorm(n = 10000,
                 mean = ((1/var_y) * sum(x^2))^(-1) * (1/var_y) * sum(x * y),
                 sd = sqrt(((1/var_y) * sum(x^2))^(-1))
                 )

## Print out a sample of the draws

draw_samp <- sample(x = draws_2, size = 10)

cat("Sample Draw:\n\n", draw_samp)



##### Part iv #####



## Construct the HPD interval

post_hub <- 1/draws_2*3.09e19/(60^2*24*365)

hpd_int <- as.vector(HPDinterval(as.mcmc(post_hub), prob = 0.95))

cat("HPD Interval: ", "(", hpd_int[1], ",", hpd_int[2], ")")
