## Load in the necessary packages

library(gamair)

## Load in the data

data("hubble")

## Define the posterior distribution

errors <- rnorm()

post_dist <- lm(y ~ x, data = hubble)
