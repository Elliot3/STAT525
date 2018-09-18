## Load necessary packages

library(coda)
library(asbio)



########## Problem 1 ##########



##### Part i #####



num_draws <- 10^4

mean_c <- 1.013
mean_t <- 1.173

## Draw from the scaled inverse chi-square

n_1 <- 32
df <- n_1 - 1

samp_i_1 <- rinvchisq(num_draws, df, scale=1/df)

## Draw from the normal distribution

samp_i_2 <- rnorm(num_draws, mean = mean_c, sqrt(samp_i_1))



##### Part ii #####

## Draw from the scaled inverse chi-square

n_2 <- 36
df <- n_2 - 1

samp_ii_1 <- rinvchisq(num_draws, df, scale=1/df)

## Draw from the normal distribution

samp_ii_2 <- rnorm(num_draws, mean = mean_t, sqrt(samp_ii_1))



##### Part iii #####

## Generate the histogram of mean differences

hist(samp_ii_2 - samp_i_2,
     xlab = "Mean Difference",
     ylab = "",
     main ="Histogram: Mean Difference")



##### Part iv #####

## Construct the HPD interval of mean differences

interval <- as.vector(HPDinterval(as.mcmc(samp_ii_2 - samp_ii_1), prob = 0.95)[1, 1:2])

## Output the result

cat("HPD Interval: ", "(", interval[1], ",", interval[2], ")")



