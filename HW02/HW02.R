
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
    geom_line(aes(x = theta, y = post_beta_dens), linetype = 1) +
    geom_line(aes(x = theta, y = dbeta(x = theta, shape1 = alpha, shape2 = beta)), linetype = 2) +
    geom_vline(aes(xintercept = (y / n)), linetype = 3) +
    labs(x = "Theta", y = "Density", title = "Beta Density",
         subtitle = "Smooth Line - Posterior Distribution   /   Dashed Line - Prior Distribution   /   
         Dotted Line - MLE")



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



## Sample mean and variance

samp_mean <- ((1 * 21) + (2 * 43) + (3 * 23) + (4 * 13)) / (21 + 43 + 23 + 13)
samp_var <- var(c(rep(1, 21), rep(2, 43), rep(3, 23), rep(4, 13)))

## Compute the gamma parameters

find_gamma_params <- function(mu, sigma_sq) {
    
    beta <- mu / sigma_sq
    alpha <- beta * mu
    return(c(alpha, beta))
    
}

gamma_params <- round(find_gamma_params(samp_mean, samp_var), 4)



### Part ii ###



## Set theta values

theta <- seq(from = 1, to = 5, by = 0.01)

## Set the sample values

y_1 <- 217
y_2 <- 66
n_1 <- 111
n_2 <- 44
lambda_1 <- 217 / 111
lambda_2 <- 66 / 44

## Generate samples from the posterior distributions

post_dens_y1 <- dgamma(x = theta, shape = y_1 + gamma_params[1], rate = gamma_params[2] + n_1)
post_dens_y2 <- dgamma(x = theta, shape = y_2 + gamma_params[1], rate = gamma_params[2] + n_2)
prior_dens <- dgamma(x = theta, shape = gamma_params[1], rate = gamma_params[2])

## Generate the plot

ggplot() +
    geom_line(aes(x = theta, y = post_dens_y1), linetype = 1) +
    geom_line(aes(x = theta, y = post_dens_y2), linetype = 2) +
    geom_line(aes(x = theta, y = prior_dens), linetype = 3) +
    geom_vline(aes(xintercept = lambda_1), linetype = 4) +
    geom_vline(aes(xintercept = lambda_2), linetype = 6) +
    labs(x = "Theta", y = "Density", title = "Gamma Density",
         subtitle = "Smooth Line - y_1 Posterior   /   Dashed Line - y_2 Posterior   /   
         Dotted/Short Dash - y_1 MLE   /   Dotted/Long Dash - y_2 MLE   /   
         Dotted Line - Prior")



### Part iii ###



## Draw a sample from each posterior distribution

post_gamma_y1 <- rgamma(n = 10000, shape = y_1 + gamma_params[1], rate = n_1 + gamma_params[2])
post_gamma_y2 <- rgamma(n = 10000, shape = y_2 + gamma_params[1], rate = n_2 + gamma_params[2])

## Compute the HPD intervals

hpd_int_y1 <- as.vector(HPDinterval(as.mcmc(post_gamma_y1), prob = 0.95)[1, 1:2])
hpd_int_y2 <- as.vector(HPDinterval(as.mcmc(post_gamma_y2), prob = 0.95)[1, 1:2])

## Output the results

cat("Posterior Credible Interval, theta_1: ", "(", hpd_int_y1[1], ",", hpd_int_y1[2], ")")
cat("Posterior Credible Interval, theta_2: ", "(", hpd_int_y2[1], ",", hpd_int_y2[2], ")")



### Part iv ###



## Get the difference of our random samples

post_diff_gamma <- post_gamma_y1 - post_gamma_y2

## Compute the HPD interval

diff_int <- as.vector(HPDinterval(as.mcmc(post_diff_gamma), prob = 0.95)[1, 1:2])

## Output the results

cat("Posterior Credible Interval, theta_1 - theta_2: ", "(", diff_int[1], ",", diff_int[2], ")")



### Part v ###



post_prob <- mean(post_diff_gamma > 0)