
########## Prepare Workspace ##########



## Set the working directory

setwd('~/Documents/Rice_University/Fall_2018/STAT525/HW08')

## Load in the necessary packages

suppressMessages(
    suppressWarnings(
        library(rstan)
    )
)
suppressMessages(
    suppressWarnings(
        library(coda)
    )
)
suppressMessages(
    suppressWarnings(
        library(bayesplot)
    )
)

## Detect the number of core for parallel processing

options(mc.cores = parallel::detectCores())



########## Problem 1 ##########



## Set the known parameters

n <- 42
y <- 10

## Compile both of the stan models

mod_binom_beta <- stan_model('bayes_binom.stan')
mod_binom_unif <- stan_model('bayes_binom_eps.stan')

## Store the data in a list

dat_beta <- list(n = n, y = y, alpha = 1, beta = 1)
dat_unif <- list(n = n, y = y, alpha = 0, beta = 1) 

## Perform the sampling from the beta prior posterior

fit_beta <- sampling(object = mod_binom_beta,
                     data = dat_beta,
                     iter = 1000, chains = 4)

## Perform the sampling from the uniform prior posterior

fit_unif <- sampling(object = mod_binom_unif,
                     data = dat_unif,
                     iter = 1000, chains = 4)

## Print out the results for beta prior and uniform prior

print(fit_beta)
print(fit_unif)

## Run diagnostic tests for beta prior and uniform prior

check_hmc_diagnostics(fit_beta)
check_hmc_diagnostics(fit_unif)

stan_rhat(fit_beta, 'theta', bins = 10)
stan_rhat(fit_unif, 'theta', bins = 10)

## Plot the traceplot for each prior

rstan::traceplot(fit_beta, pars = c("theta"))
rstan::traceplot(fit_unif, pars = c("theta"))

## Extract the posterior draws of each of the thetas

post_theta_beta <- rstan::extract(fit_beta, "theta", permuted = TRUE)
post_theta_unif = rstan::extract(fit_unif, "theta", permuted = TRUE)

## Plot each theta in the "conventional" way

plot(as.mcmc(as.matrix(post_theta_beta[[1]])))
plot(as.mcmc(as.matrix(post_theta_unif[[1]])))

## Plot the posterior distribution for each prior

mcmc_areas(as.matrix(fit_beta), pars = c('theta'), prob = 0.8) +
    ggtitle("Posterior distribution", "with median and 80% intervals")

mcmc_areas(as.matrix(fit_unif), pars = c('theta'), prob = 0.8) +
    ggtitle("Posterior distribution", "with median and 80% intervals")
