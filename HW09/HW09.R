########## Workspace Prep ##########



## Load in the necessary packages

library(ggplot2)
suppressMessages(
    suppressWarnings(
        library(truncdist)
    )
)
library(coda)
library(tidyr)
library(knitr)
suppressWarnings(
    library(LearnBayes)
)

## No scientific notation

options(scipen = 999)



########## Problem 1 ##########



## Load in the data

data("hearttransplants")
H = length(hearttransplants$y)
y = hearttransplants$y
x = hearttransplants$e



##### Part i #####



##### Part ii #####



## Set the number of iterations

S <- 10000

## Initialize values

mu <- 2
alpha <- 1
z_0 <- 0.5

## Containers for the MCMC output

post_lambda_h <- array(0, c(S, H))
post_mu <- numeric(S)
post_alpha <- numeric(S)
post_kappa <- array(0, c(S, H))

## Run the MCMC

for (s in 1:S) {
    
    mu <- rgamma(n = 1, shape = y + 1, rate = 1)
    
    alpha <- rgamma(n = 1, shape = y + 1, rate = 1) * (z_0 / (z_0 + alpha)^2)
    
    lambda_h <- rgamma(n = H,
                       shape = alpha + y,
                       rate = (alpha / mu) + x)
    
    post_alpha[s] <- alpha
    post_mu[s] <- mu
    post_lambda_h[s,] <- lambda_h
    post_kappa[s,] = alpha / (alpha + (x * mu))
    
}

## Scale up post_kappa

post_kappa <- post_kappa * 100



##### Part iii #####



## Construct the trace plots

plot(as.mcmc(post_mu))
plot(as.mcmc(post_alpha))
plot(as.mcmc(post_lambda_h[,1]))
plot(as.mcmc(post_lambda_h[,80]))

## Calculate the effective sample sizes

eff_mu <- as.vector(effectiveSize(post_mu))
eff_alpha <- as.vector(effectiveSize(post_alpha))
eff_lambda_1 <- as.vector(effectiveSize(post_lambda_h[,1]))
eff_lambda_80 <- as.vector(effectiveSize(post_lambda_h[,80]))

## Output the results

cat("Effective Size, post_mu: ", eff_mu)
cat("Effective Size, post_alpha: ", eff_alpha)
cat("Effective Size, post_lambda, h = 1: ", eff_lambda_1)
cat("Effective Size, post_lambda, h = 80: ", eff_lambda_80)



##### Part iv #####



## Calculate the log-exposure

log_x <- log(x)

## Plot posterior means as a function of log_x

plot(x = log_x, y = colMeans(post_kappa), ylim = c(0,1),
     xlab = 'Log Exposure', ylab = 'Post_Kappa')



##### Part v #####



## Generate the plot from class

yhat_mle = (y / x)
yhat_mle_pool = (sum(y) / sum(x))

ci_h <- array(0, c(H, 3))

for (i in 1:H) {
    
    ci_h[i,1:2] <- round(as.vector(HPDinterval(as.mcmc(post_lambda_h[i,]), prob = 0.95)[1, 1:2]), 6)
    ci_h[i,3] <- round(mean(post_lambda_h[i,]), 6)
    
}

plot(log(x), yhat_mle, ylim = range(ci_h, (y /  x)),
     pch = as.character(y), cex = 1,
     main = 'Death rates by (log) Exposure')
abline(h = yhat_mle_pool, lwd=4, col='green', lty=2)

for(h in 1:H) {
    
    lines(rep(log_x[h], 2),
          ci_h[h,1:2], col='blue', lwd=1)
    
    lines(log_x[h], ci_h[h,3], 
          type='p', pch=4, lwd=1, cex = 2, col='red')
    
}



########## Problem 2 ##########



## Input the known data and parameters

y_j <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma_j <- c(15, 10, 16, 11, 9, 11, 10, 18)
sigma_j_sq <- sigma_j^2
names(y_j) <- c("A", "B", "C", "D", "E","F", "G", "H")
names(sigma_j) <- c("A", "B", "C", "D", "E","F", "G", "H")

J <- length(y_j)



##### Part i #####



### A = 100 ###



## Set some intial values

theta_j <- y_j
mu <- mean(theta_j)
sigma_theta <- sd(y_j)

A <- 100
S <- 10000

## Create the containers for the posterior draws

post_theta_j <- array(0, c(S, J))
post_mu <- numeric(S)
post_sig_theta <- numeric(S)
post_kappa_j <- array(0, c(S,J))

## Run the MCMC

for(s in 1:S) {
    
    ## Sample from mu
    
    Q_mu <- sum(1/(sigma_j^2 + sigma_theta^2))
    ell_mu <- sum(y_j/(sigma_j^2 + sigma_theta^2))
    mu <- rnorm(n = 1,
                mean = Q_mu^-1*ell_mu,
                sd = sqrt(Q_mu^-1))
    
    ## Sample from theta_j
    
    Q_theta <- 1/sigma_j^2 + 1/sigma_theta^2
    ell_theta <- y_j/sigma_j^2 + mu/sigma_theta^2
    theta_j <- rnorm(n = J,
                     mean = Q_theta^-1*ell_theta,
                     sd = sqrt(Q_theta^-1))
    
    ## Sample from sigma_theta
    
    eta_theta <- rtrunc(n = 1, 
                        'gamma',
                        a = 1/A^2,
                        b = A, 
                        shape = J/2 - 1/2,
                        rate =  sum((theta_j-mu)^2)/2)
    sigma_theta <- 1/sqrt(eta_theta)
    
    post_mu[s] <- mu
    post_theta_j[s,] <- theta_j
    post_sig_theta[s] <- sigma_theta
    
    # Get the shrinkage parameter
    
    post_kappa_j[s,] <- sigma_j^2/(sigma_theta^2 + sigma_j^2)
    
}

## Function for V inverse

v_inv <- function(sigma_j_sq, sigma_theta) {
    
    sum(1 / (sigma_j_sq + sigma_theta^2))
    
}

## Function for mu hat

mu_hat <- function(y_j, sigma_j_sq, sigma_theta) {
    
    sum(1 / (sigma_j_sq + sigma_theta^2) * y_j) / v_inv(sigma_j_sq, sigma_theta)
    
}

## Function for post sigma_theta

post_sigma_theta <- function(y_j, sigma_j_sq, sigma_theta, mu_h) {
    
    v_inv(sigma_j_sq, sigma_theta)^(-(1 / 2)) * prod((sigma_j_sq + sigma_theta^2)^(-(1 / 2) * exp(-(y_j - mu_h)^2 / (2 * sigma_j_sq + sigma_theta^2))))
    
}

## Set a sequence for sigma_theta

seq_sigma_theta <- seq(from = 0, to = 30, length.out = S)

## Set a container for plotting the posterior of sigma_theta

post_sigma_theta_plot <- numeric()

## Loop to sample from post sigma_theta and plot the results

for (i in 1:length(seq_sigma_theta)) {
    
    mu_h <- mu_hat(y_j, sigma_j_sq, seq_sigma_theta[i])
    post_sigma_theta_plot[i] <- post_sigma_theta(y_j, sigma_j_sq, seq_sigma_theta[i], mu_h)
    
}

plot(x = seq_sigma_theta, y = post_sigma_theta_plot, type = 'l',
     xlab = 'sigma_theta', ylab = 'p(sigma_theta|y)', main = 'Marginal Posterior Density')

## Plot the posterior draws of sigma_theta

hist(post_sig_theta)

## Plot the posterior draws of mu

labels = c("A", "B", "C", "D", "E", "F", "G", "H")

ggplot(gather(data.frame(post_theta_j)), aes(value, linetype=key)) + geom_freqpoly(binwidth=2) + 
    scale_linetype_discrete(name = "School", labels = labels)

par(mfrow = c(2, 4))

for (col in 1:ncol(post_theta_j)) {
    
    hist(post_theta_j[,col], main = paste("Hist of School" , labels[col]))
    
}

par(mfrow = c(1, 1))

## Calculate the probability table

max_vals <- apply(post_theta_j, 1, which.max)
prob_max <- table(max_vals) / S

prob_table <- matrix(NA, J, J)

for (i in 1:(J-1)) {
    
    for (j in (i+1):J) {
        
        prob_table[i,j] = sum(post_theta_j[,i] > post_theta_j[,j]) / S
        prob_table[j,i] = 1 - prob_table[i,j]
        
    }
    
}

prob_table[is.na(prob_table)] = "---"
schools = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(prob_table) <- schools
kable(cbind(schools, prob_max, prob_table))



### sigma_theta = Infinity ###



## Set some intial values

theta_j <- y_j
mu <- mean(theta_j)
sigma_theta <- sd(y_j)

S <- 10000

## Create the containers for the posterior draws

post_theta_j <- array(0, c(S, J))
post_mu <- numeric(S)
post_sig_theta <- numeric(S)
post_kappa_j <- array(0, c(S,J))

## Run the MCMC

for(s in 1:S) {
    
    ## Sample from mu
    
    Q_mu <- sum(1/(sigma_j^2 + sigma_theta^2))
    ell_mu <- sum(y_j/(sigma_j^2 + sigma_theta^2))
    mu <- rnorm(n = 1,
                mean = Q_mu^-1*ell_mu,
                sd = sqrt(Q_mu^-1))
    
    ## Sample from theta_j
    
    Q_theta <- 1/sigma_j^2 + 1/sigma_theta^2
    ell_theta <- y_j/sigma_j^2 + mu/sigma_theta^2
    theta_j <- rnorm(n = J,
                     mean = Q_theta^-1*ell_theta,
                     sd = sqrt(Q_theta^-1))
    
    ## Sample from sigma_theta
    
    # eta_theta <- rtrunc(n = 1, 
    #                     'gamma',
    #                     a = 1/A^2,
    #                     b = A, 
    #                     shape = J/2 - 1/2,
    #                     rate =  sum((theta_j-mu)^2)/2)
    # sigma_theta <- 1/sqrt(eta_theta)
    
    sigma_theta <- 99999
    
    post_mu[s] <- mu
    post_theta_j[s,] <- theta_j
    post_sig_theta[s] <- sigma_theta
    
    # Get the shrinkage parameter
    
    post_kappa_j[s,] <- sigma_j^2/(sigma_theta^2 + sigma_j^2)
    
}

## Function for V inverse

v_inv <- function(sigma_j_sq, sigma_theta) {
    
    sum(1 / (sigma_j_sq + sigma_theta^2))
    
}

## Function for mu hat

mu_hat <- function(y_j, sigma_j_sq, sigma_theta) {
    
    sum(1 / (sigma_j_sq + sigma_theta^2) * y_j) / v_inv(sigma_j_sq, sigma_theta)
    
}

## Function for post sigma_theta

post_sigma_theta <- function(y_j, sigma_j_sq, sigma_theta, mu_h) {
    
    v_inv(sigma_j_sq, sigma_theta)^(-(1 / 2)) * prod((sigma_j_sq + sigma_theta^2)^(-(1 / 2) * exp(-(y_j - mu_h)^2 / (2 * sigma_j_sq + sigma_theta^2))))
    
}

## Plot the posterior draws of mu

labels = c("A", "B", "C", "D", "E", "F", "G", "H")

ggplot(gather(data.frame(post_theta_j)), aes(value, linetype=key)) + geom_freqpoly(binwidth=2) + 
    scale_linetype_discrete(name = "School", labels = labels)

par(mfrow = c(2, 4))

for (col in 1:ncol(post_theta_j)) {
    
    hist(post_theta_j[,col], main = paste("Hist of School" , labels[col]))
    
}

par(mfrow = c(1, 1))

## Calculate the probability table

max_vals <- apply(post_theta_j, 1, which.max)
prob_max <- table(max_vals) / S

prob_table <- matrix(NA, J, J)

for (i in 1:(J-1)) {
    
    for (j in (i+1):J) {
        
        prob_table[i,j] = sum(post_theta_j[,i] > post_theta_j[,j]) / S
        prob_table[j,i] = 1 - prob_table[i,j]
        
    }
    
}

prob_table[is.na(prob_table)] = "---"
schools = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(prob_table) <- schools
kable(cbind(schools, prob_max, prob_table))



##### Part ii #####



##### Part iii #####