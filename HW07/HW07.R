########## Workspace Prep ##########



########## Problem 1 ##########



## Load in the lh data

data(lh)
y <- as.numeric(lh)
Y <- length(y)

## Initialize the mu function

mu <- function(x) {
    
    dnorm(x = x, mean = 0, sd = 10)
    
}

## Initialize the phi function

phi <- function(x) {
    
    dnorm(x = x, mean = 0, sd = 10)
    
}

## Initialize the tau function

tau <- function(x) {
    
    dgamma(x = x, shape = 0.01, rate = 0.01)
    
}

## Initialize the posterior density

post_func <- function(theta) {
    
    res <- (mu(theta[1]) * phi(theta[2]) * tau(theta[3]))
    
    y_t <- numeric()
    
    for (i in 1:(Y - 1)) {
        
        y_t[i] <- dnorm(x = y[i + 1],
                        mean = theta[1] + theta[2] * (y[i] - mu(theta[1])),
                        sd = sqrt(theta[3]^-1))
        
    }
    
    y_post <- prod(y_t * res)
    
    return(y_post)
    
}

## Optimize for the posterior density

opt <- optim(par = c(1,1,1),
            fn = post_func,
            control = list(fnscale = -1),
            hessian = TRUE,
            method = "L-BFGS-B", lower = c(-10, -10), upper = c(100, 100))

## Esimtate the posterior mode

post_mode <- opt$par

## Estimate the posterior variance

post_var <- solve(-opt$hessian)



##### Part i #####



## Initialize theta

theta <- post_mode

## Propose a variance

prop_var <- 2 * post_var

## Initialize an acceptance counter

count_accept <- 0

## Store the posterior simulations

post_theta = array(0, c(S, 3))

## Loop through and generate the posterior simulations

for (s in 1:S) {
    
    ##
    
}



##### Part ii #####



##### Part iii #####



##### Part iv #####



##### Part v #####



########## Problem 2 ##########



##### Part i #####



##### Part ii #####



##### Part iii #####


