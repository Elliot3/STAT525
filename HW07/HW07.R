########## Workspace Prep ##########



########## Problem 1 ##########



## Load in the lh data

data(lh)
y <- as.numeric(lh)
Y <- length(y)



##### Part i #####



## Initialize an acceptance counter

count_accept <- 0

## Store the posterior simulations

post_theta = array(0, c(S, 3))

## Loop through and generate the posterior simulations

for (s in 1:S) {
    
    mu_star <- rnorm(n = 1, mean = 0, sd = 10)
    phi_star <- rnorm(n = 1, mean = 0, sd = 10)
    tau_star <- rgamma(n = 1, shape = 0.01, rate = 0.01)
    
    theta_star <- mu_star * phi_star * tau_star
    
}



##### Part ii #####



##### Part iii #####



##### Part iv #####



##### Part v #####



########## Problem 2 ##########



##### Part i #####



##### Part ii #####



##### Part iii #####


