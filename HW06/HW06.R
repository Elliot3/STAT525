########## Problem 1 ##########



##### Part i #####



##### Part ii #####



##### Part iii #####



##### Part iv #####



########## Problem 2 ##########



##### Part i #####



##### Part ii #####



nu <- 4
S <- 10000

post_gamma <- numeric()

for (s in 1:S) {
    
    post_gamma[s] <- rgamma(n = 1, shape = (nu / 2), rate = (nu / 2))
    
}

post_norm <- dnorm(x = post_gamma, mean = 0, sd = sqrt(1 / post_gamma))

t_dens <- dt(x = post_gamma, df = nu)

par(mfrow = c(1, 2))

hist(post_norm)
hist(t_dens)



##### Part iii #####



########## Problem 3 ##########



##### Part i #####



##### Part ii #####



##### Part iii #####



##### Part iv #####



##### Part v #####













