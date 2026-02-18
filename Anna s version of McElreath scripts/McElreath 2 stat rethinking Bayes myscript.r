## CODE TO GO WITH BOOKCLUB
## ON Richard McElreath's 'Statistical Rethinking'
# Anna Dornhaus, 2025

########### PART II: Regressions #####################

## LIBRARIES  -----------------------------
library(rethinking)
library(scales)

## WHAT THIS SCRIPT CONTAINS ---------------------
## This script is the first part of chapter 4 and video 3,
## 2nd edition of the book, and videos from 2023
## This part is about starting on regressions, but first just modeling the 
## mean and standard deviation of a single variable. 

## CHAPTER 4 & VIDEO 3 -----------------------------

### Discussion of Gaussian distributions ------------------
# What's with Gaussian error and the central limit theorem

# The point is, if you add a bunch of uniformly distributed
# random numbers, the *sums* you get are normally distributed.
# It is also true that for any distribution, if you sample it and calculate
# the average of the sample, the different averages you get from different sample
# sets will have a Gaussian distribution. (this is the central limit theorem)

# Here is a little simulation to show this:
# Let 1000 people take a number of random steps between 0 and 1 meters long
# each, but in a random direction forward or backward:
steps <- 50
people <- 1000
pos <- replicate(people, sum(runif(steps,-1,1)))
# The vector pos gives the positions of all the people at the end.
dens(pos, norm.comp=TRUE)
# Or:
hist(pos)
plot(density(pos))
# Voila! Here is a Gaussian (='normal') distribution. 

# Why does this work with any distribution (here we used a uniform
# distribution for the steps)?
# This is really a property of what it means to be a 'mean':
# by definition, the mean of a distribution is the result of adding all
# possible values (and then dividing by the number). So the 
# average position, after taking all possible steps from that distribution,
# is at the mean! This is what it means to be the mean.
# Which means if we take a random selection of steps, the expected
# mean is still the mean, just with some error around it.
# The error is symmetrical, because ultimately the 'steps' possible
# in each direction have to cancel out all the steps possible in the 
# other direction. 
# So we have a distribution with a mean of zero and that is symmetrical.
# There will also be a higher density around zero, because there are more ways
# to realize a sum where left and right steps approximately cancel each other 
# out than there are cases where the sum leads far away from zero -
# the latter requires a 'run' of lots of values in one direction. 

# This is even true for calculating a product (instead of a sum),
# as long as the numbers are small (close to 1):
growth <- replicate( 10000 , prod( 1 + runif(12,0,0.1) ) )
dens(growth, norm.comp=TRUE)

# Or for products of big numbers on a log scale
log.big <- replicate( 10000 , log(prod(1 + runif(12,0,0.5))) )
dens(log.big, norm.comp=TRUE )

# There is also an inferential argument for Gaussian distribution:
# it is the least information-containing way to describe a distribution
# that contains a mean and a variance. So it assumes the least amount of things.

# Moreover, as a 'machine' for estimating a mean and variance, is works well
# even for non-normal distributions.


### !Kung people height & weight data ------------------
# Example we'll be thinking about:
# Anthropological data on !Kung people growth
d <- data("Howell1")
# At first we'll only include the adults. 
d2 <- Howell1[Howell1$age>=18,]
# Summary of these data:
precis(d2)

# This is the data we might want to fit a model to, i.e. 
# estimate the parameters that describe its distribution,
# or even more precisely, derive a posterior distribution given this
# data for the general way in which population weights depend on height.

### 'Modeling' height: finding the mean and standard deviation of height -------------

# Now it might be tempting to say we can just calculate the mean and standard deviation
# of the sample. But the sample is just that: a sample of a larger population. What
# we actually want it to derive what for the entire population, or even population of 
# possible humans in these conditions, the average height is and how variable it is.

# So we are going to fit a 'model', in other words we assume a function that fits
# the data (here a normal distribution) and then fit the parameters that characterize
# that function (the mean and standard deviation). So we derive posterior probability
# distributions for the population mean for height, and a posterior probability 
# distribution for the population variation in height. Which means we now have two
# parameters and two posteriors (when before we had only one, the proportion of water).

# To derive the posterior, we need a prior and the likelihood of the data given different
# parameter values. 

# Let's first look at the data for HEIGHT
dens(d2$height)

#### Thinking in more detail about priors --------------------------

# McElreath believes in fairly informative, but still 'soft', priors.
# More importantly, he strongly advises simulating data according to the priors
# first, to see what the priors actually suggest about our current beliefs.

# For example, for height, we might assume a prior for 
# mu, average height, of 178cm +/- 20cm
curve(dnorm(x, 178, 20), from=100, to=250)
# Note that this is just the prior for the mean, and thus the 20 is basically 
# our (prior) uncertainty about the mean. 

# For sigma, the SD of height across people in the population, a uniformly 
# distributed prior that is positive and smaller than 50cm.
# Note that this is a different type of variation, namely how far individual
# people are from the mean. 
curve(dunif(x, 0, 50), from=-10, to=60) 
# (We're plotting from -10 but the prior assumes positive sigma)

# Now we can simulate data (heights) given these priors.
# That is, we generate possible distributions of heights (the priors form
# a distribution of possible distributions).
sample_mu <- rnorm(1e4, 178, 20) # 10000 possible means
sample_sigma <- runif(1e4, 0, 50) # 10000 possible sigmas
prior_h <- rnorm(1e4, sample_mu, sample_sigma) 
# 10000 possible heights drawn from 10000 combinations of possible mu and sigma
dens(prior_h)

##### Doing the Bayesian model (to later check effect of priors) --------------------------

# Grid approximation is going to be increasingly hard with complex models, 
# but we're doing it again here:
# Make a list of mu's
mu.list <- seq(from=150, to=160 , length.out=1000)
# Make a list of sigmas
sigma.list <- seq(from=7 , to=9 , length.out=1000)
# Make a list of all combinations (our 'grid' is now a 2-D grid not just a list)
post <- expand.grid(mu=mu.list, sigma=sigma.list)

# Likelihood function
# How often will you get this height given each possible mu and sigma
post$LL <- sapply(1:nrow(post), function(i) sum(
  dnorm(d2$height, post$mu[i], post$sigma[i], log=TRUE)))
# Note this can take a little while to calculate! Give it a minute. 

# Product of prior and likelihood
post$prod <- post$LL + dnorm(post$mu, 178, 20, TRUE) +
  dunif(post$sigma, 0, 50, TRUE)

# Standardization to 1
post$prob <- exp(post$prod - max(post$prod))

# There are many combinations of mu and sigma, and we just calculated a 
# posterior probability for each combination (10000 x 10000).

###### Several ways of plotting the 2-dimensional posterior --------------------------

# Four different ways of plotting:
# (1)
plot(post$prob~post$mu
     , xlab = "Possible population means"
     , ylab = "Posterior probabilities")
plot(post$prob~post$sigma
     , xlab = "Possible variation in population height around mean"
     , ylab = "Posterior probabilities")
# Both the above take a while to plot, and each only shows one of the parameters.
# Better is a 2-dimensional plot that shows the two of them together (and they
# are not independent). Both the below are contour plots, showing the height of
# the probability distribution either as elevation lines -
# (2)
contour_xyz(post$mu, post$sigma, post$prob
            , xlab = "Posterior distribution for mu"
            , ylab = "Posterior distribution for sigma")

# or showing it in colors - 
# (3)
image_xyz(post$mu, post$sigma, post$prob
          , xlim = c(152, 158)
          , ylim = c(7, 9)
          , xlab = "Posterior distribution for mu"
          , ylab = "Posterior distribution for sigma"
          )

# To understand the posterior distribution, rather than calculating
# things from it, it is useful to sample it first:
sample.rows <- sample(1:nrow(post), size=1e4, replace=TRUE ,
                       prob=post$prob )
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]

# Plot the samples:
# (4)
plot(sample.mu, sample.sigma, 
     cex=0.5, 
     pch=16 , 
     col=col.alpha(rangi2,0.1) 
     )
# (5) - this is essentially the same as (1) but from the samples rather than 
# the posterior directly - 
dens(sample.mu)
dens(sample.sigma)

# Remember that these need not be Gaussian distributions.
# In particular, since sigma must be positive, it is not symmetrical
# and thus not Gaussian, especially at small expected sigma or high uncertainty.


#######################################-----------------
