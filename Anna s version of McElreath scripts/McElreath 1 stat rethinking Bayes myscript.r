## CODE TO GO WITH BOOKCLUB
## ON Richard McElreath's 'Statistical Rethinking'
# Anna Dornhaus, 2025

########### PART I: Bayes, Posteriors, Sampling, Predicting #####################

## INSTALLATION -----------------------------
# - does not need to be repeated after first time, which is why I've commented out
# the install commands

# Install rstan first, for which you need rtools

library(rstan)
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

## Then install 'rethinking' package for the book-specific functions.

# From 2019 book:
#install.packages(c("coda","mvtnorm","devtools"))
library(devtools)
#devtools::install_github("rmcelreath/rethinking",ref="Experimental")
#devtools::install_github("rmcelreath/rethinking", force=TRUE)

# From R Documentation:
#devtools::install_github("stan-dev/cmdstanr")
#install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
#devtools::install_github("rmcelreath/rethinking")

## LIBRARIES  -----------------------------

# Do I need these? I don't think so
#library(cmdstanr)
#library(posterior)
#library(bayesplot)
#color_scheme_set("brightblue")


library(rethinking)
library(scales)

## WHAT THIS SCRIPT CONTAINS ---------------------
## This script is roughly chapters 1-3 and videos 1 & 2 from 2023,
## 2nd edition of the book. 
## All about globe tossing

## VIDEO 1 & 2 and roughly chapter 1 & 2: globe tossing ---------------

# Notes on his original code on github:
# Garden of forking data graphics work with minor correction in first function call
# Globe tossing scripts don't work because it can't find
# function 'quartz', likely because it is a macOS specific graphics driver

# Note that I didn't code the '4 sided globe'. This is the general function for a real globe.

### WHAT IS BAYES? Calculating a posterior probability -----------
# Calculating the likelihood of the data
# if prob is the probability of the outcome 'success' 
# on each sample
no_of_samples <- 9
no_of_successes <- 6
# (6 W in 9 samples)

# Prefix d in R means probability density function. The “d” in dbinom stands 
# for density. Functions named in this way almost always have corresponding 
# partners that begin with “r” for random samples and that begin with “p” 
# for cumulative probabilities. (i.e. here rbinom(), pbinom() )

# Here this is the probability density for a binomial function
# which is what you need when drawing with replacement
dbinom(no_of_successes, size=no_of_samples, prob=0.5)

## QUOTE:
# "And this is Bayes’ theorem: It says that the probability of any particular value of p, considering
# the data, is equal to the product of the relative plausibility of the data, conditional on
# p, and the prior plausibility of p, divided by this thing Pr(W, L), which I’ll call the average
# probability of the data. In word form:
#    Posterior =
#    Probability of the data × Prior / Average probability of the data   "

# The thing we are dividing by is just to standardize the integral over the posterior to 1.
# QUOTE AGAIN:
# "The average probability of the data, Pr(W, L), can be confusing. It is commonly called the
# “evidence” or the “average likelihood,” neither of which is a transparent name. The probability
# Pr(W, L) is literally the average probability of the data. Averaged over what? Averaged
# over the prior."
# This standardization is actually not terribly important, so if it is confusing, feel free
# to forget about it now; the important point is that probability of data / average probability
# of data is the 'likelihood'. 

## But in general, this is Bayes' theorem. 
# Why? Because for each specific value of p, the number of paths
# through the garden of forking data is the product of the prior number of paths and the new
# number of paths. Multiplication is just compressed counting. The average probability on
# the bottom just standardizes the counts so they sum to one.

## 'Motors' are the numerical techniques that approximate this, because
# in most real applied situations the above cannot be exactly computed/derived.

# E.g. grid approximation, quadratic approximation, MCMC are such 'motors'. 

## Note that most of the first part of the book just uses the quadratic approximation
# function from the rethinking package, quap(). 

# So probability of data/avg prob of data is a measure of surprise -
# The relative number of ways that a value p can produce the data is usually called
# a likelihood. It is derived by enumerating all the possible data sequences that
# could have happened and then eliminating those sequences inconsistent with the
# data.

### CODING THE GLOBE TOSSING PROBLEM with three motors ---------------------- 
# We define three functions to estimate the posterior in three ways.

#### (1) Simple Grid approximation: -------------
grid_approx <- function(W, L, grid_size=10) {
# define grid of possible hypotheses/values of thing to be estimated
p_grid <- seq(from=0 , to=1 , length.out=grid_size)
# define prior across all hypotheses
prior <- rep(1 , grid_size)
# Alternative priors:
#prior <- ifelse( p_grid < 0.5 , 0 , 1 )
#prior <- exp( -5*abs( p_grid - 0.5 ) )
# Standardizing just for the sake of the graph
prior <- prior/sum(prior)
# compute likelihood at each value in grid (i.e. data given each hypothesis)
# Note that here we have a clear mathematical function to do this; usually we don't
likelihood <- dbinom(W, size=W+L, prob=p_grid )
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)
# Note that there is a bit of a sleight of hand here in that before he called the ratio
# not just the numerator the 'likelihood'. Don't worry about it if you didn't notice.

# Plot all three
plot(p_grid , posterior 
     , xlab="Proportion of water on globe" 
     , ylab="Posterior probability"
     , ylim=c(0,1)
     , col = "black"
     , pch = 19
     , yaxt='n'
     #, type="b" 
)
points(p_grid, prior, type="b" , col = "blue")
points(p_grid, likelihood, type="b" , col = "red")

return(posterior)
}
# Note that we thus defined a function that doesn't yet do anything, but it takes
# the number of waters (W) and lands (L) we got as arguments, and then does both
# the grid-based calculation of a posterior (which it will output to the console)
# and the corresponding plot. 

# Let's try it:
grid_approx(5, 3)
# You should see the prior (blue), the likelihood (red) and the posterior (black dots).


#### (2) Quadratic approximation using quap(): -----------------
globe.qa <- function(W, L) {
  summary_quadratic_approx <- 
    # precis() is his 'summarize model results' function
    precis(
      # quap() is his actual model fitting function, which takes a list of model 
      # assumptions and the data as its two arguments. 
    quap(
      alist(
        W ~ dbinom(W+L ,p) , # binomial likelihood
        p ~ dunif(0,1) # uniform prior
        ),
      data=list(W=W,L=L)
      )
  )
  # display the summary of quadratic approximation
  summary_quadratic_approx
  
  # Then add curves/lines to the presumed existing graph of prior and likelihood.
  # precis() returns the mean, sd, 5.5% and 94.5% percentile of the posterior.
  # So we can plot the mean of the posterior as a vertical line:
  abline(v=summary_quadratic_approx[1,1], col= alpha("seagreen", 0.5), lwd = 4)
  # And/or we can plot the interval from the lower to the higher percentile:
  segments(summary_quadratic_approx[1,3]
           ,0.5
           ,summary_quadratic_approx[1,4]
           ,col="seagreen1"
           ,lwd = 3
  )
  # And/or we can plot just the interval from mean-sd to mean+sd:
  segments(summary_quadratic_approx[1,1]-summary_quadratic_approx[1,2]
           ,0.5
           ,summary_quadratic_approx[1,1]+summary_quadratic_approx[1,2]
           ,col="seagreen3"
           ,lwd = 6
  )
  # Or, we plot an estimate of the entire posterior distribution by fitting
  # a Gaussian function with the same mean and standard deviation:
  curve((dnorm(x
               , summary_quadratic_approx[1,1] 
               , summary_quadratic_approx[1,2])
         )/max(dnorm(x, summary_quadratic_approx[1,1]
                     , summary_quadratic_approx[1,2])
               ) 
        , add=TRUE 
        , col="black"
        , lwd = 2)
  
  # Analytical calculation from book - note that the entire 
  # point of the approximations above is that this cannot usually be done
  # for real problems
  curve( (dbeta( x , W+1 , L+1 ))/max(dbeta( x , W+1 , L+1 )) , from=0 , to=1 , add=TRUE, lty=2, col = "purple")
  
  # Note that for both 'curves' I normalized by their maximum, so y-values
  # do not correspond to actual pdf but the peak is standardized at y=1
}

# Let's try it:
# (Note this only runs if you did in fact plot something previously with the grid function above)
globe.qa(5, 3)
# The purple line is the (cheating) mathematical known approximation of the posterior
# (hard to see because essentially exactly on the black).
# The total height of the curves is not informative here, because I just standardized
# the maximum to 1 to be able to fit all three approximation methods into the same graph.
# The black line is the approximated, Gaussian posterior. 
# The green segments show the mean, mean+/-SD, and the 5.5%-94.5% interval of the 
# estimated posterior from the precis() output. 

#### (3) MCMC approximation: --------------------------
# This is on page 46, the end of chapter 2. He just shows it here but doesn't explain
# how or why it works; that's for later. 
globe_mcmc <- function(W, L, n_samples=1000, n_mcmc_reps=5) {
  for (j in 1:n_mcmc_reps) {
  mcmc_points <- rep( NA , n_samples )
  mcmc_points[1] <- 0.5
  for ( i in 2:n_samples ) {
    p_new <- rnorm( 1 , mcmc_points[i-1] , 0.1 )
    if ( p_new < 0 ) p_new <- abs( p_new )
    if ( p_new > 1 ) p_new <- 2 - p_new
    q0 <- dbinom( W , W+L , mcmc_points[i-1] )
    q1 <- dbinom( W , W+L , p_new )
    mcmc_points[i] <- ifelse( runif(1) < q1/q0 , p_new , mcmc_points[i-1] )
  }
  # Plot the pdf as its own graph
  #dens( p , xlim=c(0,1) )
  # Or add into previous graph
  par(new=TRUE)
  dens(mcmc_points
       , xlim=c(0,1)
       , lwd=3
       , yaxt = 'n'
       , xaxt = 'n'
       , col=adjustcolor("grey", alpha.f = 0.4)
       , ann=FALSE)
  # Ok graphical design can be improved, again y-axis scale is not 
  # correct for probability density function, just to be able to overlap with prior graph.
  }
  return(mcmc_points)
}

#### Graph results of all methods -------------------------------

# Now that we've defined all three methods ('motors') as functions, we can play around
# with different data to our heart's content and plot the outcome of all three
# each time. 

# Let's say this is the data:
no_of_samples <- 5
observations_of_water <- 3
# And these are our decisions
grid_size <- 10 # for grid approximation
n_samples_mcmc <- 1000 # for MCMC
n_mcmc_reps <- 5 # for MCMC

## Graph everything in one graph
grid_approx(W=observations_of_water, L=no_of_samples - observations_of_water, grid_size=10)
globe.qa(W=observations_of_water, L=no_of_samples - observations_of_water) 
globe_mcmc(W=observations_of_water, L=no_of_samples - observations_of_water, n_samples_mcmc, n_mcmc_reps)
mtext(paste(grid_size, "points in grid,", n_samples_mcmc, "points in MCMC"))
# Remember y-axis scaling is not meaningful here across methods; I just plotted
# them in a way that they can all be shown in the same graph. 
# The grid method gives just estimates for a limited set of points; the quadratic
# method gives us an entire, smooth posterior; and the MCMC method essentially gives
# us a whole bunch of wiggly lines that circumscribe the 'true' posterior.
# But they all agree pretty well on the most probable proportion of water on the globe.
# Side note: a uniform probability density distribution has the value 1 at all
# x values from 0 to 1, thus summing up to an area of 1.


## VIDEO 2; this is also in chapter 3, but much later ---------------

### SIMULATION AND TESTING (of the model/approach, not the hypothesis) ------------
#### Simulate your data collection procedure to generate results for a given hypothesis --------------
sim_globe <- function(p, N) {
  sample(c("W", "L"), size=N, prob = c(p, 1-p), replace = TRUE)
}
# Define hypothesis
prop_water_on_globe <- 0.7
# Define parameters of the experiment
no_of_replicates <- 10 # of the experiment
no_of_samples # <- 4 # measurements per experiment
# Simulate once (one experiment)
sim_globe(prop_water_on_globe, no_of_samples)
# Replicate, ie. simulate a bunch of such experiments
sim_replicates <- replicate(sim_globe(prop_water_on_globe, no_of_samples), n=no_of_replicates)
# Now we'll see what the outcomes of our replicated, simulated experiments were.
# Here, we just simply count how many waters we got (given an experimental design
# and an assumed given actual proportion of water on the globe).
test_results <- c()
for (i in 1:no_of_replicates) {
  test_results[i] <- (sum(sim_replicates[,i] == "W"))
}
# Plot these counts as histogram:
simplehist(test_results
           , xlab="Simulated water count"
           , ylab="Frequency across replicates"
           )
mtext(paste("simulations given p=", prop_water_on_globe))
# This shows us the distribution of outcomes we could get with this
# sample size. 

#### Testing how the model works given this simulated dataset -----------------
# Re-run the model fitting as if this entire set of simulations were our empirical results
total_n <- no_of_replicates*no_of_samples
total_observations_of_water <- sum(sim_replicates == "W")
# And these are our decisions
grid_size <- 10
n_samples_mcmc <- 1000 # for MCMC
n_mcmc_reps <- 5 # for MCMC

## Graph everything in one graph
grid_approx(W=total_observations_of_water, L=total_n - total_observations_of_water, grid_size)
globe.qa(W=total_observations_of_water, L=total_n - total_observations_of_water) 
globe_mcmc(W=total_observations_of_water, L=total_n - total_observations_of_water, n_samples_mcmc, n_mcmc_reps)
mtext(paste(grid_size, "points in grid,", n_samples_mcmc, "points in MCMC"))

### Benefit of simulation: ---------------------------
# Test whether method can work on known answer
# Explore the effect of sample size and parameter settings
# Confirm that larger sample size leads to more correct answer, and
# small size leads to correct indication of lack of confidence

### Interpreting the entire posterior and not controlling Type I error ------------
# More from video 2:
# McElreath emphasizes also that we should always interpret the entire 
# posterior distribution. One can communicate summaries of the distribution
# (e.g. the mode, or a 50% interval), but the real result is the overall
# distribution. (E.g. he is opposed to focusing on the 95% interval as if it 
# was somehow more informative.)

# "Bayesian statistics does not control Type I error, which is the entire point
# of NHST and 95% intervals." - indeed, it is not trying to give us a unique
# decision or outcome (with known or unknown Type I error) - it is trying to better
# represent our certainty around the outcome. This could be a good thing or a bad 
# thing depending on your goals, and it one reason why it is better to be really 
# clear on what you are trying to achieve with your statistical method. 
# Are you trying to show beyond doubt that some previously unknown or unbelieved
# phenomenon exists? Then perhaps controlling Type I error would be better. 
# Are you trying to estimate as best we can a parameter? Then clearly you want a
# numerical best estimate with an understanding of confidence, and Type I error is
# not particularly relevant nor applicable (nor is Type II error).

### Storytelling and the hypothesis-model-prediction relationship ----------------
## More from chapter 2
# McElreath calls 'storytelling' what is usually referred to as translating a 
# verbal hypothesis into quantitative predictions conditioned on a specific
# method. In any statistical analysis, figuring out what precisely is meant
# /implied by a verbal hypothesis in quantitative terms is important.

# Equally, what is important is to realize that hopefully the results do not depend
# on the details of the 'story' in the sense that the same verbal hypothesis
# can lead to many 'stories'; nonetheless, the justification and specific analysis
# depend on exactly what the 'story' is. And this is indeed a danger in any statistical
# approach. 

# In his words (here the word 'model' is embodying the 'story' of data generation
# primarily, and also the Bayesian inference; and the 'small world' is again
# the world of the data story, whereas the 'large world' is the real world):
# "For now, note that the goal is not to test the truth value of the model’s 
# assumptions. We know the model’s assumptions are never exactly right, in the 
# sense of matching the true data generating process. Therefore there’s no 
# point in checking if the model is true. Failure to conclude that a model is 
# false must be a failure of our imagination, not a success of the model. 
# Moreover, models do not need to be exactly true in order to produce highly 
# precise and useful inferences. All manner of small world assumptions about 
# error distributions and the like can be violated in the large world, but a 
# model may still produce a perfectly useful estimate."

# "Instead, the objective is to check the model’s adequacy for some purpose."
# (!!)


## CHAPTER 3 -------------------------------
# Note: this is still video 2, starting at about the 1 hour mark

### SAMPLING FROM POSTERIOR -----------------------------
# Drawing samples from a calculated distribution, or calculating
# a distribution by drawing samples from some representation (MCMC)

# How to draw samples from a distribution in R. 
# The point of this is to represent/illustrate the posterior distribution by plotting
# or otherwise describing a bunch of points sampled from it (rather than the actual
# continuous distribution, which for complex problems we cannot calculate anyway).

# Here is our data again
observations_of_water
no_of_samples # Remember this is how often we tossed the globe

# First, re-create a posterior distribution (we'll use grid approximation, but
# with a much finer grid to get a smoother posterior)
grid_size <- 1000
W <- observations_of_water
L <- no_of_samples - observations_of_water
posterior <- grid_approx(W, L, grid_size)

# Now we sample from it (a lot)
n_samples <- 10000 # This is how many points we will sample from the posterior
p_grid <- seq(from=0 , to=1 , length.out=grid_size)
# sample() is the actual function sampling from a given distribution, here the 
# posterior distribution. 
samples <- sample(p_grid, prob = posterior, size = n_samples, replace = TRUE)
# or alternatively use beta distribution
#samples <- rbeta(n_samples, W+1, L+1)

# Now we can plot all these sampled points.
plot(samples
     , col = alpha("slateblue1", 0.1)
     , pch = 19
     , ylim = c(0,1)
     , xlab = "Sample number"
     , ylab = "Proportion of water on globe")
# Note that the x-axis is just the row number, and the y-axis is the x-axis value
# from the posterior distribution. 

# We can turn this by 90 degrees again so that we have the hypothesis/parameter/proportion
# of water on the x-axis, and on the y-axis how frequently we got this point in our
# samples from the posterior. 
dens(samples
     , lwd = 4
     , col = "slateblue1"
     , xlab = "Proportion water on globe"
     , adj = 0.1
     , xlim = c(0,1)
     )
# And we can add the beta distribution curve on it (as an approximation of true
# probability)
curve(dbeta(x, W+1, L+1)
      , add=TRUE
      , lty=2
      , lwd=3
      )

### SUMMARIZING a set of samples ----------------------
# (either simulated data or samples from the posterior just to understand the posterior)

# Why do we want samples instead of just the distribution?
# Because by and large, we want to 'interpret' the actual posterior
# distribution, and this is most intuitive by summarizing characteristics
# of a group of samples.

# For example, maybe we want to know how likely it is that the proportion of water
# on the globe is between interval_top and interval_bottom.
# To do this, we just define a function that counts all the sample points in that
# interval, and checks what proportion this is of the total number of samples.
interval_probability <- function(top, bottom, post, resolution, n) {
  p_grid <- seq(from=0 , to=1 , length.out=resolution)
  samples <- sample(p_grid, prob = post, size = n, replace = TRUE)
  in_interval <- subset(samples, samples<top & samples>bottom)
  length(in_interval) / length(samples)
}

interval_top <- 0.75
interval_bottom <- 0.5
interval_probability(interval_top, interval_bottom, posterior, grid_size, n_samples)
# Note: this much simpler version should give the same result...
sum(samples < interval_top & samples > interval_bottom) / n_samples

# Or, we might want to know which values of p (proportion water) are within an interval
# of a certain size, e.g. in the middle 90% interval of probability.
# He calls it 'compatibility interval' (instead of confidence interval)...
compatibility_interval <- function(interval, post, resolution, n) {
  p_grid <- seq(from=0 , to=1 , length.out=resolution)
  samples <- sample(p_grid, prob = post, size = n, replace = TRUE)
  bottom <- (1-interval) / 2
  top <- 1-bottom
  quantile(samples, c(bottom, top))
}

# So note that the result is two values of p (proportion of water), such that 
# 80% of the area of the posterior is between these two values. 
compatibility_interval(0.8, posterior, grid_size, n_samples)

# Should be the same as 
PI(samples, 0.8)

# McElreath actually argues that what we want is not the above, but instead the highest
# probability density interval of that size, i.e. the region where the probability density function 
# is highest over some area:
HPDI(samples, 0.8)

# Other summarizing statistics he suggests are the maximum
p_grid[which.max(posterior)]
# or mean or median
mean(samples)
median(samples)

### LOSS FUNCTIONS ------------------------------------
# We can assume some kind of function for the cost of estimating the actual
# value of p incorrectly. For example, if we say the cost is proportional
# to the distance of our estimate d from p (i.e. abs(d-p)), then
# for a particular true p, we can calculate the cost of every decision d. 
# If we now average this over all possible p but weighted by the posterior
# function, we get an overall loss function of d based on our information
# about p. 
loss <- sapply(p_grid, function(d) sum(posterior*abs(d - p_grid)))
# The minimum of this function is actually the median of the posterior distribution.
p_grid[which.min(loss)]

points(p_grid, loss
       , col = alpha("purple4", 0.1)
       , pch = 19
       , cex = 0.5
       )
abline(v = p_grid[which.min(loss)]
       , lwd = 3
       , lty = 2
       , col = "purple")

# If we defined the cost as (d-p)^2 (quadratic loss), we would get the mean
# of the posterior as the minimum loss. 

# This is important in applied contexts, where the true goal is to make a decision
# based on p rather than actually finding the true p.
# However, one can also frame the decision as 'whether or not to accept the hypothesis',
# and this also (implicitly or explicitly) comes with a loss function for the
# wrong decision. 



## VIDEO 2 (starting at about 1 h) AND CHAPTER 3 -----------------------

# Note that in representing uncertainty about a parameter estimate, we have
# two different sources or types of uncertainty:
# Randomness in choosing paths in the garden of forging data, i.e. sampling 
# (maybe our sample is an odd case for this hypothesis/parameter value)
# VS
# Uncertainty about p (embodied by posterior)

### PREDICTIONS (Posterior predictive distributions) -----------------------
# For example, we might want to know how likely it is that we get 
# 5 water out of 9 samples. What we know is the posterior distribution
# of how likely any proportion of water on the globe is.
# We can sample from the posterior distribution different proportions
# of water, then for each proportion of water calculate a predictive
# distribution for our new outcome; but then we add those across
# different samples from the posterior, to get an overall prediction
# that takes the uncertainty in the posterior into account.

# Recap: note that there are two different types of 'sampling' here:
no_of_samples # how often we toss the globe
samples_from_posterior <- samples # how many samples we draw from the posterior

# Now we're going to generate a posterior predictive distribution. 
PPDist <- c()
for (i in 1:length(samples_from_posterior)) {
  PPDist[i] <- sum(sim_globe(samples_from_posterior[i], no_of_samples) == "W")
}

# What is a PPD? It is a prediction of what outcomes we might get given a particular
# posterior. I.e. it simulates both types of uncertainty, the uncertainty about
# the true value of p (the parameter represented by the posterior) as well as 
# the uncertainty given that even for a given p your experimental measurements 
# (samples) might have different values. 
# Here is an illustration of our PPD:
simplehist(PPDist
           , xlab="Simulated water count"
           , ylab="Frequency across samples from posterior"
)

# McElreath argues here that a creative model checker should come up with other 
# ways to check than the obvious one - e.g. measuring other aspects of the PPDist than
# than which was used to generate the posterior. 
# E.g. in this example, instead of counting the number of water, we could count the
# length of 'runs' of the same result, or the number of switches between W and L, or
# something else. 


### Modeling more complex data collection processes -----------------

## What if there are errors in the data? 
#  We can model these explicitly, and then use that additional information
#  to improve our inference and the estimate of our uncertainty.

# Start with actually modeling the generative process. 
# p is the parameter we want to find out, here the proportion of water
# N is the number of samples we draw
# x1 is the false positive rate, i.e. getting water when it is land
# x2 is the false negative rate, i.e. getting land when it is water
sim_globe2 <- function(p, N, x1, x2) {
  # The actual sample function is the same as before
  truesample <- sample(c("W", "L"), size=N, prob = c(p, 1-p), replace = TRUE)
  # Now we simulate what we observe, which is with a potential error
  obs_sample <- ifelse(truesample=="W",
    ifelse(runif(N)<x1, # this is 'random uniform' distribution of length N
           "L", "W"), # if false negative, return L, otherwise W
    # if the true sample was L to begin with
    ifelse(runif(N)<x2,
           "W", "L")
  )
  return(obs_sample)
}
# McElreath argues that precisely doing this model of the generative process is
# straightforwardly linked to your empirical knowledge. Thus, different from 
# complex 'premade' statistical tests, you can easily incorporate what you know
# about how the data are generated into your statistical model with this method.

## Grid approximation of our new model, i.e. calculating posterior knowing that
# there are error rates
grid_approx_w_error <- function(W, L, grid_size, x) {
  # define grid of possible hypotheses/values of thing to be estimated
  p_grid <- seq(from=0, to=1, length.out=grid_size)
  # define prior across all hypotheses
  prior <- rep(1, grid_size)
  # Standardizing just for the sake of the graph
  prior <- prior/sum(prior)
  # Total samples
  Z <- L+W
  # THIS IS WHERE THE WORK OF THE MODEL HAPPENS:
  # compute likelihood at each value in grid (i.e. data given each hypothesis)
  likelihood <- ((((p_grid*(1-x)) + ((1-p_grid)*x))^W) * ((((1-p_grid)*(1-x)) + (p_grid*x))^L) )/ Z
  # Formula in minute 1:31 in video 2
  
  # compute product of likelihood and prior
  unstd.posterior <- likelihood * prior
  # standardize the posterior, so it sums to 1
  posterior <- unstd.posterior / sum(unstd.posterior)

  return(posterior)
}

# Ok let's do the whole thing. [The owl!]
# (1) simulate data!
# Assume a world to simulate:
prop_water_on_globe #<- 0.7
error_rate_fp <- 0.1
error_rate_fn <- 0.1
# Design experiment:
no_of_replicates <- 10 # of the experiment
no_of_samples <- 5 # measurements per experiment
# Simulate once to test it works
sim_globe2(prop_water_on_globe, no_of_samples, error_rate_fp, error_rate_fn)
# Replicate
(sim_replicates <- replicate(sim_globe2(prop_water_on_globe, no_of_samples, error_rate_fp, error_rate_fn), n=no_of_replicates))
test_results <- c()
for (i in 1:no_of_replicates) {
  test_results[i] <- (sum(sim_replicates[,i] == "W"))
}
# Check that data produced look reasonable
simplehist(test_results
           , xlab="Simulated water count"
           , ylab="Frequency across replicates"
)

# (2) Test (fit model, i.e. calculate posterior) on simulated data
# Re-run everything as if these were our empirical results
total_no_of_samples <- no_of_replicates*no_of_samples
total_observations_of_water <- sum(sim_replicates == "W")
grid_size <- 100
W <- total_observations_of_water
L <- total_no_of_samples - total_observations_of_water

# Or don't do the simulated data but test out specific situations
#W <- 6
#L <- 3

# Now we're doing the actual model, i.e. posterior calculation - we'll do it 
# twice, once with our new model above and once WITH THE MODEL AS IF THERE WAS NO ERROR 

# As if there was no error (this one included a graph so that's what we get)
posterior <- grid_approx(W, L, grid_size)
# Our new model with error
posterior_with_error <- grid_approx_w_error(W, L, grid_size, error_rate_fn)

# And we need this for the graph
p_grid <- seq(from=0 , to=1 , length.out=grid_size)
transparency_of_points <- 15/grid_size + 0.05
y_max <- 1.1 * max(posterior, posterior_with_error)

# (3) Now we plot posterior with both models
plot(p_grid, posterior_with_error 
     , pch = 19
     , type = "b"
     , lwd = 2
     , col = alpha("slateblue1", transparency_of_points)
     , xlab="probability of water" 
     , ylab="posterior probability"
     , ylim=c(0, y_max)
     #, yaxt='n'
)
points(p_grid, posterior 
     , pch = 19
     , col = alpha("red", transparency_of_points)
     , type = "b"
     , lwd = 2
)
# If we want we can add the prior and likelihood from the with-error model.
# I have to repeat the code here (not very elegant) because we encapsulated it in the 
# function and so the variables are not accessible outside of it.
prior <- rep(1, grid_size)
prior <- prior/sum(prior)
Z <- L+W
likelihood <- ((((p_grid*(1-error_rate_fn)) + ((1-p_grid)*error_rate_fn))^W) * ((((1-p_grid)*(1-error_rate_fn)) + (p_grid*error_rate_fn))^L) )/ Z

points(p_grid, prior
       , pch = 19, col = "purple4"
       , type = 'l'
       , lwd = 3)
points(p_grid, likelihood * 0.02/max(likelihood)
       , pch = 19, col = "darkseagreen"
       , type = 'l'
       , lwd = 3)
# As expected, since the prior is flat, the likelihood is very similar to the posterior
# with error (which is the likelihood I used). Also note that I rescaled the likelihood so
# its pattern would be visible. 



### 4 lines of code version --------------------
# of the entire model with grid approx from chapter 4
w <- 6; n <- 9;
p_grid <- seq(from=0,to=1,length.out=100)
posterior <- dbinom(w,n,p_grid)*dunif(p_grid,0,1)
posterior <- posterior/sum(posterior)

plot(posterior)




