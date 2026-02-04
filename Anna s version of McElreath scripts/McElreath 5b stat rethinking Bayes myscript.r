## CODE TO GO WITH BOOKCLUB
## ON Richard McElreath's 'Statistical Rethinking'
## Spring 2025

########### PART Vb - doing the full multiple reg modeling #####################
# LIBRARIES  -----------------------------
library(rethinking)
library(scales)
library(dagitty) # for drawing DAGs
library(viridis)
# Resetting graphics parameters ------------------
par(mar=c(5,5,1,1), oma=c(0,0,0,0), mfrow=c(1,1)) #b,l,t,r
# Defining some color scales
colorscale100 <- viridis(100)
colorscale50 <- magma(50)
taxoncolors <- viridis(5)

# High-level procedure for proper multiple regression with testing and confounds ---------------------------
## (I) Question: define a set of variables and which relationship we are
##           interested in
## (II) Hypotheses: draw possible DAGs representing the different possible
##           relationships, not just the one we are interested in, but 
##           allowing for the 'interesting' one to be present or not
##           THIS IS WHAT HE CALLS THE SCIENTIFIC MODEL
## (III) Going towards 'predictions':
##           The scientific model needs to be translated into a statistical
##           model, that is, into a formula such that the different 
##           possible DAGs are represented by different quantitative
##           parameters, and identifying the one specific parameter related to
##           the key relationship we want to estimate (test). 
## (IV) Causal Inference: From the drawn DAG with the effect we want to estimate, 
##           use the 'backdoor criterion' to figure out which confounds we want 
##           to stratify by (or not). All paths that end in X are 'backdoor' paths;
##           these have to be closed (unless they are closed to begin with).
##           In forks, we DO want to include the factor at the base of the fork;
##           in colliders and pipes, we do NOT want to include it. This 
##           is also a shortcut that prevents us from having to use the 'Full
##           Luxury Bayes' strategy, which is to separately include different 
##           models in a single fitting procedure (for each arrow).
## Actual predictions:
##           What do we expect to happen? We can develop our intuition as
##           well as testing both our model and coding by running simulations. 
##           Or we can shorten the process by thinking about confounds and 
##           which other factors we will need to stratify by. However, since
##           we need to simulate anyway, it seems simulation is always a 
##           good way to go first. 
##           In a sense this is different from NHST 'predictions', since there 
##           is no 'null' and we also don't have a yes/no answer in the 
##           Bayesian modeling strategy as he advocates here. On the other hand,
##           it is similar in that before we run anything on real data, we 
##           figure out whether our model would correctly distinguish between 
##           the hypotheses, and whether it would produce an estimate for our 
##           estimand that is accurate enough for our purposes. 
##           So as in all strong inference, we commit to what we will conclude
##           after defining what method we are using, but before knowing what
##           result this method will yield with the actual data. 
##  (V)  Generative models to achieve simulations but also as a structure for the
##    statistical model we will fit:
##           We code several (!!) generative models, i.e. we use the formula(e)
##           developed above to generate fictive datasets in which either one or
##           another relationship exists, matching the relevant hypotheses/DAGs
##           from above.
##           The statistical model should have the same structure, but different
##           parameters (which, when =0, may mean a particular factor has no
##           effect or a particular relationship does not exist).
##           It makes things more efficient to work with scaled and centered data
##           both for the real data and the generative models.
##   (VI) Testing analysis on generated (fake) data: We define a Bayesian model
##           and priors, very closely matching the structure of the generative 
##           models. This should recover the 
##           parameters as we know them to be from how we generated the data in
##           the first place. If the analysis cannot distinguish between different
##           generative models, then the information we seek may just not be in 
##           the data.
## (VII) Run fitting procedure on actual data. We get fitted parameter values 
##           that tell us which relationships the model 'sees' and how strong they
##           are. 
## (VIII) Illustrate and interpret results
##           Several things should happen here. First, we can look at the parameter
##           estimates (means and sds).
##           Second, by plotting the relationships of the posterior estimates
##           (e.g. the different slope estimates for different parameters, as 
##           in the right-leg, left-leg example), we can see whether there is 
##           structure in the data not illustrated just by the mean and SD of a
##           parameter overall. 
##           Third, by using 'counterfactual plots'. Here we simulate data again,
##           in which we artificially vary a factor as if we were doing an experimental
##           manipulation (breaking any causal associations going into the factor,
##           i.e. 'backdoor' paths). Then we use fitted parameters from the original
##           full model (from real data) to determine the response variable. This
##           is 'testing' whether internally the model assumes a causal association
##           between the manipulated factor and the response. 

# Let's do it! ----------------------------
## (I) The QUESTION: The 'milk' data -----------------------------
## Our goal is to understand what (evolutionarily) drives the calorie content of
## milk in primates. 
data("milk")
milk_comp <- milk[complete.cases(milk),]
d <- data.frame(
  milk_cal = scale(milk_comp$kcal.per.g),
  taxon = as.numeric(milk_comp$clade), # categorical factor
  bodysize = scale(milk_comp$mass),
  rel_brain = scale(milk_comp$neocortex.perc)
)
# Just checking the scaling worked:
par(mar=c(5,5,1,1), oma=c(0,0,0,0), mfrow=c(1,1)) #b,l,t,r
boxplot(d, ylim=c(-3, 3))
points(c(mean(d$milk_cal),mean(d$bodysize),mean(d$rel_brain))~c(1,3,4), col = "darkblue", pch=19)
abline(h=0, lty = 1, lwd = 2)
abline(h=1, lty = 2)
abline(h=-1, lty = 2)

# (II) The HYPOTHESES -----------------------
## Ok now we have a response and three factors, all scaled and centered.
## Now we think about the DAG(s):
par(mfrow=c(2,2))
dag_hyp1 <- dagitty( "dag {
taxon -> rel_brain
taxon -> milk_cal
taxon -> bodysize
rel_brain -> milk_cal
}")
dag_hyp2 <- dagitty( "dag {
taxon -> rel_brain
taxon -> milk_cal
taxon -> bodysize
}")
dag_hyp3 <- dagitty( "dag {
taxon -> bodysize
bodysize -> rel_brain
bodysize -> milk_cal
}")
drawdag(dag_hyp1)
drawdag(dag_hyp2)
drawdag(dag_hyp3)
par(mar=c(5,5,1,1), oma=c(0,0,0,0), mfrow=c(1,1)) #b,l,t,r

# (III) PREDICTIONS/STATISTICAL MODEL ---------------------
# What is our statistical model?
# Other than what is written in the DAG, we assume things like:
# - there are unmeasured factors that affect body size, milk, and brain size
#    independently of taxon: their effects are in the respective error terms
# - we assume each of the three continuous variables has a linear and additive
#    effect on any other factors it has an effect on.
# - otherwise we assume normal distributions: so the effects of parameters are
#   on the means of the effects, and the error terms drive the spread around them
# - we're going to need priors for a bunch of parameters; generally we assume
#    everything to be normally distributed around a mean of 0, and any arrows
#    reflect linear, independent effects.
# - for 'taxon', we need a whole matrix of parameters of how each taxon affects
#   each of the slope sets.

# (IV) CAUSAL INFERENCE -------------------------
## From this we decide what to include in the model. 
## Backdoor paths are the ones that end in X, which in our case is mostly
## relative brain size (the most interesting arrow being the one from brain
## to milk); but at some level we also want to estimate the size of the effect
## of taxon and body size.
## In all hypotheses, taxon is at the base of a fork on relative brain size
## and milk. So we need to include taxon. In hyp1, this complicates things as 
## taxon->brain->milk also form a pipe. In that case the inclusion of taxon 
## reduces the effect of brain that we will find, but there may not be a way around
## that. 

# (V) GENERATIVE MODEL ---------------------

## Parameter assumptions ----------------
### SLOPES ('effects') ------------
# We assume a slope each of the effect of body size and relative brain size on
# milk energy content:
slope_bs <- 2
slope_rb <- 0.5
# We also assume an effect of taxon on milk energy content; this is categorical
# so we have a separate parameter for each taxon. (4 total)
# I addition though, taxon may also affect the other two parameters, body size
# and relative brain size. So since we have another 4 parameters for each of
# those, we get a total of 4 taxa x 3 parameters = 12 parameters. I'll put them
# in a matrix for ease of handling. 
# We will here assume that taxon 1 has no effect on means, taxon 2 increases
# all three means by 1, taxon 3 all three means by 2, etc.
taxon_factors <- matrix(0:3, nrow=3, ncol=4, byrow=TRUE) 
# Lastly, our hypotheses also allow for a possible effect of body size on 
# relative brain size, so we need an effect size for this:
slope_bsrb <- 0.1

### STANDARD DEVIATIONS ('error terms') --------------
E_bs <- 1
E_rb <- 2
E_mc <- 1

### INTERCEPTS ----------------------------
# Arguably with centered data, these intercepts aren't necessarily useful.
# Moreover, one might just see the taxon-specific means as 'intercepts'. 
I_bs <- 0
I_rb <- 0
I_mlk <- 0

## Actual model ---------------------------------
# Here comes the actual model:

# There must be a way to do this with merge() or something, but I'm just
# going to write a lookup function:
taxonmean <- function(tx,par, tx_matrix) {
  if(tx==1) return(tx_matrix[par,1])
  if(tx==2) return(tx_matrix[par,2])
  if(tx==3) return(tx_matrix[par,3])
  if(tx==4) return(tx_matrix[par,4])
}

generative_model <- function(tx_matrix, slope1, slope2, slope12, E1, E2, E3, n) {
  # First we pick a taxon for each of our simulated observations
  taxon <- sample(c(1,2,3,4), n, replace=TRUE)
    # or draw simulated taxa from original categories
    #sample(as.numeric(d$taxon), N, replace = TRUE) 
  
  # Body size is normally distributed around the mean given by taxon (and we
  # simulate N observations).
  bodysize = rnorm(n, I_bs+apply(as.matrix(taxon), 1, function(x) taxonmean(x,1,tx_matrix)), E1)
  
  # Relative brain size is also not affected by anything except perhaps taxon.
  rel_brain = rnorm(n, I_rb+slope12*bodysize+apply(as.matrix(taxon), 1, function(x) taxonmean(x,2,tx_matrix)), E2)
  
  # Milk on the other hand could be affected by taxon, body size, and relative
  # brain size. How much the latter two affect it is determined by the slopes.
  taxonfac <- apply(as.matrix(taxon), 1, function(x) taxonmean(x,3,tx_matrix))
  bsfac <- slope1*bodysize
  rbfac <- slope2*rel_brain
  milk_cal = rnorm(n, I_mlk+taxonfac+rbfac+bsfac, E3)
  
  # Assemble all into a data frame:
  generated_data <- data.frame(
    taxon = taxon,
    bodysize = bodysize,
    rel_brain = rel_brain,
    milk_cal = milk_cal
  )
  return(generated_data)
}

N <- 50 # number of simulated observations
sim_d <- generative_model(taxon_factors, slope_bs, slope_rb, slope_bsrb, E_bs, E_rb, E_mc, N)

### A little plotting to check our simdata -------------------
# Just to see what we wrought, let's plot it:
plot(bodysize~rel_brain
     , data = sim_d
     , col = taxoncolors[taxon+1]  
     , pch = 19
)
# Ok so we clearly have larger body size and relative brain size as a result 
# of taxon.
boxplot(milk_cal~taxon, data=sim_d, col = taxoncolors[2:5])
# And we have an effect of taxon on milk energy content. 
# Now together:
plot(milk_cal~bodysize
     , data = sim_d
     , col = taxoncolors[taxon+1]  
     , pch = 19
)
taxon1 <- subset(sim_d, taxon==1)
taxon2 <- subset(sim_d, taxon==2)
taxon3 <- subset(sim_d, taxon==3)
taxon4 <- subset(sim_d, taxon==4)
abline(lm(milk_cal~bodysize, data=taxon1)
       , lwd = 2
       , col = taxoncolors[2])
abline(lm(milk_cal~bodysize, data=taxon2)
       , lwd = 2
       , col = taxoncolors[3])
abline(lm(milk_cal~bodysize, data=taxon3)
       , lwd = 2
       , col = taxoncolors[4])
if(length(taxon4$milk_cal)>1) 
  abline(lm(milk_cal~bodysize, data=taxon4)
       , lwd = 2
       , col = taxoncolors[5])
# We see that we successfully simulated milk contents that depend both on body
# size and taxon. 
# Remember the specific slopes etc. that we simulated are arbitrary; we just want
# to see that our statistical procedure can recover them. 
# Same for brains:
plot(milk_cal~rel_brain
     , data = sim_d
     , col = taxoncolors[taxon+1]  
     , pch = 19
)
abline(lm(milk_cal~rel_brain, data=taxon1)
       , lwd = 2
       , col = taxoncolors[2])
abline(lm(milk_cal~rel_brain, data=taxon2)
       , lwd = 2
       , col = taxoncolors[3])
abline(lm(milk_cal~rel_brain, data=taxon3)
       , lwd = 2
       , col = taxoncolors[4])
if(length(taxon4$milk_cal)>1) 
  abline(lm(milk_cal~rel_brain, data=taxon4)
       , lwd = 2
       , col = taxoncolors[5])
# Ok, same, and given the parameter values we chose a little more spread. 

# (VI) BAYESIAN MODEL AND PRIORS ----------------------
# This 'full luxury Bayes' model includes all four variables and implicitly
# all the possible arrows from our DAGs:
# taxon affects all three other variables, and both body size and relative brain
# size may affect milk energy content. 
# Each of these also have intercepts and error terms.
# The priors are included in this model code, one could also define them 
# first, but there are quite a number of parameters.
# I left out intercepts, i.e. assuming they are = 0
full_model <- function (dat) {
  quap(
    alist(
      bodysize ~ dnorm(mu_bs, sigma_bs),
      rel_brain ~ dnorm(mu_rb, sigma_rb),
      milk_cal ~ dnorm(mu_mlk, sigma_mlk),
      mu_bs <- tx_bs[taxon],
      mu_rb <- tx_rb[taxon] + fit_S_bsrb*bodysize,
      mu_mlk <- tx_mlk[taxon] + fit_S_bs*bodysize + fit_S_rb*rel_brain,
      tx_bs[taxon] ~ dnorm(0,2),
      tx_rb[taxon] ~ dnorm(0,2),
      tx_mlk[taxon] ~ dnorm(0,2),
      fit_S_bsrb ~ dnorm(0,1), # slope rb~bs
      fit_S_bs ~ dnorm(0, 1), # slope  mlk~bs
      fit_S_rb ~ dnorm(0, 1), # slope  mlk~rb
      sigma_bs ~ dexp(1), # error term bs
      sigma_mlk ~ dexp(1), # error term mlk
      sigma_rb ~ dexp(1)   # error term rb
    )
    , data = dat
  )
}

## TEST ON SIMULATED DATA ------------------
result_sim <- full_model(sim_d)
precis(result_sim, depth=2)
## So far, we have 17 parameters, and we fit them using a large number of data points.
## Even with 1000 or 2000 points, I had runs where the model didn't converge, 
## and when it doesn't converge, some of the parameters can be pretty far off, 
## or even undetermined, although others may be pretty close (we know what the 
## correct values are because this is the simulated dataset where we simply 
## assumed some relationships).
## Let's see how close we are:
par(mfrow=c(1,1), mar=c(5,8,1,1), oma=c(0,0,0,0))
plot(coeftab(result_sim)
     , prob=0.95 # Length of interval shown, relative to posterior distribution
     , col.ci = "darkblue" # Color of the interval line
     , cex=0.5 # Size of the parameter labels 
     , pt.cex = 1 # Size of opints
     , pch=19 # Type of points
     , col="purple" # Color of the name of the model fit as well as the points
     , xlab="Posterior distribution of parameter value"
)
# These are the 'correct' values we put into the simulation:
points(c(16, 13,10,7,1,4)~c(slope_bsrb, slope_bs,slope_rb,E_bs,E_rb,E_mc), col="darkred", pch=1, cex=2)
for (i in 1:4) {
  points(c((19-i)*3-2, (15-i)*3-2, (11-i)*3-2) ~ taxon_factors[,i]
         , col="darkred", pch=1, cex=2)
}

## What this tells us is that we may or may not be able to estimating this many 
## parameters with the real dataset. We can try; or we can try to get additional 
## information: better priors perhaps, or exclude some parameters (perhaps we 
## can independently measure or estimate their values).

## REEXAMINING PRIORS -----------------------
# Let's plot the real&simulated data
plot(milk_cal~bodysize
       , data = sim_d
       , col = alpha(taxoncolors[taxon+1], 1)
       , pch = 19
       )
points(milk_cal~bodysize
       , data = d
       , col = "black"
       , pch = 19)
# It is not required that the simulated and real data are at all similar, except
# in that we want the same priors to work for both. The real data are n=17, while
# you can simulate any number of data points you want. 

# The function extract.prior() automatically extracts combinations of the fitted parameters 
# drawn according to the prior (i.e. prior not posterior values). The default
# is 1000 of them. 
prior <- extract.prior(result_sim)
# This defines two body size values we'll use as the end points of our lines
# (remember these are in standardized, i.e. z-score, units):
bs_seq <- c(-6,6)
# The link() function from rethinking is explained in the overthinking box on
# page 110 (end of text section 4.4.3.4 Plotting regression intervals);
# it does what code section 4.58 does:
# link() takes the 'data' given and plops them into the function (model)
# used to calculate the mean response in the model (here: 'result'), assuming the 
# parameter values given in 'post' (here actually the priors not posteriors).
prior_mus <- link(result_sim, post=prior, data=list(bodysize=bs_seq, taxon=c(1,1), rel_brain=c(0,0)))
# So the error terms are ignored here. Also, since in the prior all taxa are 
# identical, I am just feeding in one taxon (the priors for the others will be the same);
# and I am assuming rel_brain = 0 (which should be the mean) so I can just plot
# the priors for the milk_cal~bodysize relationship.
for (i in 1:50) 
  lines(bs_seq, prior_mus$mu_mlk[i,] 
        , col = alpha("black",0.3) 
  )
# Ok, relative to the real data lines are all over the place, but relative to 
# the simdata we see that priors barely reach the actual slopes. 
# So not a lot of room to improve priors.  
# Priors as plotted here also don't cover the higher intercepts of some taxa.

# One of the sources of high variation here is taxon: 
# since we're not assuming anything about
# differences between taxa on any dimension, a lot of variance is introduced this
# way. No clear way around it if figuring out the uncertainty around taxonomic
# differences is a goal here. 
# But if it is not, i.e. we already know and accept that taxa differ in average
# body size for example, we may not have to include that as a to-be-fitted parameter
# in the model. But the alternatives all depend on specific assumptions.

# It might be time to just try on real data...

# (VII) TEST ON REAL DATA --------------------------------
result <- full_model(d)
precis(result, depth=2)
# Ok, in my version this model did in fact converge, so we have 17 fitted parameters
# despite only having 17 data points.
# McElreath talks repeatedly about sample size: importantly, in the Bayesian 
# method, we calculate just what our best guess is about the true parameter values,
# which is valid at any sample size. For small sample sizes or uninformative data, 
# our posterior parameter estimates will just be pretty close to the priors. 
# See e.g. 'Rethinking: sample size and reliable inference' on p30 (section 2.2.2).

# (VIII) ILLUSTRATE & INTERPRET -------------------------------
# McElreath's entire approach is about not just taking a single number as a result,
# but to try and examine it from different angles to develop an understanding
# of what it means.
#### Understanding posterior of parameters -------------------
# Our first job is to understand the estimated parameter means, the output from 
# precis(result, depth=2).
# A plot of the 'confidence intervals' for each parameter is common:
par(mfrow=c(1,1), mar=c(5,8,1,1), oma=c(0,0,0,0))
plot(coeftab(result)
     , prob=0.95 # Length of interval shown, relative to posterior distribution
     , col.ci = "darkblue" # Color of the interval line
     , cex=0.5 # Size of the parameter labels 
     , pt.cex = 1 # Size of opints
     , pch=19 # Type of points
     , col="purple" # Color of the name of the model fit as well as the points
     , xlab="Posterior distribution of parameter value"
     )
# In this case, this shows a high intercept for body size for taxon 1, a low
# intercept for brain size for taxa 2 and 4, but all taxa overlap in their 
# milk-energy intercept intervals. 
# Also, both the effect of body size and of relative brain size on milk overlap
# zero.

# But this is just plotting the means and distribution for each parameter averaged
# (or rather added) over the distributions of all the other parameters. However,
# parameter estimates may well depend on each other. To examine this, we extract
# a large number of samples (combinations of parameter values) from the posterior:
no_samples <- 10000
post <- extract.samples(result, no_samples)

# Calculating some overview outcomes for the taxon effects
for (i in 1:no_samples) {
  post$tx_bs_var[i] <- sd(post$tx_bs[i,])
  post$tx_bs_mu[i] <- mean(post$tx_bs[i,])
  post$tx_rb_var[i] <- sd(post$tx_rb[i,])
  post$tx_rb_mu[i] <- mean(post$tx_rb[i,])
  post$tx_mlk_var[i] <- sd(post$tx_mlk[i,])
  post$tx_mlk_mu[i] <- mean(post$tx_mlk[i,])
}
post$rel_bs_rb_var <- post$tx_bs_var / post$tx_rb_var

par(mfrow=c(1,1), mar=c(5,5,1,1), oma=c(0,0,0,0))
plot(fit_S_bs ~ fit_S_rb
     , data = post
     , col = colorscale50[as.numeric(cut(post$tx_mlk_mu,breaks = 50))]
     , pch = 19
)
abline(lm(fit_S_bs ~ fit_S_rb, data=post), lty=2,lwd=2,col=rangi2)
# It is clear that the estimated effect of body size
# and the estimated effect of relative brain size negatively correlate (i.e. 
# if one is estimated larger, the other is estimated smaller). 

# However, it is less clear how the taxon effects are linked to these slopes.
# We can experiment with a few other relationships:

plot(tx_mlk_var ~ tx_mlk_mu
     , data = post
     , col = colorscale50[as.numeric(cut(post$tx_mlk_mu,breaks = 50))]
     , pch = 19
     )
abline(lm(tx_mlk_var ~ tx_mlk_mu, data=post), lty=2,lwd=2,col=rangi2)
plot(tx_bs_var ~ fit_S_bs
     , data = post
     , col = colorscale50[as.numeric(cut(post$tx_bs_mu,breaks = 50))]
     , pch = 19
)
abline(lm(tx_bs_var ~ fit_S_bs, data=post), lty=2,lwd=2,col=rangi2)
plot(tx_bs_mu ~ fit_S_bs
     , data = post
     , col = colorscale50[as.numeric(cut(post$tx_bs_var,breaks = 50))]
     , pch = 19
)
abline(lm(tx_bs_mu ~ fit_S_bs, data=post), lty=2,lwd=2,col=rangi2)
plot(tx_bs_var ~ fit_S_rb
     , data = post
     , col = colorscale50[as.numeric(cut(post$tx_rb_var,breaks = 50))]
     , pch = 19
)
abline(lm(tx_bs_var ~ fit_S_rb, data=post), lty=2,lwd=2,col=rangi2)

#### Comparing predictions to original data  ----------------------
# This is pretty much exactly his code, from fig. 5.5.
mu_ <- link(result)
pred_obs_ <- sim(result)
# Since no data are given, this is for 17 observations as in the original data
mu_mean <- apply(mu_$mu_mlk, 2, mean)
mu_PI <- apply(mu_$mu_mlk, 2, PI)
pred_obs_PI <- apply(pred_obs_, 2, PI)

# plot
par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(0,0,0,0)) #b,l,t,r
scale_SD <- 2.5
transparency <- 0.1
plot(mu_mean ~ d$milk_cal
     , xlab = "Observed milk energy content"
     , ylab = "Predicted milk energy content"
     , ylim = c(-scale_SD, scale_SD)
     , xlim = c(-scale_SD, scale_SD)
)
# A perfect fit would look like this:
abline(a=0, b=1, lty=2)
# Code for plotting error bars on points!
for (i in 1:nrow(d)) 
  lines(rep(d$milk_cal[i],2), mu_PI[,i], col=rangi2)
# What we see is that (a) the fit isn't terribly good and 
# (b) the confidence intervals on data points are large.
# Perhaps also (c) the 'actual' pattern may be a bit decelerating 
# (if there is a relationship at all).


#### Posterior prediction plots with more information --------------
# But McElreath calls this 'golem-speak': instead, we want to visualize what this
# means in terms of points on a real data plot.

# Remember there are two types of uncertainty: inferential uncertainty and outcome
# variability (check out https://www.pnas.org/doi/10.1073/pnas.2302491120).
# McElreath calls them both 'uncertainty', but the first is uncertainty about parameter
# values, i.e. the spread in the posterior, our uncertainty about whether relationships
# exist and how strong they are. The second is the noise in the actual individual 
# measurements, described by the various 'sigma' parameters. I.e. even if we 
# know exactly the effects and effect sizes, we cannot predict the precise value
# of each individual measurement: the data spread around a predicted mean. 

# When we visualize a posterior prediction, we can visualize the lines that are
# the predicted mean outcomes for a given factor level (the 'regression line'),
# or we can visualize the entire distribution of likely outcomes, which includes the 
# spread of the data points around these lines. We can do that either as some kind of
# shading for a defined interval (e.g. shading the interval covering 95% of data points)
# or we can simply sample a bunch of points from the predicted distribution and plot
# those. 
# Note that even the 'regression line' itself is not a single line in a Bayesian frame-
# work, since an uncertainty is attached to the slope (and possibly the intercept as well).
# So we might visualize a bunch of lines illustrating the uncertainty, or the distribution
# of slope estimates. 

# Number of samples
n_samples <- 10000
# Now we extract samples from the posterior (i.e. parameter estimates, 
# reflecting inferential uncertainty):
PostParameterSamples <- extract.samples(result, n=n_samples) # 10000 rows (samples) x 17 parameters

# With the parameter samples in hand, we can simulate data points according to the
# model, i.e. posterior prediction points (i.e. reflecting both types of uncertainty).
# Since we already wrote it, we can simply use our generative model (which has the
# same structure as the statistical model we fit):
PostPredPoints <- matrix(nrow = 0, ncol = 4)
colnames(PostPredPoints) <- c("taxon", "bodysize", "rel_brain", "milk_cal")
for (i in 1:n_samples) {
  s1 <- PostParameterSamples$fit_S_bs[i]
  s2 <- PostParameterSamples$fit_S_rb[i]
  s12 <- PostParameterSamples$fit_S_bsrb[i]
  e1 <- PostParameterSamples$sigma_bs[i]
  e2 <- PostParameterSamples$sigma_rb[i]
  e3 <- PostParameterSamples$sigma_mlk[i]
  tx_m <- rbind(PostParameterSamples$tx_bs[i,], PostParameterSamples$tx_rb[i,], PostParameterSamples$tx_mlk[i,])
  PostPredPoints <- rbind(PostPredPoints, generative_model(tx_m,s1,s2,s12,e1,e2,e3,1))
}

# We would like to also show just the inferential uncertainty, i.e. we want to 
# calculate for each taxon the mlk~rb relationship, i.e. the mean predicted mlk
# value for each taxon by rb combination, without the 'sigma' error. 
# We'll do it the same way:
PostMeans <- matrix(nrow = 0, ncol = 4)
colnames(PostMeans) <- c("taxon", "bodysize", "rel_brain", "milk_cal")
for (i in 1:n_samples) {
  s1 <- PostParameterSamples$fit_S_bs[i]
  s2 <- PostParameterSamples$fit_S_rb[i]
  s12 <- PostParameterSamples$fit_S_bsrb[i]
  e1 <- 0
  e2 <- 0
  e3 <- 0
  tx_m <- rbind(PostParameterSamples$tx_bs[i,], PostParameterSamples$tx_rb[i,], PostParameterSamples$tx_mlk[i,])
  PostMeans <- rbind(PostMeans, generative_model(tx_m,s1,s2,s12,e1,e2,e3,1))
}
# We can't just prescribe an array of rel_brain over which to calculate these
# means, because rel_brain itself may be a result of taxon, i.e. 'confounded'. 
# If we want to preserve that effect, we can't independently determine it (although
# see later for counterfactual plots).

# Ok so now we have posterior prediction points and posterior means.
# We can plot these directly as points or summarize them by to calculate an 
# area of shading.
# Here, I'm going to stick with the points themselves. An area of shading makes
# more sense when we actually have an independent variable on the x-axis, so 
# that we can calculate the expected mean for each x value. If you want a smoother
# graph, you can increase the number of samples (in line 559) and decrease the 
# transparency when plotting the points below. 

# Add the priors as well if you'd like: 
# We'll extract priors exactly as above but this time with relative brain size
# on the x-axis. 
rb_range <- c(-6,6)
priors <- extract.prior(result)
prior_mus <- link(result, post=priors, data=list(bodysize=c(0,0), taxon=c(1,1), rel_brain=rb_range))
# Note that as above, the second kind of uncertainty, the sigma, are ignored here
# (although they are part of the prior).
# Also note that I only extracted priors for taxon 1: but since the priors are
# identical across taxa, this is irrelevant. I'm also just showing the priors 
# as they would be for bodysize=0, i.e. the mean of the original data.

## GRAPH ---------------------------------
## Actual plot setup is best in reverse order, to put important info on top
# frame ----------------
par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(0,0,0,0)) #b,l,t,r
y_max <- max(d$milk_cal)
scale_SD <- 5
transparency <- 0.1
plot(NULL
     , xlab = "Relative brain size"
     , ylab = "Milk energy content"
     , ylim = c(-scale_SD, scale_SD)
     , xlim = c(-scale_SD, scale_SD)
)

# plot priors -------------------------
effectline <- function(rb, s, i){
  mlk<- rb*s + i
  return(mlk)
}
for (i in 1:(round(n_samples/10))) {
  lines(rb_range, prior_mus$mu_mlk[i,] 
        , col = alpha("purple",0.2) 
        , lty = 2
  )
  curve(effectline(x,priors$fit_S_rb[i],0)
        , from=rb_range[1]
        , to=rb_range[2]
        , add=TRUE
        , lwd=2
        , lty=2
        , col=alpha("purple", 0.3)
  )
}
# plot uncertainty just in rb slope --------------------
# What we really said we had set out to know is whether there is
# an effect of relative brain size directly on milk energy content,
# i.e. beyond the effect driven by taxon.
# The parameter that reflects this is the rb slope. 
# We did not model any statistical interactions, so there is a single
# (no taxon specific) estimate for this slope.
for (i in 1:(round(n_samples/3))) {
  curve(effectline(x,PostParameterSamples$fit_S_rb[i],0)
        , from=rb_range[1]
        , to=rb_range[2]
        , add=TRUE
        , lwd=3
        , col=alpha("black", transparency*0.2)
  )
}

# plot ppp ------------------------------
points(PostPredPoints$milk_cal~PostPredPoints$rel_brain
       , col=alpha(taxoncolors[PostPredPoints$taxon+1], transparency)
       , pch=19
)


# plot posterior means --------------------
# These lines are regular regression lines on the points explained
# below, meaning they reflect our mean estimate of the combination
# of taxon and rb effects. 
for (i in 1:4) {
  taxonline <- subset(PostMeans, PostMeans$taxon==i)
  abline(lm(taxonline$milk_cal~taxonline$rel_brain)
         , col=alpha(taxoncolors[i+1], 1)
         , lwd=3
  )
}
# The following points are the actual outcomes *if there were no 
# error* around means at all. So the scatter reflects just the uncertainty
# in the taxon intercepts (on both mlk and rb) and the slope of mlk~rb. 
points(PostMeans$milk_cal~PostMeans$rel_brain
       , bg=alpha(taxoncolors[PostMeans$taxon+1], transparency*1.5)
       , col=alpha("black", transparency*2)
       , pch=21
)

# plot original data ------------------
points(milk_cal~rel_brain
  , data=d
  , bg=alpha(taxoncolors[taxon+1], 1)
  , col="black"
  , lwd=3
  , pch=23
  , cex=1.5
)

## end GRAPH --------------------

# Ok so now we understand more intuitively something about what 
# our model extracted from the data, namely that it is not unlikely that there is some
# positive effect of brain size on milk energy content, but we have considerable 
# uncertainty around this.

#### Comparing predictions to priors, to see what we learned from data ----------------------
# Another approach McElreath suggests is to compare the preditions to the priors,
# to see whether the model actually learned (and how much) from the data. 
# So far he has not really given a particularly systematic way of looking at this,
# but even here in our graph we can see that the prior was pretty unconstrained, 
# allowing for a large range of slopes in this relationship - while the fitted
# slopes from the model, i.e. the posterior for the rb slope parameter, is not 
# as large a range. 

# Another way, perhaps, of showing this is this:
dataslope <- lm(PostMeans$milk_cal~PostMeans$rel_brain)$coefficients[2]
abline(h=0, lty=2, col="grey")
boxplot(priors$fit_S_rb, PostParameterSamples$fit_S_rb, dataslope
        , names = c("Prior", "Posterior", "Posterior means")
        , range = 0
        , ylab = "Estimated effect of rb on mlk"
        )
# The range of the estimated slope in the posterior is clearly smaller than that
# of the prior, even if still considerable. 
# 'Posterior means' here means the simple regression slope of the (no-sigma) 
# means from the posterior samples - i.e. this includes the possibly confounding
# effects of taxon and body size on the relationship. 

#### Looking at 'counterfactual' plots -----------------
n_resolution <- 100
# Generating a counterfactual, i.e. 'manipulated' set of rb data
rb_seq <- seq(from=-4, to=4, length.out=n_resolution)
# Simulating the effect first on bodysize, then mlk_cal:
# the vars argument in sim() tells the model what order to simulate the variables in
counter_input_dat <- data.frame(rel_brain=rep(rb_seq,4)
                                , taxon=c(rep(1,n_resolution), rep(2,n_resolution), rep(3,n_resolution),rep(4,n_resolution)))
counter_pred <- sim(result, data=counter_input_dat, vars=c("bodysize","milk_cal"))
# Display counterfactual predictions
par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(0,0,0,0)) #b,l,t,r
scale_SD <- 5
transparency <- 0.2
plot(colMeans(counter_pred$milk_cal) ~ rep(rb_seq,4)
     , ylim=c(-scale_SD,scale_SD) 
     , cex= 0.5
     , pch=19
     , col=taxoncolors[1+c(rep(1,n_resolution), rep(2,n_resolution), rep(3,n_resolution),rep(4,n_resolution))]
     , xlab="manipulated rb"
     , ylab="counterfactual mlk"
     )
counter_PI <- apply(counter_pred$milk_cal,2,PI)
for (i in 1:4){
  shade(counter_PI[,((i-1)*n_resolution+1):(i*n_resolution)]
    , rb_seq
    , col=alpha(taxoncolors[1+i], transparency)
  )
}
# Ok, cool - now we have independently manipulated taxon and relative brain
# size, and allowed taxon to influence body size and body size to influence 
# mlk, as well as having brain size affect mlk - and according to the posterior,
# brain size still has an effect even within each taxon. 

# Now we can do the symmetrical thing with body size:
counter_input_dat <- data.frame(bodysize=rep(rb_seq,4)
                                , taxon=c(rep(1,n_resolution), rep(2,n_resolution), rep(3,n_resolution),rep(4,n_resolution)))
counter_pred <- sim(result, data=counter_input_dat, vars=c("rel_brain","milk_cal"))
plot(colMeans(counter_pred$milk_cal) ~ rep(rb_seq,4)
     , ylim=c(-scale_SD,scale_SD) 
     , cex= 0.5
     , pch=19
     , col=taxoncolors[1+c(rep(1,n_resolution), rep(2,n_resolution), rep(3,n_resolution),rep(4,n_resolution))]
     , xlab="manipulated bs"
     , ylab="counterfactual mlk"
)
counter_PI <- apply(counter_pred$milk_cal,2,PI)
for (i in 1:4){
  shade(counter_PI[,((i-1)*n_resolution+1):(i*n_resolution)]
        , rb_seq
        , col=alpha(taxoncolors[1+i], transparency)
  )
}
# Unlike in the book, I did this completely symmetrically, i.e. simulating both
# body size and brain size in both graphs; taxon can affect both and can affect 
# mlk directly as well. Here, the variable that is being manipulated is not affected
# by taxon in these graphs, so the effect we see is from the arrow from that variable
# to mlk, not from the fork of taxon affecting both the variable and mlk. 
# As in the book, any effect of body size on rb is also controlled for in the left 
# graph. 

#### Plotting/analyzing residuals *of factors against each other* -----------------

# Now, to understand/illustrate this, we might also want to plot 
# effects against residuals of A or M. By that we mean we first
# regress M on A, then calculate the difference between predicted and actual M.
dat <- d
rb_on_bs_for_residuals <- quap(
  alist(
    rel_brain ~ dnorm(mu_rb_bs, sigma_rb_bs),
    mu_rb_bs <- I_rb_bs + slope_rb_bs * bodysize,
    I_rb_bs ~ dnorm(0, 0.2) ,
    slope_rb_bs ~ dnorm(0, 0.5),
    sigma_rb_bs ~ dexp(1)
  ) , data = dat)
mu <- link(rb_on_bs_for_residuals)
mu_mean <- apply(mu, 2, mean)
mu_resid_rb_bs <- dat$rel_brain - mu_mean

bs_on_rb_for_residuals <- quap(
  alist(
    bodysize ~ dnorm(mu_bs_rb, sigma_bs_rb),
    mu_bs_rb <- I_bs_rb + slope_bs_rb * rel_brain,
    I_bs_rb ~ dnorm(0, 0.2) ,
    slope_bs_rb ~ dnorm(0, 0.5),
    sigma_bs_rb ~ dexp(1)
  ) , data = dat)
mu <- link(bs_on_rb_for_residuals)
mu_mean <- apply(mu, 2, mean)
mu_resid_bs_rb <- dat$bodysize - mu_mean


# Plotting the residuals against the response: 
par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(0,0,0,0)) #b,l,t,r
scale_SD <- 2
transparency <- 0.2
plot(milk_cal~mu_resid_rb_bs
     , data = dat
     , xlim=c(-scale_SD,scale_SD) 
     , ylab = "Milk energy content"
     , xlab = "Relative brain size residuals on body size"
     , col = rangi2
     , pch = 19
)
abline(lm(dat$milk_cal~mu_resid_rb_bs))
text(mu_resid_rb_bs+0.1, dat$milk_cal
     , labels = milk_comp$species
     , cex=0.6, col="black"
     , pos = 4)
# What this shows is essentially the effect of relative brain size on milk
# energy content *on top of* any effects of body size. It's another way to 
# illustrate the information contained in the model about what is the causal 
# effect, and what is an effect resulting from a 'pipe' or 'fork' (assuming
# the model has detected this).
# Here, we determined the residuals just from the raw data, ignoring taxon:
# so any shared effects of taxon on both rb and bs are also controlled for. 

plot(milk_cal~mu_resid_bs_rb
     , data = dat
     , xlim=c(-scale_SD,scale_SD) 
     , ylab = "Milk energy content"
     , xlab = "Body size residuals on relative brain size"
     , col = rangi2
     , pch = 19
)
abline(lm(dat$milk_cal~mu_resid_bs_rb))
text(mu_resid_bs_rb+0.1, dat$milk_cal
     , labels = milk_comp$species
     , cex=0.6, col="black"
     , pos = 4)
# Unlike in fig. 5.4, we do actually find that both variables have an effect on 
# the response, mlk, even when we control for the other respective factor. 

## --------------------------------


