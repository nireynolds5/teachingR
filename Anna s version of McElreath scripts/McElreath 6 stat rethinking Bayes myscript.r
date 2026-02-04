## CODE TO GO WITH BOOKCLUB
## ON Richard McElreath's 'Statistical Rethinking'

########### PART VI #####################

library(rethinking)
library(scales)
library(dagitty) # for drawing DAGs
opar <- par()
## Chapter 5 (still)
## VIDEO 4 (still)

data(Howell1)
d <- Howell1[ Howell1$age >= 18 , ]
str(d)

#  CATEGORICAL VARIABLES

# Interesting point made here. If the categorical variable
# is treated like any other factor, i.e. as a factor multiplied
# by a slope parameter, this generates more uncertainty
# for whichever actual category state is '1' instead of '0',
# since the '1' state has two parameters (intercept and slope),
# whereas the '0' only has the intercept. 

# Also unclear how this would work with categories that have 
# more than two states. 
# PS: Actually you have to turn your categories into many categories, each
# with its own yes/no parameters! This is what a regular linear regression 
# does, which is why there are these incomprehensible category comparison
# outputs and why one category is always designated the comparison category!

# Instead, one can define a (numerical) index for each state:
d$sex <- ifelse( d$male==1 , 2 , 1 )
# We now turned the categorical factor back into a number.

m5.8 <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a[sex] ,
    a[sex] ~ dnorm( 178 , 20 ) ,
    sigma ~ dunif( 0 , 50 )
    )
  , data=d
  )
# The number was used as an index, to define [i] different
# parameters for the intercept.

precis( m5.8 , depth=2 )

post <- extract.samples(m5.8)
post$diff_fm <- post$a[,1] - post$a[,2]
precis(post , depth=2)
# Posterior distribution of the difference! (in height, between
# the category states, i.e. male vs female).
dens(post$a[,1], xlim=c(min(d$height), max(d$height)), lwd=3, col=2, xlab="height")
dens(post$a[,2], lwd=3, col=4, add=TRUE)

# In addition to this, we might want to see not just how different
# the means are but how different individuals are likely to be.

# Simulating individuals:
H1 <- rnorm(1000, post$a[,1], post$sigma)
H2 <- rnorm(1000, post$a[,2], post$sigma)
dens(H1, xlim=c(min(d$height), max(d$height)), lwd=3, col=2, xlab="height")
dens(H2, lwd=3, col=4, add=TRUE)

H_contrast <- H2-H1
dens(H_contrast, lwd=3, col=1, xlab="posterior height contrast")
# This shows the distribution of differences between any two individuals
# from the two groups. Much different from means!
# We can also calculate the proportion of cases where the woman was taller:
sum(H_contrast<0)/1000
# or the man:
sum(H_contrast>0)/1000

# Quite striking!

##

data(milk)
d <- milk
d$K <- scale( d$kcal.per.g )
unique(d$clade)
# Forcing the categorical variable into a numerical one:
d$clade_id <- as.integer( d$clade )

# Does this also work with the regular linear regression 
# approach?
reg_lm <- lm(K~clade, data=d)
summary(reg_lm)
# Hm. Not sure how I can force it to do what I want.

m5.9 <- quap(
  alist(
    K ~ dnorm( mu , sigma ),
    mu <- a[clade_id],
    a[clade_id] ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ) , data=d )
labels <- paste( "a[" , 1:4 , "]:" , levels(d$clade) , sep="" )
plot( precis( m5.9 , depth=2 , pars="a" ) , labels=labels ,
      xlab="expected kcal (std)" )
# Ha! Cool!

# We're making up another category called 'house':
set.seed(63)
d$house <- sample( rep(1:4,each=8) , size=nrow(d) )

m5.10 <- quap(
  alist(
    K ~ dnorm( mu , sigma ),
    mu <- a[clade_id] + h[house],
    a[clade_id] ~ dnorm( 0 , 0.5 ),
    h[house] ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ) , data=d )
labels <- paste( "h[" , 1:4 , "]:" , c("Gryffindor", "Hufflepuff", "Ravenclaw", "Slytherin"), sep="" )
plot(precis( m5.10, depth=2 , pars="h") , labels=labels ,
      xlab="expected kcal (std)" )
precis(m5.10, depth=2)

# Quote:
# Lurking underneath this example is a more fundamental mistake
# in interpreting statistical significance: The mistake of accepting the null hypothesis. Whenever
# an article or book says something like “we found no difference” or “no effect,” this usually means
# that some parameter was not significantly different from zero, and so the authors adopted zero as the
# estimate. This is both illogical and extremely common.

# My comment:
# It is important to realize that 'accepting the null hypothesis' is not 
# concluding the estimate is zero. 
# There are two separate issues: Type II errors can be high - so failure to reject
# the null hypothesis is not strong evidence against the alternative. This should be 
# reflected in language, and possibly accompanied by attempts to quantify uncertainty,
# eg with power analysis.
# Second, what does it mean to accept the null? It means we believe that 
# the evidence from the current study wasn't enough to override our prior, not
# perhaps quantified but resulting from current state of the world, that
# the null hypothesis is correct. It is not an attempt to quantify a parameter,
# it is an attempt to make a qualitative decision about a scientific hypothesis
# possibly tenuously related to the prediction of a non-zero parameter estimate.

# In any case, reporting and interpreting *effect size estimates* 
# AND *uncertainty in effect size estimates* is good and important.


## VIDEO 4

# more on categories: different intercept AND slope
# So a different intercept essentially means there is a direct
# effect of S on W, regardless of H; 
# a different slope means H influences W differently across different S

data(Howell1)
realdat <- Howell1[ Howell1$age >= 18 , ]
realdat$S <- ifelse(realdat$male==1 , 2 , 1)
realdat$H <- scale(realdat$height)
realdat$W <- scale(realdat$weight)

# Simulating the weight-height relationship for 2 sexes
sim_HW <- function(S, b, a) {
  N <- length(S)
  H <- ifelse(S==1, 150, 160) + rnorm(N, 0, 5)
  W <- a[S] + b[S]*H + rnorm(N,0,5)
  data.frame(S,H,W)
}
S <- rbern(100) + 1 # make a hundred random males or females
simdat <- sim_HW(S, b=c(0.5, 0.6), a=c(0,10))
simdat$H <- scale(simdat$H)
simdat$W <- scale(simdat$W)
# Ok this is awkward because after scaling the data, I can't extract
# the parameter from the generating function as easily. 
# Should have centered the model anyway? Or written the sim function
# for centered data?

d <- simdat
d <- realdat

# Quick note: I am not centering on H in the model (i.e. using H 
# not H-H_avg in the model) - only because I am assuming data
# are already centered. 

m_HWS <- quap(
  alist(
    W ~ dnorm(mu, sigma),
    mu <- a[S] + b[S]*H,
    a[S] ~ dnorm(0, 0.5),
    b[S] ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ) 
  , data=d)
precis(m_HWS, depth=2)
labels <- c("a[1]", "a[2]", "b[1]", "b[2]")
plot(precis(m_HWS, depth=2, pars=c("a","b"))
     , labels=labels
     , xlab="Parameter estimate [std]")
# Separate estimates for intercept and slope in each category!

# But as always, we want to illustrate this not only in golem speak.
post <- extract.samples(m_HWS, n=1000)
# One thing we might want is the actual difference between the categories -

# CONTRASTS
# You want the mean difference not the difference of means!! Not the same!

post$diff_Intercept <- post$a[,1] - post$a[,2]
post$diff_Slope <- post$b[,1] - post$b[,2]
precis(post, depth=2)

# To better illustrate not just the difference in intercepts but how the weight
# contrast between the categories changes with H, 
hseq <- seq(-2.5, 2.5, len=50)
# I think what 'link' does here is give us the expected mean weight
# rather than a simulation over the posterior distribution?
sim_Ws_from_model_females <- link(m_HWS, data = list(H=hseq, S=rep(1,50)))
weight.PI.female <- apply(sim_Ws_from_model_females, 2, PI, prob=0.95)
sim_Ws_from_model_males <- link(m_HWS, data = list(H=hseq, S=rep(2,50)))
weight.PI.male <- apply(sim_Ws_from_model_males, 2, PI, prob=0.95)
xseq <- c(-2,2)
contrasts_weight <- sim_Ws_from_model_females - sim_Ws_from_model_males
plot(NULL, 
     xlim=xseq, ylim=c(-1, 1), 
     xlab="Height [std]", ylab="Weight contrast F-M [std]"
)
lines(apply(contrasts_weight, 2, median) ~ hseq, lwd=2)
for (p in c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99)) 
  shade(apply(contrasts_weight, 2, PI, prob=p), hseq)
abline(h=0, lty=2)

# The second thing might be to plot the range of the parameter estimates
xseq <- c(-2,2)
dens(post$a[,1], xlim=xseq, lwd=3, col=2, xlab="weight")
dens(post$a[,2], lwd=3, col=4, add=TRUE)
dens(post$b[,1], xlim=c(0,1), lwd=3, col=2, xlab="weight increase per height")
dens(post$b[,2], lwd=3, col=4, add=TRUE)

# Third, plot actual relationships
plot(NULL, 
     xlim=xseq, ylim=xseq, 
     xlab="Height [std]", ylab="Weight [std]"
     )
# Actual data points
points(d$W ~ d$H
       , col = c(alpha("violetred1",0.3), alpha("slateblue2", 0.3))[d$S]
       , lwd = 3)
# Posterior relationships for male and female
for (i in 1:100) 
  abline(a=post$a[i, 2], b=post$b[i,2], col=col.alpha("slateblue2",0.7))
for (i in 1:100) 
  abline(a=post$a[i, 1], b=post$b[i,1], col=col.alpha("violetred1",0.6))
# Prediction ranges
shade(weight.PI.female, hseq
      , col = col.alpha("violetred1", 0.08)
)
shade(weight.PI.male, hseq
      , col = col.alpha("slateblue2", 0.08)
)
lines(hseq, weight.PI.female[1,], lty=2, lwd=2
      , col=col.alpha("violetred1", 0.3))
lines(hseq, weight.PI.female[2,], lty=2, lwd=2
      , col=col.alpha("violetred1", 0.3))
lines(hseq, weight.PI.male[1,], lty=2, lwd=2
      , col=col.alpha("slateblue2", 0.3))
lines(hseq, weight.PI.male[2,], lty=2, lwd=2
      , col=col.alpha("slateblue2", 0.3))

# End of video 4 material:

# 'Full Luxury Bayes'

# Both models in a single model
data(Howell1)
d <- Howell1
d <- d[ d$age >= 18 , ]
Hbar <- mean(d$height)
d$S <- ifelse(d$male==1 , 2 , 1)
d$H <- d$height
d$W <- d$weight

m_SHW_full <- quap(
  alist(
    # weight
    W ~ dnorm(mu, sigma),
    mu <- a[S] + b[S]*(H-Hbar),
    a[S] ~ dnorm(60,10),
    b[S] ~ dunif(0,1),
    sigma ~ dunif(0,10),
    
    # height
    H ~ dnorm(nu, tau),
    nu <- h[S],
    h[S] ~ dnorm(160, 10),
    tau ~ dunif(0,10)
    
  ), data = d
)

post <- extract.samples(m_SHW_full)
n <- 1e4

with(post, {
  # simulate W for S=1
  H_S1 <- rnorm(n, h[,1], tau)
  W_S1 <- rnorm(n, 
                a[,1] + b[,1]*(H_S1-Hbar), sigma)
  # simulate W for S=2
  H_S2 <- rnorm(n, h[,2], tau)
  W_S2 <- rnorm(n, 
                a[,2] + b[,2]*(H_S2-Hbar), sigma)
  # compute contrast
  W_do_S <<- W_S2 - W_S1
})

# This is basically the counterfactual simulation again,
# we simulate both S=1 and S=2 data based on posteriors, 
# and calculate contrasts.
# This is called p(W|do(S)), meaning we 'do', i.e. intervene
# on S, and calculate the effect on W, which is thus
# the actual causal effect, implying causal through all arrows.

dens(W_do_S)
# (men are about 10 kg heavier)

# The difference to previous approach is:
# First approach was one stat model for each estimand

# Second approach is joint posterior for causal system,
# then simulate each estimand as an intervention.







## Finally,
## CHAPTER 6 and VIDEO 5 - > but put these back into 
## my previous script #5, because they are really about DAGs
## and multiple factors/multiple regression models

