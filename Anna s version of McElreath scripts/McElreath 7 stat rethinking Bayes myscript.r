## CODE TO GO WITH BOOKCLUB
## ON Richard McElreath's 'Statistical Rethinking'

########### PART VII #####################

library(rethinking)
library(scales)
library(dagitty) # for drawing DAGs
#dev.off()
par(mar=c(5,5,1,1), oma=c(0,0,0,0)) #b,l,t,r

# CHAPTER 6
# starting 6.2 Post-treatment bias
# Video 6 (Bad and good controls)

# Fungicide example from video 5 on Post-treatment bias, i.e. 
# 'Pipes'.
set.seed(71)
# number of plants
N <- 100
# simulate initial heights
h0 <- rnorm(N,10,2)
# assign treatments and simulate fungus and growth
treatment <- rep(0:1 , each=N/2)
fungus <- rbinom(N, size=1 , prob=0.5 - treatment*0.4)
h1 <- h0 + rnorm(N, 5 - 3*fungus)
# compose a clean data frame
d <- data.frame( h0=h0 , h1=h1 , treatment=treatment , fungus=fungus )
precis(d)

# Thinking more about priors
# logarithic priors to ensure the parameter is positive
sim_p <- rlnorm( 1e4 , 0 , 0.25 )
precis( data.frame(sim_p) )

# The point of post-treatment bias is that we should not include fungus
# in the model

m6.8 <- quap(
              alist(
                h1 ~ dnorm( mu , sigma ),
                mu <- h0 * p,
                p <- a + bt*treatment,
                a ~ dlnorm( 0 , 0.2 ) ,
                bt ~ dnorm( 0 , 0.5 ),
                sigma ~ dexp( 1 )
              ), data=d )

precis(m6.8)

post <- extract.samples(m6.8, n=1000)
dens(post$bt)

# Didn't actually write a plot, maybe will at some point
plot(NULL 
#    , xlim=xseq, ylim=c(-1, 1) 
     , xlab="Treatment", ylab="Growth of plant"
)

# Keep in mind that model selection doesn't solve these problems (post treatment bias
# or collider bias).
# Quote "As argued in Chapter 1, prediction and causal inference are just not
# the same task. No statistical procedure can substitute for scientific knowledge 
# and attention to it."



# Causal Inference and 'DO-CALCULUS'

# Important point in video 6:
# The coefficient of X on Y is not sufficient to estimate the effect
# of X on Y. Instead, the 'do(X)' mechanism should be used, i.e. 
# calculating the effect of X on Y averaged over the distribution of 
# U, the confound that affects both. 

# 'Backdoor criterion'
# 1 Identify all paths connecting X to Y (follow arrows in either direction)
# 2 Paths with arrows entering X are backdoor paths (which you want to cut)
# 3 Find adjustment set that closes backdoor paths
#    i.e. try to find control variables that block the backdoor paths
dag_backdoor <- dagitty( "dag {
X -> Y
U -> Y
U -> Z
Z -> X
}")
coordinates(dag_backdoor) <- list( x=c(X=0,Z=1,U=2,Y=3) , y=c(X=0,Z=1,U=2,Y=0))
drawdag(dag_backdoor)
# The backdoor path goes through U and Z. If we can measure Z, 
# then that is sufficient to 'close' the backdoor path, and it's
# ok if we can't measure U. 
par(mar=c(5,5,1,1), oma=c(0,0,0,0)) #b,l,t,r

# So to estimate X -> Y, we can stratify by Z, and block all the influence
# of U (because all of it goes through Z).

# 'Stratify' means we calculate x->Y for each value of Z, and then
# average over the distribution of z's.
# This is the same as a linear model
# mu = a + bX + cZ,  Y~Normal(mu, sigma)

# Pre-defined happiness simulation
d <- sim_happiness( seed=1977 , N_years=1000 )
precis(d)

# page 190 in the book: general recipe to close paths and statistically achieve
# a controlled trial 



