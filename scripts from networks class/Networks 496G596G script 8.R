# Header-------------------------
## R script for ECOL496G 596G 
## Complex Systems and Networks
## UA class Spring 2025
## written by Anna Dornhaus

## This copy belongs to:

# Last R script for the class! 
# Before you read any more script stuff, clean up your final project ---------------------------------
# Roughly in this order:
# Check that you have a nice exploratory visualization of each of your networks.
#      - check out https://kateto.net/network-visualization
#      - did you use color deliberately to convey information or emphasis?
#      - did you remove labels that convey no interesting information, add labels
#         that do? You can use text() or mtext() to add a variety of labels.
#      - make sure you inspect the visualization to see if there may be errors.
#        Does it have the correct number of nodes and edges (roughly)? Does it 
#        seem to show the correct density or other patterns? If it looks suspicious,
#        make sure to do a little detective work to make sure you imported and 
#        plotted the right thing. 
#      - always have figure captions explaining what everything is but also what
#        your interpretation is and why. 
# Check that for each network, you have identified *some* interesting pattern
# or question, i.e. something that was unexpected or interesting when you looked
# at the above, or something you wondered based on the biological context.
#      - check that you have an additional graph, perhaps a data graph (i.e. 
#        boxplot or regression or histogram) that depicts quantitatively the 
#        interesting pattern
#      - if possible, think of a statistical test, perhaps by comparing to another
#        network or a random network, that makes it more convincing that your 
#        pattern is unexpected
# Make sure you are clear in your head and in writing on the page what *exactly*
#        the edges show. Are they yes/no or measuring a continuous quantity? What
#        is the temporal or spatial scale on which they are measured? Do they 
#        measure the quantity of interest directly or are they a proxy (e.g. 
#        proximity as a proxy for information exchange)? Never forget what the 
#        actually measured quantity is, and make sure you interpret any results
#        in that light, i.e. with appropriate uncertainty about what it means for
#        your conclusion. 
#        If your network data are from a database and no detailed original paper
#        describing these things is available, that's fine; but point it out.
# Make sure your script, submitted with the final project, has your name on it,
#        runs in order without error even if you first clear the environment and the plots, 
#        and has a comment at least every couple of lines and before every major
#        function or calculation stating what you mean to do. If you do have errors,
#        point out where they are, and add information on what class the variables
#        and functions you use are supposed to be, and what they should contain.
#        Also delete everything from the script that isn't necessary to make your
#        graphs or calculate the measures you need/use. (frequently save your file
#        under a distinctive name)
# When you consider the network measures to use, make sure you think about the
#        type of network you have. 'Betweenness' makes a lot of sense with unipartite,
#        unweighted networks. 'Flow' makes a lot of sense with weighted networks.
#        'Nestedness' (see below) makes a lot of sense with bipartite networks. A lot of other
#        combinations don't really work: if you have disconnected components, don't
#        use measures that calculate paths from every node to every other node. 
#        If you have bipartite networks, ask yourself whether using any measure 
#        that implies flow on a multi-step path makes any sense. Always try to link
#        network measures with their biological meaning.
#        If you coerce a weighted network into an unweighted one, or a directed into
#        undirected, or a temporal one into a static one, you are always not only 
#        losing information but also possibly introducing 'flow' where there isn't
#        any. Be cautious about interpretations in this case. 
# If you are not comparing multiple networks, often you just have a single outcome
#        (a number or a distribution) for your network measure. What does this mean?
#        Always think of a comparison. Is there a random network you can compare
#        to? A theoretically expected number? A different network? The same network
#        with something removed or added? Start from thinking about your biological
#        system to come up with good comparisons. Programming your own little 
#        simulation of how the network comes about can be a lot of help in figuring
#        out what you would 'expect' under different assumptions (see below).
# Machine learning methods can be a good way to detect pattern, e.g. groupings or
#        relevant traits, if you don't have ideas from visually inspecting the 
#        data or from first principles. See code no. 7. 

# --------------------------------

# Activate the Network Analysis packages ------------------
if ("circlize" %in% (.packages()))  detach("package:circlize", unload = TRUE)
library(sna)
library(intergraph)
library(igraph)
# And for color/graphing:
library(viridis)
library(scales)

# Importing dolphins, organizations, and zebras -----------------
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Lusseau dolphins/"
dolphins_df <- read.table(paste(path, "dolphinnetwork.csv", sep = ""), sep = ",")
dolphins_nodeinfo <- read.table(paste(path, "dolphinvertexinfo.csv", sep = ""), header = TRUE, sep = ",")
dolphins <- graph_from_adjacency_matrix(as.matrix(dolphins_df), mode = "undirected", weighted = TRUE)
V(dolphins)$sex <- dolphins_nodeinfo$dolphin.sex
V(dolphins)$disappeared <- dolphins_nodeinfo$dolphin.disappeared

data(emon)
emon1 <- emon[[1]]
organizations <- asIgraph(emon1)
V(organizations)$name <- V(organizations)$vertex.names

path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Rubenstein etal zebras/"
wildasses_raw <- read.table(paste(path, "wildasses.csv", sep = ""), header=T, sep = ",")
wildasses_indiv <- data.frame(ID = unique(wildasses_raw$Individual.ID))
for (i in 1:nrow(wildasses_indiv)) {
  groups_list <- list()
  for (r in 1:nrow(wildasses_raw)) {
    if(wildasses_raw$Individual.ID[r]==wildasses_indiv$ID[i]) {
      groups_list <- append(groups_list, wildasses_raw$Group.No[r])
    }
  }
  wildasses_indiv$groups[[i]] <- groups_list
}  
wildasses_interactions_m <- matrix(0, nrow = nrow(wildasses_indiv), ncol = nrow(wildasses_indiv))
colnames(wildasses_interactions_m) <- wildasses_indiv$ID
rownames(wildasses_interactions_m) <- wildasses_indiv$ID
for (r in 1:nrow(wildasses_indiv)) {
  for (c in 1:nrow(wildasses_indiv)) {
    weight <- 0
    weight <- length(intersect(wildasses_indiv$groups[[r]], wildasses_indiv$groups[[c]]))
    wildasses_interactions_m[r,c] <- weight
  }
  wildasses_interactions_m[r,r] <- 0 # Remove the loops
}
wildasses <- graph_from_adjacency_matrix(wildasses_interactions_m, weighted=TRUE, mode = "undirected")
wildassN <- delete_vertices(wildasses, "03-005T")

warnings()
color_closeness <- viridis(100)



# Generative models -----------------------------
# To understand what a type of analysis can and cannot do, it is often useful
# to 'simulate' data in a 'generative model'. The idea is that you can, for your
# simulated dataset, know exactly what the true answers/mechanisms/patterns
# are, because that's how you generated the data. Then you run your analysis to 
# see if your analysis correctly identifies the 'true' answers.

# In the context of networks, a lot of empirical networks suffer from being
# relatively sparse, that is, they only contain a small subset of the actually
# occurring interactions. 

## SIMBEES -----------------------------------
# Let's simulate a social network among worker bees. 
# Let's say that what we really want to capture is how information about a food
# source spreads among potential foragers, i.e. information transfer via 
# some kind of signals or cues exchanged between bees; and let's further assume
# these signals are not long-distance signals, i.e. the receiver has to come close
# to the sender to perceive them.

# Now we have to make a bunch of assumptions about how the system really works
# to simulate it; that's fine, right now we don't care what *actual* bees are 
# doing, we just want to see if our approach to analysis would recover whatever
# assumptions we're putting in.

# (1) Individuals have a small chance of independently searching for and finding 
# a food source.
# Let's say we have 1000 bees.
n <- 1000
# The chance of searching for and finding a food source with no information, 
# each morning, is
p_discovery <- 0.01
# (2) Individuals' enthusiasm for recruiting others is exponentially distributed.
individuals <- data.frame(recr_enthus = rexp(n, rate = 1))
# We created a data frame with one column, and filled that column with random
# numbers, exponentially distributed. See:
hist(individuals$recr_enthus)
# So most individuals have low enthusiasm for recruiting, but some have a high one.

# (3) The probability of an encounter leading to recruitment is proportional
# to the recruiters enthusiasm and the quality of the food source.
# Ok, so let's see who actually found something:
individuals$found_something <- ifelse(runif(n)<p_discovery, TRUE, FALSE)
summary(individuals$found_something) # This shows the successful bees for today.
# Now let's say the food source quality is normally distributed:
individuals$res_qual <- rnorm(n, mean = 1, sd = 0.5)
individuals$res_qual[individuals$res_qual<0] <- 0
# We just gave everyone a food resource quality, even the unsuccessful ones, but 
# that's fine, because:
adjustment_factor <- 1/10
individuals$recruiting_prob <- adjustment_factor * individuals$recr_enthus * individuals$res_qual * individuals$found_something
hist(individuals$recruiting_prob)
# All the 'FALSE', when multiplied, turned the value into zero.
# For the bees who found something, we now have a probability that is proportional
# to the resource quality they discovered and their individual enthusiasm.

# (4) Bees have a chance of recruiting another if they encounter them. 
# Encounters are spatially structured, i.e. you will only encounter individuals
# close to you.
# One could simulate all this in 2 dimensions and whatnot, making it pretty complex.
# But here, for simplicity, I just want to introduce a little spatial structure.
# So let's say the bees are arrayed on a 1-dimensional line as they are in the 
# data frame, and they can encounter everyone within
bee_density <- 100
# places of themselves in the line.
# To record who encounters whom, we now make an adjacency matrix.
bee_encounter_matrix <- matrix(data = NA, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    bee_encounter_matrix[i,j] <- ifelse(j<i+bee_density,
                                        ifelse(j>i-bee_density,
                                               ifelse(j==i, FALSE, TRUE),
                                               FALSE), 
                                        FALSE)
  }
}

bee_recruitment_matrix <- matrix(data = 0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    if(runif(1) < (bee_encounter_matrix[i,j]*individuals$recruiting_prob[i])){
      bee_recruitment_matrix[i,j] <- bee_recruitment_matrix[i,j] + 1
      # We'll also assume the recruits now learn the quality of the recruiter's resource:
      individuals$found_something[j] <- TRUE
      individuals$res_qual[j] <- individuals$res_qual[i]
    }
  }
}
# This is how many recruitment events we had:
sum(bee_recruitment_matrix)
# Now we can make a graph (i.e. a network) from this matrix:
bee_net_1stround <- graph_from_adjacency_matrix(bee_recruitment_matrix, mode = "directed")

# Ok but now the recruits have found the resource and can recruit themselves:
individuals$recruiting_prob <- adjustment_factor * individuals$recr_enthus * individuals$res_qual * individuals$found_something
for (i in 1:n) {
  for (j in 1:n) {
    if(runif(1) < (bee_encounter_matrix[i,j]*individuals$recruiting_prob[i])){
      bee_recruitment_matrix[i,j] <- bee_recruitment_matrix[i,j] + 1
      # We'll also assume the recruits now learn the quality of the recruiter's resource:
      individuals$found_something[j] <- TRUE
      individuals$res_qual[j] <- individuals$res_qual[i]
    }
  }
}
bee_net_2ndround <- graph_from_adjacency_matrix(bee_recruitment_matrix, mode = "directed")

# And let's assume in one day there is a third round of recruiting:
individuals$recruiting_prob <- adjustment_factor * individuals$recr_enthus * individuals$res_qual * individuals$found_something
for (i in 1:n) {
  for (j in 1:n) {
    if(runif(1) < (bee_encounter_matrix[i,j]*individuals$recruiting_prob[i])){
      bee_recruitment_matrix[i,j] <- bee_recruitment_matrix[i,j] + 1
      # We'll also assume the recruits now learn the quality of the recruiter's resource:
      individuals$found_something[j] <- TRUE
      individuals$res_qual[j] <- individuals$res_qual[i]
    }
  }
}
bee_net_3rdround <- graph_from_adjacency_matrix(bee_recruitment_matrix, mode = "directed")

bee_net <- bee_net_3rdround

ecount(bee_net_1stround)
ecount(bee_net_2ndround)
ecount(bee_net_3rdround)
plot(bee_net)
# Ta-da, now we have a network of who recruited whom, with quite a lot of recruiting
# happening.
par(mfrow=c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0))
color_quality <- viridis(round(10*max(individuals$res_qual)))
# For the plot, I am deleting the unconnected vertices - simply because not
# deleting them creates a really odd graph in which all the connected nodes are 
# squished together. (try it)
plot(delete_vertices(bee_net, v = V(bee_net)[degree(bee_net, V(bee_net))==0])
     # We could for example also color by resource quality found:
     #, vertex.color = alpha(color_quality[round(10*individuals$res_qual)], 0.3)
     # Or by location of bee:
     , vertex.color = alpha(viridis(n), 0.5)
     , vertex.size = 5
     , vertex.label.cex = 0.3
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0.5
     , edge.width = 2
     )
# Or a heatmap - but note that this is very sparse, and a large heatmap...
color_heatmap <- viridis(max(bee_recruitment_matrix))
heatmap(bee_recruitment_matrix
                 , col = c("slateblue4", "magenta")
                 , margins=c(8,8)
                 , scale = "none"
                 , cexRow = 0.5
                 , cexCol = 0.5
                 #, Rowv = TRUE
                 #, Colv = TRUE
                 , revC = TRUE
                 #, RowSideColors = magma()[degree(bee_net)]
                 #, ColSideColors = 
)
par(mfrow=c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))
hist(degree(bee_net, V(bee_net)))
## So far so good. We have artificially created a recruitment network. 

### Now the point of this exercise was to see how our sampling would affect our 
### understanding of the network.

## Network measurement --------------------------------

# Let's say we actually only measure proximity, not actual communication. 
# Then instead our measured network is like the encounter matrix:
proximity_net <- graph_from_adjacency_matrix(bee_encounter_matrix)
hist(degree(proximity_net, V(proximity_net)))
# Well that's really not great!

# Ok let's say instead that we do measure communication events, but we only 
# see about half of them. 
sampling_intensity <- 0.01
measured_matrix <- matrix(ncol = n, nrow = n)
for (i in 1:n) {
  for (j in 1:n) {
    measured_matrix[i,j] <- ifelse(runif(1)<sampling_intensity, 
                                   bee_recruitment_matrix[i,j],
                                   0)
  }
}
measured_net <- graph_from_adjacency_matrix(measured_matrix)
plot(delete_vertices(measured_net, v = V(measured_net)[degree(measured_net, V(measured_net))==0])
     # We could for example also color by resource quality found:
     #, vertex.color = alpha(color_quality[individuals$res_qual], 0.3)
     # Or by location of bee:
     , vertex.color = alpha(viridis(n), 0.5)
     , vertex.size = 5
     , vertex.label.cex = 0.3
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0.5
     , edge.width = 2
)
# Looks similar, although note there seems to be more of a separation between 
# the sections.
# The degree distribution still looks roughly the same, although it
# goes only half as far. 
hist(degree(measured_net, V(measured_net)))

## Network analysis -------------------------

# Now we should think about what it is we wanted to learn from studying this
# network. 
plot(degree(bee_net, V(bee_net))~individuals$recr_enthus
     , col = alpha(color_quality[round(10*individuals$res_qual)], 0.8)
     , pch = 1
     )
points(degree(measured_net, V(measured_net))~individuals$recr_enthus
       , col = alpha(color_quality[round(10*individuals$res_qual)], 0.5)
       , pch = 19
       )
legend("topleft"
       , c(paste("High (", 1.8, ")", sep=""), paste("Medium (", 0.9, ")", sep=""), paste("High (", 0.1, ")", sep=""))
       , col = c(color_quality[18], color_quality[9], color_quality[1])
       , pch = 19 
       )
# It's hard to see with both types of points in the same figure, but overall 
# I'd say we see two important things here:
# most individuals have the exact same color, i.e. same resource quality, i.e.
# are probably foraging from the same resource.
# Second, it looks like there are two major resources, and while there is a
# correlation of degree with inherent individual traits, it seems the better
# resource has 'recruited' more bees (degrees are higher for individuals foraging
# on the better resource).
# Both of the datasets show this if we have access to the individual trait, the 
# network, and the resource quality.
# However if we only had the degree, even with resource quality, that would not
# be a great estimator of individual propensity to recruit. 

# We can actually look at this quantitatively though:
individuals$disovered_res_qual <- individuals$res_qual*individuals$found_something
model_fit1 <- lm(degree(measured_net, V(measured_net), mode = "out")~ individuals$disovered_res_qual)
model_fit2 <- lm(degree(measured_net, V(measured_net), mode = "out")~individuals$recr_enthus)
summary(model_fit1)
summary(model_fit2)
# What we see is that both factors are clearly significant; but that 'resource quality'
# only has an R2 of about 3%, while individual 'recruitment enthusiasm' explains
# 28% of the variation in degree (i.e. number of successful recruitment events).
# 

## Thorough exploration of method robustness -----------------------------
# Go back and change the 'sampling intensity' to see how much this changes the result.
# In this example, it seems you will not get a qualitative change of the conclusion
# (that individual propensity to recruit matters more than what resource they discover)
# even if your sampling is a quite low fraction of the interactions that are taking place.

# You could even cycle over all the above several times and plot the robustness of 
# the effect size estimates as well as the R2 values against the sampling intensity...
summary(model_fit1)$adj.r.squared  # this is the R2 value
summary(model_fit1)$coefficients[2,1] # this is the slope, assuming one factor in the model
## ----------------------------

# Nestedness in bipartite networks ---------------------
# I copied the code from https://ieubascomptelab01.uzh.ch/bio365-fs22/session-03-23.html
# There are also packages that include this function, see below. 
compute_nestedness <- function(bipartite_network, mode = "total"){
  # Get number of rows and columns
  nrows <- nrow(bipartite_network)
  ncols <- ncol(bipartite_network)

  # Compute nestedness of rows
  nestedness_rows <- 0
  for(i in 1:(nrows-1)){
    for(j in (i+1): nrows){
      # Number of interactions shared by i and j (multiplying the two rows)
      c_ij <- sum(bipartite_network[i,] * bipartite_network[j,])      
      k_i <- sum(bipartite_network[i,])               # Degree of node i
      k_j <- sum(bipartite_network[j,])               # Degree of node j
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected
      o_ij <- c_ij / min(k_i, k_j)    # Overlap between i and j
      nestedness_rows <- nestedness_rows + o_ij
    }
  }
  
  # Compute nestedness of columns
  nestedness_cols <- 0
  for(i in 1: (ncols-1)){
    for(j in (i+1): ncols){
      c_ij <- sum(bipartite_network[,i] * bipartite_network[,j])      # Number of interactions shared by i and j
      k_i <- sum(bipartite_network[,i])               # Degree of node i
      k_j <- sum(bipartite_network[,j])               # Degree of node j
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected.
      o_ij <- c_ij / min(k_i, k_j)    # Overlap between i and j
      nestedness_cols <- nestedness_cols + o_ij         
    }
  }
  
  # Compute nestedness of the network
  nestedness <-(nestedness_rows+nestedness_cols)/((nrows*(nrows-1)/2)+(ncols*(ncols-1)/2))
  nestedness_rows <- nestedness_rows / (nrows*(nrows-1)/2)
  nestedness_cols <- nestedness_cols / (ncols*(ncols-1)/2)
#  if(mode=="rows") {return(nestedness_rows)}
#  if(mode=="cols") {return(nestedness_cols)}
#  return(nestedness)
  return(c(nestedness, nestedness_cols, nestedness_rows))
}

## Let's import the pollination network from Memmott again: ----------
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Memmott pollinator network"
pollination_df <- read.table(paste(path, "/memmott_1999.csv", sep = ""), header = TRUE, sep = ",")
pollination_df <- pollination_df[-1]
pollination_nodeinfo1 <- read.table(paste(path, "/plants.csv", sep = ""), sep = ",")
pollination_nodeinfo2 <- read.table(paste(path, "/pollinators.csv", sep = ""), header = TRUE, sep = ",")
colnames(pollination_df) <- pollination_nodeinfo1$V2
rownames(pollination_df) <- pollination_nodeinfo2$species
pollination <- graph_from_biadjacency_matrix(as.matrix(pollination_df), directed = FALSE, weighted = TRUE)

## Calculating nestedness! ------------------------------

# First we'll make an unweighed version of the pollination network, since the
# nestedness formula assumes this. 
unweighed_poll <- matrix(0, nrow = nrow(pollination_df), ncol = ncol(pollination_df))
unweighed_poll[pollination_df>0] <- 1
colnames(unweighed_poll) <- colnames(pollination_df)
rownames(unweighed_poll) <- rownames(pollination_df)

compute_nestedness(unweighed_poll, "cols")
compute_nestedness(unweighed_poll, "rows")
compute_nestedness(unweighed_poll, "total")
# What does this tell us? 
# Nestedness essentially calculates, for each node, the fraction of links that 
# the node with fewer links shares with each node with more links. 
# Roughly, if nodes overlap more, nestedness is higher, and if they overlap less,
# nestedness is lower.
# But more specifically than that, since it is only the node with fewer links we
# care about, the nodes with fewer links may share all their links with a node
# with many links, without the reverse being true. This is what 'nested'ness means:
# we're trying to see if all the low-connected nodes' networks are 'nested' inside
# the networks of the well-connected nodes. 

# Illustrating nestedness -----------------------
# Let's use the code from script 2 again to make a bipartite plot:
poll_web_as_matrix <- t(unweighed_poll)
# Or try:
# poll_web_as_matrix <- motten1982
# We order both plants and pollinators by degree:
plant_order <- rownames(poll_web_as_matrix)[order(rowSums(poll_web_as_matrix), decreasing = FALSE)] 
poll_order <- colnames(poll_web_as_matrix)[order(colSums(poll_web_as_matrix), decreasing = TRUE)] 
# We define what colors we'll use:
max_number_of_links_per_node <- max(colSums(poll_web_as_matrix))
intensity_color <- viridis(max_number_of_links_per_node)
plotweb(poll_web_as_matrix
        , method = "normal"
        , sequence = list(seq.high=poll_order, seq.low=plant_order)
        , text.rot = 90
        , col.high = intensity_color[colSums(poll_web_as_matrix)]
        , col.low = viridis(max(rowSums(poll_web_as_matrix)))[rowSums(poll_web_as_matrix)]
        , col.interaction = intensity_color[colSums(poll_web_as_matrix)]
        , bor.col.interaction	= intensity_color[colSums(poll_web_as_matrix)]
        , text.high.col = "white" # the pollinator labels are uninformative
        , empty = FALSE
)
# If the network was extremely 'nested', we would not see the diagonal links from 
# the bottom left to the top right; instead we would see all the specialist pollinators
# visiting wild carrot, and all the specialist plants being visited just by the 
# generalist pollinators on the left. This is true to some degree but not universal.



## --------------------------------

# Scale-free property ------------------
# Scale-free networks are defined as ones with a power-law degree distribution, 
# meaning that if you plot the degree distribution (degree on x-axis, number of
# nodes that have it on y) on a log-log scale, you get a straight (decreasing) line.
# This would show that the degree distribution conforms to freq(deg) = deg^(-gamma)
# where gamma is an exponent that is also of interest. 
# Generally it does not make sense to even calculate or suspect this unless you have
# a very large network in which the degree of different nodes in fact varies over
# orders of magnitude.
# If this pattern is present, various other properties can be derived mathematically
# (e.g. see https://networksciencebook.com/chapter/4#power-laws).
# For the empiricist, the most important element here (IMHO) is that some (though
# not all!) networks have a few extremely well-connected nodes, and many poorly
# connected nodes, and this can have consequences for robustness, the small-world
# -property, etc. 
# We already looked at degree distributions of smaller networks previously. 
library(networkdata)
data(cent_lit)
data(netsci)
data(powergrid) # This is quite large (>4K nodes)
data(sn_auth)

# Here is a larger one:
large_net <- powergrid
deg_dist <- hist(degree(large_net, V(large_net)))
deg_dist <- data.frame(freq = deg_dist$counts, degr = deg_dist$mids)
plot(freq~degr
     , data = deg_dist
     , col = "purple"
     , pch = 19
     , log = "xy"
)
# If you want to add a linear fit, it's easier to calculate the logs yourself:
deg_dist$logfreq <- log(deg_dist$freq)
deg_dist$logdegr <- log(deg_dist$degr)
deg_dist$logdegr[deg_dist$logdegr==-Inf] <- NA
deg_dist$logfreq[deg_dist$logfreq==-Inf] <- NA
plot(logfreq~logdegr
     , data = deg_dist
     , col = "purple"
     , pch = 19
     , xlab = "Log(Degree of nodes)"
     , ylab = "Log(Number of nodes with this degree)"
)
data_for_model <- deg_dist[complete.cases(deg_dist),]
fit <- lm(logfreq~logdegr
          , data = data_for_model)
abline(fit, xpd=FALSE, lwd = 2, col = "darkblue")
summary(fit)
# Note that this approach is valid if you want to fit an exponent for the power
# -law of degree distribution; but most log-log distributions look very linear, 
# so eyeballing this or fitting a line is not a test of whether a power-law is
# the best-fit relationship here. I.e. this is not a test whether or not the
# relationship is strictly scale-free.


## -----------------------------------

# Time-ordered or dynamic networks ------------------
#install.packages("timeordered")
library(timeordered)
data("ants")
dynamic_net <- generatetonetwork(ants)
plottonet(dynamic_net)
static_net <- generatetimeaggregatednetwork(dynamic_net, 0, 1500)
plot(static_net)
# Slicing the total duration into 16 intervals
timeslices <- data.frame(from = seq(from = 0, to = 1500, length.out = 16), to = seq(from = 100, to = 1600, length.out = 16))
# Alternatively, you could make cumulative time intervals
timeslices <- data.frame(from = rep(0, length.out = 16), to = seq(from = 100, to = 1600, length.out = 16))

sliced_nets <- generatenetworkslices(dynamic_net, timeslices)
par(mfrow = c(4,4), mar = c(0,0,0,0), oma = c(0,0,0,0))
total_nodes <- vcount(static_net)
nodecolors <- viridis(total_nodes)
V(static_net)$color <- nodecolors
diam <- numeric(0)
avg_deg <- numeric(0)
avg_betw <- numeric(0)
for (i in 1:16) {
  plot(sliced_nets[[i]]
       , vertex.color = alpha(nodecolors, 0.5)
       , vertex.size = 10
       , vertex.label.cex = 0.3
       , vertex.label.color = "black"
       , vertex.label.family = "Arial"
       , edge.arrow.size = 0.1
       , edge.width = 2
       , layout = layout_in_circle(sliced_nets[[i]])
  )
  diam[i] <- diameter(sliced_nets[[i]])
  avg_deg[i] <- mean(degree(sliced_nets[[i]], v = V(sliced_nets[[i]])))
  avg_betw[i] <- mean(betweenness(sliced_nets[[i]], v = V(sliced_nets[[i]])))
}
warnings()
# Other than being able to actually look at processes over time, and/or time-ordered
# communication on the network, this also allows us to see how the network measures
# will be affected by our choice of sampling. 
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0), xpd = FALSE)
plot(NULL
     , xlim = c(1, 16)
     , ylim = c(0, 30)
     , xlab = "Length of sampling period"
     , ylab = "Measure"
     )
points(diam
       , col = "darkblue"
       , pch = 19
       )
points(avg_betw
       , col = "purple"
       , pch = 19
)
points(avg_deg
       , col = "darkgreen"
       , pch = 19
)
legend("bottomright"
       , c("Diameter", "Avg. betweenness", "Avg. degree")
       , col = c("darkblue", "purple", "darkgreen")
       , pch = 19 
)


