# Header-------------------------
## R script for ECOL496G 596G 
## Complex Systems and Networks
## UA class Spring 2025
## written by Anna Dornhaus

## This copy belongs to:

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

## Degree distributions -------------
# These are the three networks we prepared above:
hist(degree(dolphins, V(dolphins)))
hist(degree(organizations, V(organizations)))
hist(degree(wildassN, V(wildassN)))
# Now we'll do n random networks, and compare the degree distributions
n <- 100
net <- as.undirected(dolphins)
no_nodes <- length(V(net))
random_degs <- data.frame(matrix(nrow = no_nodes, ncol = 0))
random_net <- sample_gnm(vcount(net), ecount(net)
                         , directed = FALSE, loops = FALSE)
random_degs <- cbind(random_degs, degree(random_net, V(random_net)))
rownames(random_degs) <- V(net)$name
random_deg_dist <- data.frame(matrix(data = NA, ncol = 1+no_nodes, nrow = n))
colnames(random_deg_dist) <- seq(from = 0, to = no_nodes)
for (i in 1:(n-1)) {
  random_net <- sample_gnm(vcount(net), ecount(net)
                           , directed = FALSE, loops = FALSE)
  random_degs <- cbind(random_degs, degree(random_net, V(random_net)))
  for (j in 0:no_nodes) {
    random_deg_dist[i,(j+1)] <- length(which(random_degs[,i]==j))
  }
}
# This is the real degree distribution:
real <- hist(degree(net, V(net)), breaks = seq(0,no_nodes,l=no_nodes+1))
# This is the distribution of node-distributions in the random networks:
boxplot(random_deg_dist
        , range = 0
        , xlab = "Degree of nodes in network"
        , ylab = "Frequency")
# Adding the real one for comparison:
points(real$counts
       , pch = 19
       , col = alpha("darkblue", 0.5)
)

### Closeness distributions ----------------------
# You can do this with any other node-based information - let's do closeness:
n <- 100
net <- as.undirected(organizations)
no_nodes <- length(V(net))
close_breaks <- round(seq(0, 0.1,l=100), digits = 3)
no_breaks <- length(close_breaks)
random_close <- data.frame(matrix(nrow = no_nodes, ncol = 0))
random_net <- sample_gnm(vcount(net), ecount(net)
                         , directed = FALSE, loops = FALSE)
random_close <- cbind(random_close, round(closeness(random_net, V(random_net)), digits = 3))
rownames(random_close) <- V(net)$name
random_close_dist <- data.frame(matrix(data = NA, ncol = no_breaks, nrow = n))
colnames(random_close_dist) <- close_breaks
for (i in 1:(n-1)) {
  random_net <- sample_gnm(vcount(net), ecount(net)
                           , directed = FALSE, loops = FALSE)
  random_close <- cbind(random_close, 
                        round(closeness(random_net, V(random_net)), digits = 3))
  for (j in 1:no_breaks) {
    random_close_dist[i,j] <- length(which(random_close[,i]==close_breaks[j]))
  }
}
# This is the real closeness distribution:
real <- hist(round(closeness(net, V(net)), digits = 3), breaks = close_breaks)
# This is the distribution of node-distributions in the random networks:
boxplot(random_close_dist
        , range = 0
        , xlab = "Closenesses of nodes in network"
        , ylab = "Frequency")
# Adding the real one for comparison:
points(real$counts
       , pch = 19
       , col = alpha("darkblue", 0.5)
)


## Bootstrapping in general ----------------------
# The idea of bootstrapping is simply to delete one or more datapoints randomly
# from your dataset, and recalculate the value of whatever you are measuring.
# By doing this many times (each time using a different datapoint or datapoints),
# you can see how sensitive your measure is to the particular set of points you
# used. 

# For example, we can check how accurate our estimation of the mean or median is:
testdata <- runif(100, min = 0, max = 100)
#testdata <- rnorm(100, mean = 50, sd = 20)
#testdata <- rexp(100, rate = 1) #if using this make the range of the boxplot 0-2
means <- numeric(0)
medians <- numeric(0)
for (i in 1:200) {
  testsample <- sample(testdata, 50) # Here we are dropping fully 50% of the samples
  means <- c(means, mean(testsample))
  medians <- c(medians, median(testsample))
}
boxplot(testdata, means, medians
        , col = c("grey", "darkred", "darkseagreen")
        , names = c("Data", "Means", "Medians")
        , ylim = c(30,70)
        , range = 0
)
abline(h=mean(testdata))
abline(h=median(testdata))
# Despite dropping a lot of data, our means and medians remain pretty accurate. 
## Bootstrapping network measures -----------------------
### Bootstrapping 'by hand' ----------------------------------

n <- 100
boot_prop <- 0.6
net <- dolphins
no_boot_nodes <- round(length(V(net))*boot_prop, digits = 0)
if (no_boot_nodes==length(V(net))) {no_boot_nodes<-length(V(net))-1}
# Ok now we know how many nodes we want in our 'bootstrapped' networks
# Which means we have to remove:
to_remove <- length(V(net)) - no_boot_nodes

boot_diam_dist <- numeric(0)
for(i in 1:n) {
  to_delete <- sample(V(net), to_remove, replace = FALSE, prob = NULL)
  boot_net <- delete_vertices(net, v = to_delete)
  diam <- diameter(boot_net)
  boot_diam_dist <- c(boot_diam_dist, diam)
}
boxplot(boot_diam_dist
        , xlab = "Distribution of diameters in bootstrapped data"
        , ylab = "Diameter"
        , ylim = c(0, 1+max(boot_diam_dist))
        , range = 0
        , col = "darkseagreen4"
        )
abline(h=diameter(net), col = alpha("darkblue",0.5), lwd = 3)
# What this does is it gives you an idea how much the estimate / measure of
# diameter depends on the particular sample of nodes you have. 
# This is implying that you don't necessarily actually have all the nodes
# or all the interactions correctly or completely, which is the case in many
# empirical networks. 
# If the bootstrapped measure deviates a lot from the overall, or has a large
# spread, it may be worth investigating which nodes, when removed, have the 
# largest effect, perhaps leading into an assessment of robustness. 


### Bootstrapping with package --------------------
library(snowboot)
net <- as.undirected(dolphins, mode = "collapse")
real <- hist(degree(net, V(net)), breaks = seq(0,vcount(net),l=vcount(net+1)))
plot(boot_dd(lsmi_dd(net = igraph_to_network(net)), B = 100))
# Here, the package snowboot automatically uses bootstrapping to estimate 
# the underlying degree distribution (green boxplot) compared to the actual 
# empirical degree distribution. 
points(real$counts/vcount(net)
       , pch = 19
       , col = alpha("darkblue", 0.5)
)
abline(v=mean(degree(net))
         , col = alpha("darkblue", 0.5)
       , lwd = 2
       )
# Another package that can do this with more different parameters (not just degree)
# and more options is 'bootnet'. An accompanying tutorial paper is 
# 'Estimating psychological networks and their accuracy: A tutorial paper'
# Behav Res 50, 195–212 (2018), by Epskamp, S., Borsboom, D. & Fried, E.I. 

## Robustness curves ---------------------------------
## The idea behind 'robustness curves' is to see how some network measure we
## are interested in changes or breaks down as nodes are removed. 
## The curves start at x=0 with the original network, and then as x increases,
## more and more nodes are removed, with y showing how network functionality
## changes. 
### Importing some real examples to use later ----------
# If you haven't installed the package 'networkdata' yet, use this:
#install.packages("remotes")
#remotes::install_github("schochastics/networkdata")
library(networkdata)
#data(package = "networkdata")
data(cent_lit)
data(movie_103)
data(netsci)
data(powergrid) # This is quite large (>4K nodes)
data(sn_auth)
### Plotting robustness of several network measures ------------
## Roughly, what we do here is we initialize a vector for each measure we
## care about; then we use a for-loop, and each time we go through the loop,
## we first calculate the network measure for the current network, and then
## remove an additional node from the network. 
## This way, at the end, we have a vector that shows the value for the network
## measure for each number of nodes removed.

#start_net <- powergrid
start_net <- cent_lit
red_net <- start_net
diam <- vector(mode = "numeric", length = vcount(start_net))
eff <- vector(mode = "numeric", length = vcount(start_net))
shortest_dist <- vector(mode = "numeric", length = vcount(start_net))
avg_degree <- vector(mode = "numeric", length = vcount(start_net))
comp_size <- vector(mode = "numeric", length = vcount(start_net))
for(i in 1:(vcount(start_net)-1)) {
  diam[i] <- diameter(red_net)
  eff[i] <- global_efficiency(red_net)
  shortest_dist[i] <- mean_distance(red_net)
  avg_degree[i] <- mean(degree(red_net))
  comp_size[i] <- max(components(red_net)$csize)
  vertex_deleted <- sample(V(red_net), 1)
  if (vertex_deleted==13) (dude_node <- i)
  red_net <- delete_vertices(red_net, v = vertex_deleted)
}
## We've collected the data we want - note the above may take a while for a
## big network, and you may want to just calculate one measure to being with
## in that case.
## Now we plot the robustness curves:
cols <- viridis(5)
plot((shortest_dist/max(shortest_dist, na.rm=T))
     , xlab = "Number of nodes removed +1"
     , ylab = "Measure"
     , ylim = c(0, 1)
     , pch = 19
     , col = cols[1]
)
points((diam/max(diam, na.rm = T))
       , pch = 19
       , col = cols[2]
)
points((eff/max(eff, na.rm = T))
       , pch = 19
       , col = cols[3]
)
points((avg_degree/max(avg_degree, na.rm = T))
       , pch = 19
       , col = cols[4]
)
points((comp_size/max(comp_size, na.rm = T))
       , pch = 19
       , col = cols[5]
)
legend("bottomleft"
       , legend = c("shortest dist", "diameter", "efficiency", "avg degree", "largest component")
       , col = cols, pch = 19
       )
abline(v = dude_node)
## Something interesting can be seen here, which is that different measures
## show quite different patterns. This is what you want to see from a robustness
## curve.
## If the curve essentially goes linearly from top left to bottom right, then
## the network measure just decreases the same amount for each node removed.
## If the curve drops steeply first, staying under the diagonal, then the measure
## is very sensitive to even a few nodes being removed.
## If the curve stays above the diagonal, perhaps even staying pretty flat
## until a breakdown towards the right of the graph, then we see a pretty robust
## network. 

## For the measures I included above however, another pattern emerges: some curves
## actually increase from the starting point. This illustrates a problem that can
## emerge with many network measures, namely that they are either not properly 
## defined for a disconnected network, or that their values are not comparable
## across networks of different sizes. 
## As you remove nodes, eventually you get several disconnected subnetworks.
## With the problems above, you either get NAs - which in the plot above I am 
## simply ignoring - or you get values (e.g. for avg shortest distance) that are
## just calculated for subnetworks. 

### Plotting robustness of network *performance* measures ------------
## Given the above, it would be better to think about what is actually a 
## relevant performance measure in the biological system you are interested in.
## Is it information transmission? Is it just being connected with everyone
## at all? Is it subgroup size? Or something else? 
## Then try to identify the operational measure that relates to the phenomenon
## that matters, and think about how you want to treat cases of disconnected 
## nodes or subnetworks. 

## If you google-scholar 'network robustness measures', you'll find a series
## of papers proposing different robustness measures, often condensing the 
## process we did above into a single measure (e.g. measuring how many nodes 
## you can remove before the network breaks into multiple components). 

## If you are interested in bipartite networks, I recommend this article on 
## robustness of ecological networks: https://www.r-bloggers.com/2024/09/estimating-ecological-network-robustness-with-r-a-functional-approach/
## Another useful tutorial is in: https://ieubascomptelab01.uzh.ch/bio365-fs24/index.html

start_net <- cent_lit

# One common measure is total 'flow' through the network, which is measuring
# the number of distinct paths one can go through to get from each node to 
# every other node. 
flow_measure <- function(a_network) {
flow_matrix <- matrix(data = NA, nrow = vcount(a_network), ncol = vcount(a_network))
for (i in (1:vcount(a_network))) {
  for (j in (1:vcount(a_network))) {
    if (i!=j) {
      flow_matrix[i,j] <- max_flow(a_network
                                   , V(a_network)[i], V(a_network)[j]
                                   #, capacity = E(start_net)$weight
                                   )$value
    }
  }
}
total_flow <- sum(flow_matrix, na.rm = T)
return(total_flow)
}
# For a network with ~100 nodes, this takes a while to calculate, so doing the
# entire robustness graph will take even longer.

# Another measure of overall network functionality is whether the network has
# broken down into unconnected components. This of course only makes sense to
# calculate if the graph is not disconnected to begin with; and even if it is
# not, it is a useful measure only for fairly 'dense' networks. 
start_net <- wildassN
# Density measures the number of edges over the number of possible edges.
edge_density(start_net, loops = FALSE)
# Vertex connectivity measures how many nodes you have to remove to make
# the network disintegrate into disconnected components.
vertex_connectivity(start_net, source = NULL, target = NULL, checks = TRUE)
# This function counts the number of disconnected components in the network:
count_components(start_net, mode = "weak")

### Plotting robustness curves while bootstrapping of node removal order -------------
## In all the above, we just randomly removed nodes until there were none
## left. However, how quickly a network or network measure breaks down may depend
## critically on which nodes are removed first, or more generally the order
## of network removal. To make a general statement about robustness, we would 
## therefore want to randomize this order, and see the distribution of resulting
## curves. 
## Here, I am putting everything we did above into another loop, which generates
## a new random order of node removal each time. All results are saved, we thus
## get a matrix, instead of a vector, of results, with each column reflecting 
## a different order of node removal, and each row containing the network measure
## as it is 

start_net <- cent_lit
# Remember that below you are doing several matrices with n * no of nodes entries.
# So make sure either the network is not too large or n is small or you don't
# calculate all the measures below at once
red_net <- start_net
n <- 10
comp_number <- matrix(data = NA, nrow = vcount(start_net), ncol = n)
comp_size <- matrix(data = NA, nrow = vcount(start_net), ncol = n)
#flow <- matrix(data = NA, nrow = vcount(start_net), ncol = n)
for (j in 1:n) {
  removal_order <- sample(V(start_net)$name, vcount(start_net), replace = FALSE)
  for(i in 1:(vcount(start_net)-2)) {
    comp_number[i,j] <- count_components(red_net, mode = "weak")
    comp_size[i,j] <- max(components(red_net)$csize)
#    flow[i,j] <- flow_measure(red_net)
    red_net <- delete_vertices(start_net, v = removal_order[1:i])
    print(paste(i, ",", j))
  }
}
point_color_darkness <- 0.1
plot(NULL
     , xlab = "Number of nodes removed +1"
     , ylab = "Measure"
     , ylim = c(0, 1)
     , xlim = c(0, vcount(start_net))
)
for (j in 1:n) {
#  points((flow[,j]/max(flow, na.rm = T)) ~ seq(1:vcount(start_net))
#         , pch = 19
#         , col = alpha(cols[1], point_color_darkness)
#  )
  points((comp_number[,j]/max(comp_number, na.rm = T)) ~ seq(1:vcount(start_net))
         , pch = 19
         , col = alpha(cols[2], point_color_darkness)
  )
  points((comp_size[,j]/max(comp_size, na.rm = T))~ seq(1:vcount(start_net))
         , pch = 19
         , col = alpha(cols[4], point_color_darkness)
  )
}
legend("topright"
       , legend = c("No. of components", "Size of largest component")
       , col = c(cols[2], cols[4]), pch = 19
)

# Robustness under targeted removal  -----------------

start_net <- cent_lit
red_net <- start_net
# Initialize vectors
cn <- count_components(start_net, mode = "weak")
cs <- max(components(start_net)$csize)
comp_number <- rep(cn, length.out = vcount(start_net))
comp_size <- rep(cs, length.out = vcount(start_net))
# Define a list of all the nodes in the network but in the order we want to 
# remove them
removal_order <- V(start_net)[order(betweenness(start_net))]
  # This would be random order: 
  # sample(V(start_net)$name, vcount(start_net), replace = FALSE)
for(i in 1:(vcount(start_net)-1)) {
  # Remove the next node from the list
  red_net <- delete_vertices(start_net, v = removal_order[1:i])
  comp_number[i+1] <- count_components(red_net, mode = "weak")
  comp_size[i+1] <- max(components(red_net)$csize)
}
cols <- viridis(5)
point_color_darkness <- 0.5
plot(NULL
     , xlab = "Number of nodes removed"
     , ylab = "Measure"
     , ylim = c(0, 1)
     , xlim = c(0, vcount(start_net))
)
points((comp_number/max(comp_number, na.rm = T)) ~ seq(0:(vcount(start_net)-1))
       , pch = 19
       , col = alpha(cols[1], point_color_darkness)
)
points((comp_size/max(comp_size, na.rm = T))~ seq(0:(vcount(start_net)-1))
         , pch = 19
         , col = alpha(cols[3], point_color_darkness)
)
legend("topright"
       , legend = c("No. of components", "Size of largest component")
       , col = c(cols[1], cols[3]), pch = 19
)
# There is a whole field of network science concerned with how best to pick 
# nodes to remove to do maximal damage to a network. 
# One alternative strategy to the above might be to calculate the betweenness
# on the reduced network, i.e. to pick whatever node NOW has the highest betweenness
# rather than in the order of the ORIGINAL betweennesses. 

start_net <- cent_lit
red_net <- start_net
# Initialize vectors
cn <- count_components(start_net, mode = "weak")
cs <- max(components(start_net)$csize)
comp_number <- rep(cn, length.out = vcount(start_net))
comp_size <- rep(cs, length.out = vcount(start_net))
for(i in 1:(vcount(start_net)-1)) {
  # Determine which node to remove
  to_remove <- V(red_net)[which.max(betweenness(red_net))]
  # Remove the next node from the list
  red_net <- delete_vertices(red_net, v = to_remove)
  comp_number[i+1] <- count_components(red_net, mode = "weak")
  comp_size[i+1] <- max(components(red_net)$csize)
  print(i)
}
cols <- viridis(5)
point_color_darkness <- 0.5
plot(NULL
     , xlab = "Number of nodes removed"
     , ylab = "Measure"
     , ylim = c(0, 1)
     , xlim = c(0, vcount(start_net))
)
points((comp_number/max(comp_number, na.rm = T)) ~ seq(0:(vcount(start_net)-1))
       , pch = 19
       , col = alpha(cols[1], point_color_darkness)
)
points((comp_size/max(comp_size, na.rm = T))~ seq(0:(vcount(start_net)-1))
       , pch = 19
       , col = alpha(cols[3], point_color_darkness)
)
legend("topright"
       , legend = c("No. of components", "Size of largest component")
       , col = c(cols[1], cols[3]), pch = 19
)


# -----------------------------


