# Header-------------------------
## R script for ECOL496G 596G 
## Complex Systems and Networks
## UA class Spring 2025
## written by Anna Dornhaus

## This copy belongs to:

# Activate the Network Analysis packages ------------------
library(sna)
library(intergraph)
library(igraph)
# And for color/graphing:
library(viridis)
library(scales)

# Importing dolphins & organizations -----------------
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Lusseau dolphins/"
dolphins_df <- read.table(paste(path, "dolphinnetwork.csv", sep = ""), sep = ",")
dolphins_nodeinfo <- read.table(paste(path, "dolphinvertexinfo.csv", sep = ""), header = TRUE, sep = ",")
dolphins <- graph_from_adjacency_matrix(as.matrix(dolphins_df))
V(dolphins)$sex <- dolphins_nodeinfo$dolphin.sex
V(dolphins)$disappeared <- dolphins_nodeinfo$dolphin.disappeared

data(emon)
emon1 <- emon[[1]]
organizations <- asIgraph(emon1)
V(organizations)$name <- V(organizations)$vertex.names

# Actually analyzing and visualizing results -------------

# We already know how to calculate some network statistics, both for the whole
# network, like
diameter(dolphins) 
# and for individual nodes, e.g.
closeness(dolphins)
# But what do these mean?

## Regular statistical tests ----------------
# Sometimes we have other biological context that means we want to relate these
# numbers to something else. E.g. we might wonder if closeness is related to sex
# in dolphins. 
# In that case, we can treat the network statistic (here: closeness) just like 
# any other empirical measure. In this example, we're probably comparing the 
# closeness across two sexes, so we need a two-group statistical comparison 
# (typically a Mann-Whitney-U-Test) and a boxplot to illustrate.

### Comparing median trait across two groups: -----------------
# Confusingly, the U-Test is actually called 'wilcox.test' in R - 
# (for the purpose of this example, I am ignoring the dolphins of 'unknown' sex)
dolphins_malefemale_only <- as.factor(V(dolphins)$sex[V(dolphins)$sex=="M" | V(dolphins)$sex=="F"])
closeness_malefemale_only <- closeness(dolphins)[V(dolphins)$sex=="M" | V(dolphins)$sex=="F"]
wilcox.test(closeness_malefemale_only ~ dolphins_malefemale_only)
# The p-value is 0.09, which means there is a 9% chance of getting medians this different
# or more different even if dolphin sex has no impact on their 'closeness'.
# We call this 'not significant', and accept the null hypothesis that sex has
# no impact. 
# Then we already know how to do the plot: 
par(mfrow=c(1,1), oma = c(0,0,0,0), mar = c(2,6,1,1)
    , mgp=c(4, 1, 0), las=1) # bottom, left, top, right
colors_dolphinsex <- c("darkseagreen", "darkmagenta")
y_max <- 0.010
y_axis_offset <- y_max * 0.95 # This puts sample size numbers 5% below max
Nice_Plot <- boxplot(closeness_malefemale_only ~ dolphins_malefemale_only
                     , xlab = ""
                     , ylab = "Closeness [calculated over entire\n network regardless of sex]" 
                     , names = c("Males", "Females")
                     , col = alpha(colors_dolphinsex, 0.5)
                     , ylim = c(0, y_max)
)
nbGroup <- nlevels(as.factor(Nice_Plot$names))
text( 
  x=c(1:nbGroup), 
  y=y_axis_offset,
  cex = 1,
  paste(Nice_Plot$n,sep="")  
)
stripchart(closeness_malefemale_only ~ dolphins_malefemale_only
           , add = TRUE
           , pch = 19
           , col = colors_dolphinsex
           , method = "jitter"
           , jitter = 0.2
           , vertical = TRUE,
)

### Relating two continuous measures: -----------------
# Alternatively, we might have information about nodes that is a continuous measure.
# In that case, we use a regression: For example, the 'Paid.Staff' in the
# organizations dataset. 
stat_model <- lm(closeness(organizations) ~ V(organizations)$Paid.Staff)
summary(stat_model)
# This is also not significant (the p-value is 0.65). Let's illustrate this:
par(mfrow=c(1,1), oma = c(0,0,0,0), mar = c(4,4,1,1)
    , mgp=c(3, 1, 0), las=1, xpd=FALSE) # bottom, left, top, right
color_closeness <- viridis(max(100*closeness(organizations), na.rm=T))
plot(closeness(organizations) ~ V(organizations)$Paid.Staff
     , pch = 19 # set point shape
     , col = color_closeness[100*closeness(organizations)]
     # This is a little longer than it would be if we just had used the subsets of data. I like it
     # better because it keeps flexibly using the entire dataset. 
     , cex = 1.5 # point size - 1 is default
     , xlab = "Number of paid staff"
     , ylab = "Closeness in network"
     , ylim = c(0, 0.1)
)
text(V(organizations)$Paid.Staff+10, closeness(organizations)+0.003, labels=V(organizations)$name, cex= 0.7)
# For non-significant regressions, we normally don't add a line. But if it had been 
# significant, we might have added
abline(stat_model
       , lwd = 2
       , col = "magenta"
)

## Comparing between different networks --------------
# Again, this pretty much works as you expect: your network measures are just
# another thing you could have measured empirically. Maybe you want to compare
# between a treatment and a control, or between two species.
#### Importing the zebra data again -------
# Let's use the zebra datasets. This is just the code to generate the network objects,
# copied from an earlier script:
# Make sure you change 'path' to match your computer
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Rubenstein etal zebras/"
grevys_raw <- read.table(paste(path, "grevys.csv", sep = ""), header=T, sep = ",")
wildasses_raw <- read.table(paste(path, "wildasses.csv", sep = ""), header=T, sep = ",")
grevys_indiv <- data.frame(ID = unique(grevys_raw$Individual.ID))
wildasses_indiv <- data.frame(ID = unique(wildasses_raw$Individual.ID))
for (i in 1:nrow(grevys_indiv)) {
  groups_list <- list()
  for (r in 1:nrow(grevys_raw)) {
    if(grevys_raw$Individual.ID[r]==grevys_indiv$ID[i]) {
      groups_list <- append(groups_list, grevys_raw$Group.No[r])
    }
  }
  grevys_indiv$groups[[i]] <- groups_list
}  
grevys_interactions_m <- matrix(0, nrow = nrow(grevys_indiv), ncol = nrow(grevys_indiv))
colnames(grevys_interactions_m) <- grevys_indiv$ID
rownames(grevys_interactions_m) <- grevys_indiv$ID
for (r in 1:nrow(grevys_indiv)) {
  for (c in 1:nrow(grevys_indiv)) {
    weight <- 0
    weight <- length(intersect(grevys_indiv$groups[[r]], grevys_indiv$groups[[c]]))
    grevys_interactions_m[r,c] <- weight
  }
  grevys_interactions_m[r,r] <- 0 # Remove the loops
}
grevys <- graph_from_adjacency_matrix(grevys_interactions_m, weighted=TRUE, mode = "undirected")
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
# Ok now we have two network objects, one for grevys and one for wild asses.
# Actual two-network comparison -------------
wilcox.test(closeness(grevys), closeness(wildasses))
# The p-value is 0.01, which means there is only a 1% chance of getting medians
# this or more different if the zebra species have the same median 'closeness'.
# We call this 'significant', and accept the alternative hypothesis (that 
# grevy zebras have generally lower closeness). 
par(mfrow=c(1,1), oma = c(0,0,0,0), mar = c(2,6,1,1)
    , mgp=c(4, 1, 0), las=1) # bottom, left, top, right
colors_zebraspecies <- c("sienna2", "sienna4")
y_max <- 0.20
y_axis_offset <- y_max * 0.95 # This puts sample size numbers 5% below max
Nice_Plot <- boxplot(closeness(grevys), closeness(wildasses)
                     , xlab = "Zebra species"
                     , ylab = "Closeness for each observed individual" 
                     , names = c("Grevy's", "Wild Ass")
                     , col = alpha(colors_zebraspecies, 0.5)
                     , ylim = c(0, y_max)
                     , range = 0 # I like my whiskers to go to max values
)
nbGroup <- nlevels(as.factor(Nice_Plot$names))
text( 
  x=c(1:nbGroup), 
  y=y_axis_offset,
  cex = 1,
  paste(Nice_Plot$n,sep="")  
)
stripchart(closeness(grevys)
           , add = TRUE
           , pch = 19
           , col = colors_zebraspecies[1]
           , method = "jitter"
           , jitter = 0.2
           , vertical = TRUE
)
stripchart(closeness(wildasses)
           , at = 2
           , add = TRUE
           , pch = 19
           , col = colors_zebraspecies[2]
           , method = "jitter"
           , jitter = 0.2
           , vertical = TRUE
)
# Not bad as boxplots go, but there are two things very obvious here that make it 
# not ideal.
# One is that the Grevy's Zebras have 4 points that are far away from the rest of the
# distribution. That makes us wonder if there is something wrong with what we 
# calculated, or if not, what is special about these individuals.
# Second, because of these outliers, the difference in medians is actually not well
# illustrated here: the distributions are so squished their difference is hard to see.

# One could approach these two problems different ways. One way to deal with very 
# skewed distributions is to 'transform' them, e.g. to use a logarithmic y-axis here:
boxplot(closeness(grevys), closeness(wildasses)
        , xlab = "Zebra species"
        , ylab = "Closeness for each observed individual" 
        , names = c("Grevy's", "Wild Ass")
        , col = alpha(colors_zebraspecies, 0.5)
        , range = 0 # I like my whiskers to go to max values
        , log = "y"
)
# This seems to solve some of the problem, but personally I don't like it: we
# are just not that good at interpreting logarithmic scales. E.g. did you actually
# notice that the y axis wasn't linear? And your intuition about how far away the highest
# points for Grevy's are is probably just not correct here. 

#### Plotting with an axis gap -------------------
# Another solution to that second problem is to put a gap in the axis. 
# We'll do that by effectively plotting two separate graphs on top of each other.
par(mfrow=c(2,1), oma = c(3,3,1,0), mar = c(1,3,0,1)
    , mgp=c(4, 1, 0), las=1) # bottom, left, top, right
gap_start <- 0.025
gap_end <- 0.15
y_max <- gap_start + gap_end
y_axis_offset <- y_max * 0.98
Nice_Plot1 <- boxplot(closeness(grevys), closeness(wildasses)
                      , xlab = ""
                      , ylab = "" 
                      , xaxt = "n"
                      , names = c("Grevy's", "Wild Ass")
                      , col = alpha(colors_zebraspecies, 0.5)
                      , ylim = c(gap_end, y_max)
                      , range = 0 # I like my whiskers to go to max values
)
text( 
  x=c(1:2), 
  y=y_axis_offset,
  cex = 1,
  paste(Nice_Plot1$n,sep="")  
)
stripchart(closeness(grevys)
           , add = TRUE
           , pch = 19
           , col = colors_zebraspecies[1]
           , method = "jitter"
           , jitter = 0.2
           , vertical = TRUE
)
Nice_Plot2 <- boxplot(closeness(grevys), closeness(wildasses)
                      , xlab = "Zebra species"
                      , ylab = "Closeness for each observed individual" 
                      , names = c("Grevy's", "Wild Ass")
                      , col = alpha(colors_zebraspecies, 0.5)
                      , ylim = c(0, gap_start)
                      , range = 0 # I like my whiskers to go to max values
)
stripchart(closeness(grevys)
           , add = TRUE
           , pch = 19
           , col = colors_zebraspecies[1]
           , method = "jitter"
           , jitter = 0.2
           , vertical = TRUE
)
stripchart(closeness(wildasses)
           , at = 2
           , add = TRUE
           , pch = 19
           , col = colors_zebraspecies[2]
           , method = "jitter"
           , jitter = 0.2
           , vertical = TRUE
)
par(las=0)
mtext(text = "Closeness for each observed individual", 
      side = 2, line = 1, outer = TRUE)
mtext(text = "Zebra species", 
      side = 1, line = 1, outer = TRUE)
# This is a bit jerry-rigged in that the stripchart is just completely replotted
# for the Grevy's in the top and bottom graph, but since points outside the graph
# range are ignored, it works here. 

#### Illustrating what's up with those outliers ---------
# This is very similar to how we plotted the dolphins colored by 
# their closeness and betweenness earlier. 
par(mfrow=c(1,2), mar = c(0,0,2,0))
color_closeness <- viridis(100)
mx_close <- max(closeness(wildasses), closeness(grevys), na.rm = T)
plot(grevys
     , main = "Grevy's"
     , vertex.color = color_closeness[100*closeness(grevys)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     #, vertex.size = degree(grevys)
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(wildasses
     , main = "Wild Asses"
     , vertex.color = color_closeness[100*closeness(wildasses)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
# The plot immediately illustrates why we get the outliers in the boxplot:
# closeness is only defined within connected networks, and so the closeness is
# calculated within each subgroup; and the small group is much higher than in 
# the larger groups (essentially because everyone is as close as they can be in 
# that little group of 4...).
# Whether this matters scientifically I guess comes down to what the meaning of
# these links and 'closeness' really is, but I would say here that this 
# disconnectedness actually makes comparing closeness values across the groups, 
# and the species, not meaningful. Just a single connection added from the big
# to the small subgroup in the Grevy's zebras would totally change the picture. 
# *Nonetheless*: note that here the nonparametric Mann-Whitney-U-Test came in 
# particularly handy, as it (and the group medians) are not particularly affected 
# by the 'outlying' group of 4, and the median closeness is higher in Wild Asses.

# Randomizing networks ----------------------
# Above we either treated network measures as regular empirical parameters
# to do statistics with, or we compared two (or more) networks. But 
# sometimes you might just want to know for a particular network, in what
# way is the network structure 'remarkable'? I.e. are network measures like
# closeness or degree distributions or diameter different from what we expect?
# The answer to this often depends on what we mean by 'expect' - in other
# words, what is the null hypothesis? 
# We intuitively often would like to compare to a 'random' network. But there
# are many possible versions of this, and you cannot get around thinking about
# what the correct comparison is for your purposes. 

# Let's say we want to understand more about the 'wildasses' network. 
# For the purpose of this example, I'm going to remove the single unconnected
# individual to make our lives easier. In your real network, you'll have to make
# a judgement call, and provide a better justification, if you want to remove nodes.
wildassN <- delete_vertices(wildasses, "03-005T")
color_closeness <- viridis(100)
warnings()
### Random network 1 ------------
# One way to think about it is that we just want a network with the same number
# of nodes and links. 
# In the iGraph package, generating a random graph is referred to as 'sampling
# from a network model'. For example,
newAsses1a <- sample_gnm(vcount(wildassN), ecount(wildassN)
                         , directed = FALSE, loops = FALSE)
newAsses2a <- sample_gnm(vcount(wildassN), ecount(wildassN)
                         , directed = FALSE, loops = FALSE)
newAsses3a <- sample_gnm(vcount(wildassN), ecount(wildassN)
                         , directed = FALSE, loops = FALSE)
# Now we have three new networks, all resampled from the same distribution of
# random networks with the specified number of edges and vertices.
par(mfrow=c(2,2), mar = c(1,1,2,1))
mx_close <- max(closeness(wildassN)
                , closeness(newAsses1a)
                , closeness(newAsses2a)
                , closeness(newAsses3a)
                , na.rm = T)
plot(wildassN
     , main = paste("Wild Asses, V:", vcount(wildassN), " E:", ecount(wildassN))
     , vertex.color = color_closeness[100*closeness(wildassN)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(newAsses1a
     , main = paste("V:", vcount(newAsses1a), " E:", ecount(newAsses1a))
     , vertex.color = color_closeness[100*closeness(newAsses1a)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(newAsses2a
     , main = paste("V:", vcount(newAsses2a), " E:", ecount(newAsses2a))
     , vertex.color = color_closeness[100*closeness(newAsses2a)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(newAsses3a
     , main = paste("V:", vcount(newAsses3a), " E:", ecount(newAsses3a))
     , vertex.color = color_closeness[100*closeness(newAsses3a)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
# The first thing you should notice here is that the three random networks
# look pretty similar, and not at all like the actual wildassN network:
# among other things, the closeness values seem really different. 
# This indicates that in fact the real animals don't connect 'randomly' in this
# way. 
# But as we'll see, there are different ways of interpreting what 'random' means
# (which is why mathematicians really dislike that term).

### Random network 2 ------------
# Another way to think about it would be to say that the *probability* of 
# a link existing between any two nodes should be the same. 
# That probability is the number of edges that exist divided by the number of
# edges that could exist, which is ((number of nodes)^2-number of nodes)/2 
# if we have n undirected network without self-loops, as here. 
edge_prob <- ecount(wildassN)*2/(vcount(wildassN)^2-vcount(wildassN))
newAsses1b <- sample_gnp(vcount(wildassN), edge_prob
                         , directed = FALSE, loops = FALSE)
newAsses2b <- sample_gnp(vcount(wildassN), edge_prob
                         , directed = FALSE, loops = FALSE)
newAsses3b <- sample_gnp(vcount(wildassN), edge_prob
                         , directed = FALSE, loops = FALSE)
# Now we have three new networks, all resampled from the same distribution of
# random networks with the specified number of edges and vertices.
par(mfrow=c(2,2), mar = c(1,1,2,1))
mx_close <- max(closeness(wildassN)
                , closeness(newAsses1b)
                , closeness(newAsses2b)
                , closeness(newAsses3b)
                , na.rm = T)
plot(wildassN
     , main = paste("Wild Asses, V:", vcount(wildassN), " E:", ecount(wildassN))
     , vertex.color = color_closeness[100*closeness(wildassN)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(newAsses1b
     , main = paste("V:", vcount(newAsses1b), " E:", ecount(newAsses1b))
     , vertex.color = color_closeness[100*closeness(newAsses1b)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(newAsses2b
     , main = paste("V:", vcount(newAsses2b), " E:", ecount(newAsses2b))
     , vertex.color = color_closeness[100*closeness(newAsses2b)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(newAsses3b
     , main = paste("V:", vcount(newAsses3b), " E:", ecount(newAsses3b))
     , vertex.color = color_closeness[100*closeness(newAsses3b)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
# What is different about this new way of generating random networks?

### Random network 3 -----------------
# Let's do a totally different way of getting 'random' networks:
newAsses1c <- rewire(wildassN, keeping_degseq(loops = FALSE, niter = 100))
newAsses2c <- rewire(wildassN, keeping_degseq(loops = FALSE, niter = 100))
newAsses3c <- rewire(wildassN, keeping_degseq(loops = FALSE, niter = 100))
par(mfrow=c(2,2), mar = c(1,1,2,1))
mx_close <- max(closeness(wildassN)
                , closeness(newAsses1c)
                , closeness(newAsses2c)
                , closeness(newAsses3c)
                , na.rm = T)
plot(wildassN
     , main = paste("Wild Asses, V:", vcount(wildassN), " E:", ecount(wildassN))
     , vertex.color = color_closeness[100*closeness(wildassN)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(newAsses1c
     , main = paste("V:", vcount(newAsses1c), " E:", ecount(newAsses1c))
     , vertex.color = color_closeness[100*closeness(newAsses1c)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(newAsses2c
     , main = paste("V:", vcount(newAsses2c), " E:", ecount(newAsses2c))
     , vertex.color = color_closeness[100*closeness(newAsses2c)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(newAsses3c
     , main = paste("V:", vcount(newAsses3c), " E:", ecount(newAsses3c))
     , vertex.color = color_closeness[100*closeness(newAsses3c)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
# What happened here? The rewire() function with the keeping_degseq option 
# generates a 'randomized' network differently, namely by starting from the original
# network and 'rewiring' it, i.e. moving the links around. It thus also always
# has the same number of edges and nodes as the original, but on top of that, 
# this option preserves the degree distribution across nodes. I.e. the nodes
# all keep the same number of links they had before, just now to random other
# nodes. 

#### Qualitative comparison -----------------
### Assuming we're still primarily interested in the 'closeness' measure, 
### let's see how these different network variations perform on it. 
closen_WA <- rbind(closeness(wildassN)
                   , closeness(newAsses1a)
                   , closeness(newAsses2a)
                   , closeness(newAsses3a)
                   , closeness(newAsses1b)
                   , closeness(newAsses2b)
                   , closeness(newAsses3b)
                   , closeness(newAsses1c)
                   , closeness(newAsses2c)
                   , closeness(newAsses3c)
)
closen_WA <- t(as.matrix(closen_WA))
colnames(closen_WA) <- c("WildAss", "N&L 1", "N&L 2","N&L 3"
                         ,"Prob 1","Prob 2","Prob 3"
                         ,"Rw deg 1", "Rw deg 2", "Rw deg 3")
# Now that we made a table of all the closeness values, let's plot it. 
par(mfrow=c(1,1), oma = c(0,0,0,0), mar = c(6,4,1,1)
    , mgp=c(3, 1, 0), las=2) # bottom, left, top, right
y_max <- 0.04
y_axis_offset <- 0.005 # This puts sample size numbers at the bottom
Nice_Plot <- boxplot(closen_WA
                     #, xlab = "Network"
                     , ylab = "Closeness" 
                     , col = color_closeness[100*apply(closen_WA,2,median)/max(closen_WA)]
                     , ylim = c(0, y_max)
                     , range = 0
)
nbGroup <- nlevels(as.factor(Nice_Plot$names))
text( 
  x=c(1:nbGroup), 
  y=y_axis_offset,
  cex = 0.75,
  paste(Nice_Plot$n,sep="")  
)
# What does this tell you about closeness and the different randomization
# algorithms?
# You may want to try the same thing with other measures, such as degree, or 
# betweenness. 

## Testing a measure against 'random' ------------------
# Let's remember our real problem. We have a single network (here the contact
# network for a zebra species, the wild asses), and we want to know whether it is
# somehow different from 'expectation'. 
# So we want to compare the real measure (e.g. of closeness) to that measure in 
# a random network. 
# Now that we made some random networks, we know our first decision is about which 
# random network type to compare to. Let's say we want to compare to a network that 
# preserves the degree distribution (like the last three we generated).
# You could simply generate one such random network and use a U-Test to see if the 
# median closeness was different from that of the real data.
closen_WA_df <- as.data.frame(closen_WA)
wilcox.test(closen_WA_df$WildAss, closen_WA_df$`Rw deg 1`)
# This is highly significant, not too surprising from the plot.
# The p value is given as 8.1e-10, which means
# p = 0.00000000081 .

# However, this was one arbitarily picked randomized network, and since the computer 
# can do lots of these, one would probably normally generate ~50-100 such randomized
# networks to generate a smooth expected distribution of values of interest. 
# Let's do it:
# (Before you generate a whole bunch of new data: remember that your computer
# has limited working memory. Use the broom occasionally to delete all data, 
# variables and graphs you no longer need. 
# To run this code, the only thing you do need is the wild asses data. Here we import
# it again, assuming you deleted it:)


##### Reimporting wildasses in case you cleared your environment -------------------
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
color_closeness <- viridis(100)
### Generating 50 random rewirings -------------------
n <- 100
big_closeness_dist <- data.frame(NULL)
for (i in 1:n) {
  random_net <- rewire(wildassN, keeping_degseq(loops = FALSE, niter = 100))
  big_closeness_dist <- rbind(big_closeness_dist, closeness(random_net))
}
# A good way to visualize a complex distribution is to use a density plot:
plot(density(apply(big_closeness_dist,1,median)))
## Quick note on this command: apply() with 1 as the second argument applies
## the function to each row; a 2 as the second argument would apply the function
## to each column. 

# However, on closer inspection, the medians always take one of two values:
hist(apply(big_closeness_dist,1,median))
# Interesting! It seems the degree distribution (which we have fixed) might
# have a big influence on 'closeness' across the network; the fact that there
# are not many possible values may have something to do with the limited number
# of nodes, as well. However, this is all speculation.  

# What we actually wanted was to compare the empirical median closeness to the 
# distribution we generated. While we could use a U-Test to compare simply the 
# medians of the two distributions - 
wilcox.test(as.matrix(big_closeness_dist), closeness(wildassN))
# Unsurprisingly, this has an even lower p-value than just comparing to one randomized
# network. 
# This is however not the strength of having generated a large number of randomzied
# networks to compare to. The strength of this approach is that we can design our
# own null hypothesis tests. What a null hypothesis test does is it tells us:
# is this value more extreme than 95% of values generated by the null hypothesis?
# We can now see directly what those values are. For example, 

# Is the median of the empirical data lower than 95% of randomized medians?
rand_medians_05 <- quantile(apply(big_closeness_dist,1,median), probs = c(0.05))
ifelse(median(closeness(wildassN))<rand_medians_05, "YES, signficantly smaller", "NO, not significantly smaller")
rand_medians_95 <- quantile(apply(big_closeness_dist,1,median), probs = c(0.95))
ifelse(median(closeness(wildassN))>rand_medians_95, "YES, signficantly bigger", "NO, not significantly bigger")
# (Note that these are the equivalent of one-tailed tests...)
# But we can now do this for a number of other measures. How about the variation in  
# closeness among nodes? (I'll use interquartile range as a measure of variation)
rand_IQR_05 <- quantile(apply(big_closeness_dist,1,IQR), probs = c(0.05))
ifelse(IQR(closeness(wildassN))<rand_IQR_05, "YES, signficantly smaller", "NO, not significantly smaller")
rand_IQR_95 <- quantile(apply(big_closeness_dist,1,IQR), probs = c(0.95))
ifelse(IQR(closeness(wildassN))>rand_IQR_95, "YES, signficantly bigger", "NO, not significantly bigger")

# How about checking the range of closeness values achieved for each node, i.e. 
# for a given degree?
y_max <- 0.035
plot(closeness(wildassN)~degree(wildassN)
     , col = "black"
     , pch = 18
     , ylim = c(0.01, y_max)
)
big_closeness_dist_T <- as.data.frame(t(as.matrix(big_closeness_dist)))
row.names(big_closeness_dist_T) <- V(wildassN)$name
mx_close <- y_max
for (i in 1:n) {
  points(big_closeness_dist_T[,i] ~ jitter(degree(wildassN))
         , pch = 19
         , col = color_closeness[100*big_closeness_dist_T[,i]/mx_close]
  )
}
# Ok what does this tell us?
# The black diamonds are the empirical data: wild asses, their individual degree
# (=number of connections in the network) and closeness (=a measure of centrality
# in the network, i.e. dependent on global structure).
# The colored points/distributions are the results of 100 random reshufflings of 
# all the connections in the network, while keeping the degree of each individual
# animal constant (there is just a bit of x-axis jitter added for point visibility;
# if you wanted you could make them into violin plots instead). 
# The real zebras all have lower closeness than expected; but also, their closeness
# is more even than expected from a random network even with the same degree distribution.
# The difference between expected and empirical is larger for more highly connected 
# individuals. 
# This could for example result if the low-degree individuals are more likely to 
# be connected to others that are highly connected; although whatever our mechanistic
# or even phenomenological explanation is here would have to be tested.


