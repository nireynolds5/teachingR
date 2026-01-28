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

# Importing data -----------------------
data(emon)
emon1 <- emon[[1]]
organizations <- asIgraph(emon1)
V(organizations)$name <- V(organizations)$vertex.names

# On D2L, I'm also giving you a dataset about interactions between dolphins. You
# will have to download it and put it into your working directory before running this.

# Make sure you change the path variable to match your computer, or make it empty
# path <- ""
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Lusseau dolphins/"
dolphins_df <- read.table(paste(path, "dolphinnetwork.csv", sep = ""), sep = ",")
dolphins_nodeinfo <- read.table(paste(path, "dolphinvertexinfo.csv", sep = ""), header = TRUE, sep = ",")
dolphins <- graph_from_adjacency_matrix(as.matrix(dolphins_df))
V(dolphins)$sex <- dolphins_nodeinfo$dolphin.sex
V(dolphins)$disappeared <- dolphins_nodeinfo$dolphin.disappeared

# Motifs ------------------------------
# Motifs are sort of 'micro-networks', ways in which small sets of nodes can be
# connected. Some researchers argue that these are informative about the algorithmic
# /computational building blocks of the system.

## This function gives an array which contains the number of times each possible 
## motif of size 3 occurs.
motifs(organizations, size=3) # Remember 'organizations' is one of our iGraph objects

# We can easily use 'barplot' to show us the number of motifs, 
# but it's not terribly helpful if we don't know what motif 
# each bar corresponds to. 
# Here is a code that generates a composite graph not only of the 
# barplot but also of what the motifs look like.
# First we set the layout: we need to plot one big graph and 16
# little graphs, which are all the possible motifs of size 3.
par(mar=c(0,1,0,1), oma = c(1,2,1,1), xpd=TRUE)# bottom, left, top, right
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17, 0)
              , nrow = 2, ncol = 18, byrow = T)
       , heights=c(4,1)) 
color_distributions <- "plum"
# Here is the barplot - 'organizations' is the network I am 
# plotting this for, but you can replace that with any of our iGraph
# objects. Try out UScities and dolphins, for example - did it look
# like what you expected?
barplot(motifs(organizations, size=3)
        , col = color_distributions
        , names.arg = seq(1, 16)
)
# Now, under the barplot (in the second row of the 'layout'),
# we plot 16 little graphs.
par(mar=c(0,0.5,0,0))
for (i in 0:15) {
  # This command gives the graph number i that is possible with
  # 3 nodes. 
  motif_graph <- graph_from_isomorphism_class(3, i)
  # Now plotting that:
  plot(motif_graph
       , edge.arrow.size = 0.5
       #, edge.arrow.width = 1
       , edge.color = alpha("grey27", 0.5)
       , edge.width = 2
       , vertex.label.family = "Arial"
       , vertex.label.color = color_distributions
       , vertex.label.cex = 1
  )
}
# Depending on the size of your plot window, you may have to play around
# a bit with the sizes of the arrows to make them properly visible.


# Community detection ----------------------
# A lot of research in network analysis has gone into 'community detection', 
# i.e. trying to identify 'modules' int he network that are not entirely disconnected
# from the rest but somehow still distinct. 
# There are many algorithms for finding/defining such 'communities' or 'modules',
# so if you are doing this for a research purpose make sure you find out which one
# best suits your purposes. Several functions are 
# available in the igraph package: https://r.igraph.org/reference/index.html#community-detection
# For example:
communities_dolphins <- cluster_edge_betweenness(dolphins)
# This takes an igraph object (here: dolphins) as argument, and produces a 'communities'
# object (defined by the igraph package). There are several functions that deal with 
# such communities objects. 

# Let's start by visualizing what this does.
# Here I want to compare at least two community detection functions, so we'll first
# generate a layout that we can then repeatably use. 
dolphins_coordinates <- layout_nicely(dolphins)
# Now we'll set up the layout for the graphs:
par(mfrow=c(1,2)) # plot two graphs next to each other
par(mar = c(0,0,1,0)) # only have a top margin # bottom, left, top, right
# This is the actual community detection. 
communities_dolphins <- cluster_edge_betweenness(dolphins)
# Now we'll define a color list for the detected communities
colors_groups <- viridis(length(communities_dolphins))
# And now plot the first one
plot(dolphins
     , vertex.color = colors_groups[membership(communities_dolphins)]
     , vertex.size = 8
     , vertex.label.family = "Arial"
     , vertex.label.color = "black"
     , vertex.label.cex = 0.5
     , vertex.label.dist = 1
     , edge.arrow.size = 0.2
     , edge.arrow.width = 1
     , edge.color = "grey27"
     , edge.width = 2
     , mark.groups = communities_dolphins
     , mark.col = alpha(colors_groups, 0.3)
     , mark.border = NA
     , layout = dolphins_coordinates #we've predefined the arrangement of nodes
     # so that we can reproduce it
)
# Now a second graph with a different community detection algorithm:
communities_dolphins <- cluster_leading_eigen(dolphins)
# This one may detect a different number of communities, so we update
# the color vector
colors_groups <- viridis(length(communities_dolphins))
# And now plot again
plot(dolphins
     , vertex.color = colors_groups[membership(communities_dolphins)]
     , vertex.size = 8
     , vertex.label.family = "Arial"
     , vertex.label.color = "black"
     , vertex.label.cex = 0.5
     , vertex.label.dist = 1
     , edge.arrow.size = 0.2
     , edge.arrow.width = 1
     , edge.color = "grey27"
     , edge.width = 2
     , mark.groups = communities_dolphins
     , mark.col = alpha(colors_groups, 0.3)
     , mark.border = NA
     , layout = dolphins_coordinates
)
