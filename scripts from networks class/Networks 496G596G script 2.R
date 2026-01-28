# Header-------------------------
## R script for ECOL496G 596G 
## Complex Systems and Networks
## UA class Spring 2025
## written by Anna Dornhaus

## This copy belongs to:

# Prep --------------------------
# graph objects
# color
library(viridis)
library(scales)

library(sna)
library(intergraph)
library(igraph)


# Importing data -----------------------

# UScities distances
dist_matrix <- as.matrix(UScitiesD)
dist_matrix[dist_matrix>0] <- 1
UScities <- graph_from_adjacency_matrix(dist_matrix)
dist_list <- as.list(dist_matrix)
omit_zeros_list <- dist_list[dist_list != 0]
dist_to_weight <- function(x) {return(2000/x)}
E(UScities)$distance <- lapply(omit_zeros_list, dist_to_weight)

# Dolphins
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Lusseau dolphins"
dolphins_df <- read.table(paste(path, "/dolphinnetwork.csv", sep = ""), sep = ",")
dolphins_nodeinfo <- read.table(paste(path, "/dolphinvertexinfo.csv", sep = ""), header = TRUE, sep = ",")
dolphins <- graph_from_adjacency_matrix(as.matrix(dolphins_df))
V(dolphins)$sex <- dolphins_nodeinfo$dolphin.sex
V(dolphins)$disappeared <- dolphins_nodeinfo$dolphin.disappeared
dolphins_sna <- asNetwork(dolphins)

# Pollination network
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Memmott pollinator network"
pollination_df <- read.table(paste(path, "/memmott_1999.csv", sep = ""), header = TRUE, sep = ",")
pollination_df <- pollination_df[-1]
pollination_nodeinfo1 <- read.table(paste(path, "/plants.csv", sep = ""), sep = ",")
pollination_nodeinfo2 <- read.table(paste(path, "/pollinators.csv", sep = ""), header = TRUE, sep = ",")
colnames(pollination_df) <- pollination_nodeinfo1$V2
rownames(pollination_df) <- pollination_nodeinfo2$species
pollination <- graph_from_biadjacency_matrix(as.matrix(pollination_df), directed = FALSE, weighted = TRUE)

# Others
data(emon)
data(flo)
data(coleman)
florentines <- graph_from_adjacency_matrix(flo)
emon1 <- emon[[1]]
organizations <- asIgraph(emon1)
V(organizations)$name <- V(organizations)$vertex.names

# Calculating basic network measures -------------------

# Just in case you had loaded this previously - otherwise ignore
#detach("package:circlize", unload=TRUE)
library(igraph)

## Importantly, both SNA and iGraph have many network analysis measures built in. 
## Here I am using mostly the iGraph versions since that is the package we loaded last,
## and identical functions names from SNA were overwritten.

## I'm giving you examples of commonly calculated network measures. Experiment a bit!

plot(dolphins)
degrees <- sort(degree(dolphins))
barplot(degrees, ylab="Degree")
barplot(table(degrees), ylab="Frequency", xlab="Degree of individual dolphins", col="blue")
barplot(degree_distribution(dolphins), ylab="Frequency", xlab="Degree of individual dolphins", col="blue")
boxplot(closeness(dolphins))
plot(betweenness(dolphins)~closeness(dolphins))
vcount(dolphins)
ecount(dolphins)
graph.density(dolphins)
reciprocity(dolphins)
dyad_census(dolphins)
transitivity(dolphins)  ## same as clustering coefficient
mean(degree(dolphins, mode = "in")) 

# Better network visualization -------------------

## Here are four ways to visualize the same interaction matrix:
UScitiesD
View(as.matrix(UScitiesD))
plot(UScities, edge.width = E(UScities)$how_close, edge.arrow.width = 0)
heatmap(as.matrix(UScitiesD), col = magma(2000), margins=c(8,8))



# We'll go through these as well as customizing. It is important not to get carried
# away in graphics, but to focus on what you want to illustrate for the 
# scientific point you want to make. 

# Often, color, arrangement, ordering, and subsetting are essential to illustrate
# the important properties we care about. 

# Remember that the relevant claims/results in your paper need to be supported by
# the correct statistical tests or comparisons; the figures in your paper should
# serve to *illustrate* the results. Illustrating means convincingly and honestly
# showing in an easily digestible/understandable way, and that includes not just
# significant comparisons (maybe those the least), but also descriptive/exploratory
# data and particularly non-significant comparisons are important to show 
# graphically. 
# The size of your dataset, as well as the variation, has a big influence on whether
# it is better to focus on summary information (e.g. a boxplot) or to include all
# individual data points (whether that's a raincloud plot or similar or a full 
# network plot). 

# An example:
MaleCloseness <- subset(closeness(dolphins), V(dolphins)$sex == "M")
FemaleCloseness <- subset(closeness(dolphins), V(dolphins)$sex == "F")
MaleBetweenness <- subset(betweenness(dolphins), V(dolphins)$sex == "M")
FemaleBetweenness <- subset(betweenness(dolphins), V(dolphins)$sex == "F")

plot(MaleBetweenness ~ MaleCloseness, col="blue")
## 'points()' adds (scatterplot) points to an existing graph, but otherwise has essentially
## the same syntax as 'plot()'.
points(FemaleBetweenness ~ FemaleCloseness, col="orange")

## What's missing here to make this a great visualization?

# Set a color array to use throughout your paper
colors_dolphinsex <- viridis(3) # This is a good idea to do early on in your script
# But if you just need two or three colors, you may want to pick them individually rather
# than from a color scale. E.g.
colors_dolphinsex <- c("darkseagreen", "darkmagenta", "orange")
color_distributions <- "plum"
# 'par()' can be used to set inner and outer margins - https://r-graph-gallery.com/74-margin-and-oma-cheatsheet.html
# I always leave a note in my script: bottom, left, top, right
par(mar = c(5,5,1,1) # margin of 5 on left and bottom, 1 otherwise
    , oma = c(0,0,0,0) # no outer margins
)
plot(betweenness(dolphins) ~ closeness(dolphins) # we could have stuck with the 
     # subsets above, but this is a little cleaner - note that we still assign
     # colors by the variable $sex in the node list
     , pch = 19 # set point shape
     , col = colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))]
     # This is a little longer than it would be if we just had used the subsets of data. I like it
     # better because it keeps flexibly using the entire dataset. 
     , cex = 1.5 # point size - 1 is default
     , xlab = "Dolphin Closeness Measure"
     , ylab = "Dolphin Betweenness Measure"
)
label_xoffset <- -max(closeness(dolphins)) * 0.002 # you have to play around with 
# this to see what looks good. I do it relative to the x-axis for comparability
# between plots. 
label_yoffset <- 0
text(closeness(dolphins) + label_xoffset # x coordinates of labels
     , betweenness(dolphins) + label_yoffset # y coordinates of labels
     , labels = V(dolphins)$name # text in labels
     , cex = 0.5 # size of text
     , col = colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))] # colors as before
     , pos = 4 # make the text left aligned (to the right of given coordinates)
)
legend("topleft"
       , c("Male", "Female", "Unknown")
       , col = colors_dolphinsex
       , pch = 19 # you'll normally match the shape of the scatterplot points
)
# The plot now gives us a good intuition about how closeness relates to betweenness
# in this dolphin social network. (I am just making up that you want this - for
# your own analysis, think about what you actually want.)
# In addition though, it gives a wealth of descriptive data - the spread in either measure, 
# how dolphin sex affects it, even information about individual dolphins. 

# We can even upgrade this to a multi-panel plot. For this we can use another par()
# setting (e.g. mfrow=c(2,2)), or for more control use 'layout()'.
layout(matrix(c(1,2,0,3), 2, 2, byrow = T), widths=c(1,5), 
       heights=c(5,1)) 
# This gives us a four-panel plot; the plots will be inserted into the panels in 
# order by row.
# Various other format adjustments
par(oma = c(0,0,0,0), mgp=c(3, 1, 0), las=1)
# Panel 1:
par(mar = c(4,0,1,0)) # bottom, left, top, right
boxplot(betweenness(dolphins)
        , xaxt = 'n'
        , yaxt = 'n'
        , frame = FALSE
        , range = 0
        , col = color_distributions
)
# Panel 2: The original graph
par(mar = c(4,4,1,2)) # bottom, left, top, right
plot(betweenness(dolphins) ~ closeness(dolphins) # we could have stuck with the 
     # subsets above, but this is a little cleaner - note that we still assign
     # colors by the variable $sex in the node list
     , pch = 19 # set point shape
     , col = colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))]
     # This is a little longer than it would be if we just had used the subsets of data. I like it
     # better because it keeps flexibly using the entire dataset. 
     , cex = 1.5 # point size - 1 is default
     , xlab = "Dolphin Closeness Measure"
     , ylab = "Dolphin Betweenness Measure"
)
label_xoffset <- -max(closeness(dolphins)) * 0.002 # you have to play around with 
# this to see what looks good. I do it relative to the x-axis for comparability
# between plots. 
label_yoffset <- 0
text(closeness(dolphins) + label_xoffset # x coordinates of labels
     , betweenness(dolphins) + label_yoffset # y coordinates of labels
     , labels = V(dolphins)$name # text in labels
     , cex = 0.5 # size of text
     , col = colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))] # colors as before
     , pos = 4 # make the text left aligned (to the right of given coordinates)
)
legend("topleft"
       , c("Male", "Female", "Unknown")
       , col = colors_dolphinsex
       , pch = 19 # you'll normally match the shape of the scatterplot points
)
# Panel 3: empty
# We put a 0 in the layout matrix there, so R should know we don't want this to be used
# Panel 4: boxplot of x-axis distribution
par(mar = c(0,4,0,2)) # bottom, left, top, right
boxplot(closeness(dolphins)
        , xaxt = 'n'
        #, yaxt = 'n'
        , frame = FALSE
        , range = 0
        , col = color_distributions
        , horizontal=TRUE
)

## A different example. 
## Remember that each new graph will replace the last one, but you can use the back 
## arrow to cycle though all previous graphs of this session.
## You can export graphs by clicking on the graph, then in the Export' drop-down
## menu you can save it or by copy it to the clipboard. Note that the 'graph', not
## the image, is in the clipboard; when pasting into PowerPoint, for example, you 
## need to choose 'paste as image' in the Paste Special menu. 

# Resetting par()
par(oma = c(0,0,0,0), mar = c(3,5,1,1), mgp=c(4, 1, 0), las=1) # bottom, left, top, right
par(mfrow=c(1,1))
# Default boxplot, what we want to illustrate is the difference between males and females:
boxplot(MaleCloseness, FemaleCloseness)

# How to make a great boxplot -
y_max <- 0.010
# Saving the plot into a variable allows us to access plot parameters afterwards.
Nice_Plot <- boxplot(closeness(dolphins) ~ as.factor(V(dolphins)$sex) # again a little cleaner this way
                     #, xlab="Dolphin sex"  #optional, but here I find it takes too much extra space
                     , ylab="Closeness [calculated over entire network regardless of sex]" 
                     # Always make axis descriptions as clear and comprehensive as possible
                     , names = c("Males", "Females", "Unknown")
                     , col = alpha(colors_dolphinsex, 0.5) # use same colors as elsewhere, 
                     # but slightly transparent so we can see data points
                     , ylim = c(0, y_max) # always think about the scale - starting from zero is typically better
)

# Putting sample sizes above bars
y_axis_offset <- Nice_Plot$stats[3,] * 0.96  # This puts sample size numbers 4% below median
y_axis_offset <- y_max * 0.95 # This puts sample size numbers 5% below max
nbGroup <- nlevels(as.factor(Nice_Plot$names))
text( 
  x=c(1:nbGroup), 
  y=y_axis_offset,
  cex = 1,
  paste(Nice_Plot$n,sep="")  
)

# I like representing the raw data as well, especially for mid- to low sample sizes.
stripchart(closeness(dolphins) ~ as.factor(V(dolphins)$sex)
           , add = TRUE
           , pch = 19
           , col = colors_dolphinsex
           , method = "jitter"
           , jitter = 0.2
           , vertical = TRUE,
)


# Now let's actually plot the network as well:
# (with iGraph package)

# Note that some layouts are re-determined randomly each time, and some are not.
# For some purposes, it may be better to have a regular spatial arrangement;
# and we may want a specific order of nodes.
V(dolphins)$position <- paste(V(dolphins)$sex, ifelse(degree(dolphins)<10, paste("0", degree(dolphins), sep = ""), degree(dolphins)))
sorted_nodes <- order(V(dolphins)$position)
plot(dolphins
     , vertex.color = colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , vertex.size = degree(dolphins, mode = "all")
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = layout_in_circle(dolphins, order = sorted_nodes)
     , vertex.label = paste(V(dolphins)$name, V(dolphins)$position)
)

# But we sometimes want the power of the algorithms to draw the nodes in a pattern, e.g.
dolphins_coordinates <- layout_nicely(dolphins)
plot(dolphins
     , vertex.color = colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , vertex.size = degree(dolphins)
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = dolphins_coordinates
)
# The above allows you to now do several graphs with the nodes in the same coordinates.
# E.g. 
par(mfrow=c(1,2), mar = c(0,0,2,0))
color_closeness <- viridis(100)
plot(dolphins
     , main = "Closeness"
     , vertex.color = color_closeness[100*closeness(dolphins)/max(closeness(dolphins))]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     #, vertex.size = degree(dolphins)
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = dolphins_coordinates
)
plot(dolphins
     , main = "Betweenness"
     , vertex.color = color_closeness[100*betweenness(dolphins)/max(betweenness(dolphins))]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     #, vertex.size = degree(dolphins)
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = dolphins_coordinates
)

# However if you run
dolphins_coordinates <- layout_nicely(dolphins)
# again you get a different set, as there is stochasticity built into the process.


# Bipartite or two-mode networks -----------------

## If your table is a two-mode network, i.e. rows and columns are different things,
## such as plants and their pollinators, the typical network plot doesn't work as well.
V(pollination)$color <- ifelse(V(pollination)$type, "darkgreen", "lightblue")
plot(pollination)
## Instead, you can use a 'bipartite' built-in layout
plot(pollination, layout = layout.bipartite)
# Still not great, but the two layers show the two types of nodes, making interactions
# across the layers visible. 

## The above uses the iGraph plotting function. The SNA package does not provide
## good ways of plotting or analyzing bipartite networks. We can also use
## the package 'bipartite', which is based on sna:
library(bipartite)
# https://cran.r-project.org/web/packages/bipartite/vignettes/Intro2bipartite.pdf
par(xpd=T)
## Note this package also has a built-in pollination network:
plotweb(motten1982)
visweb(motten1982)
# So these two plots are the classic 'bipartite network' plot and the equivalent as a
# heatmap. 

## But here is this plot for our network from the Memmott data:
plotweb(as.matrix(pollination_df))
# This already looks fairly similar to the published graph in the Memmott paper.
visweb(as.matrix(pollination_df))
# Nonetheless, I think this is a good example where the heatmap may actually be
# more informative than the 'network' plot. 

## Sankey diagrams
# Sankey diagrams are a specific type of plot that is intended for illustrating
# 'flows' from a set of categories into another set. 
library(networkD3)
pollination_nD3 <- igraph_to_networkD3(pollination)
Links <- pollination_nD3$links
Nodes <- pollination_nD3$nodes
sankeyNetwork(Links = Links, Nodes = Nodes, Source = 'source', Target = 'target', 
              NodeID = 'name', Value = 'value')
# As you can see they are not too different from bipartite network plots (just 
# typically from left to right instead of top to bottom). We could design this better
# but I still like the heatmap more. 


# Better bipartite plots -------------------

# Let's start by just making it legible and adding a little color.
par(oma=c(3,0,0,0), mar = c(0,0,0,0), xpd=TRUE)# bottom, left, top, right
# I get the impression that all inner margins are ignored.

pollination_num_matrix <- t(as.matrix(pollination_df)) # transposed because I want
# plants at the bottom

plotweb(pollination_num_matrix
        , text.rot = 90
        , col.high = "orange"
        , col.low = "darkseagreen"
)
# Or the heatmap version: we now already know how to do this.
intensity_color <- viridis(max(pollination_num_matrix))
poll_map <- heatmap(pollination_num_matrix
                    , col = intensity_color
                    , margins=c(8,8)
                    , scale = "none"
                    , cexRow = 0.5
                    , cexCol = 0.5
                    , revC = TRUE
)
# Already we can draw the same conclusions as the paper, although note that the public
# dataset does not provide data on respective abundances of these species; the width of bars
# in the bipartite plot is just the degree (number of connections).

pollinator_colors <- c("green4", "skyblue", "yellow2", "slateblue")
plant_colors <- "yellowgreen"
# First we make a matrix for the link colors
link_colors <- matrix(nrow=nrow(pollination_num_matrix), ncol=ncol(pollination_num_matrix), dimnames=dimnames(pollination_num_matrix))
alt_link_cols <- matrix(nrow=nrow(pollination_num_matrix), ncol=ncol(pollination_num_matrix), dimnames=dimnames(pollination_num_matrix))
intensity_color <- viridis(max(pollination_num_matrix))
for(i in 1:nrow(link_colors)) {
  for(j in 1:ncol(link_colors)) {
    link_colors[i,j] <-  pollinator_colors[as.factor(pollination_nodeinfo2$pol_order)[j]]
    alt_link_cols[i,j] <- intensity_color[pollination_num_matrix[i,j]+1]
  }
}
link_colors[is.na(link_colors)]="grey"

plant_order <- rownames(pollination_num_matrix)[order(rowSums(pollination_num_matrix), decreasing = FALSE)] 
poll_order <- colnames(pollination_num_matrix)
plotweb(pollination_num_matrix
        , method = "normal"
        , sequence = list(seq.high=poll_order, seq.low=plant_order)
        , text.rot = 90
        , col.high = pollinator_colors[as.factor(pollination_nodeinfo2$pol_order)]
        , col.low = plant_colors
        , col.interaction = t(link_colors)
        , bor.col.interaction	= t(link_colors)
        , text.high.col = "white" # the pollinator labels are uninformative
        , empty = FALSE
        #, plot.axes = TRUE # I used this temporarily to figure out where the top labels should go
)
# I'd like to label the top boxes with something, but this by-hand-spacing is obviously
# neither replicable (if you resize the window) nor particularly transparent or robust. 
# If anyone can find a way to access the coordinates of the boxes produced by plotweb()
# let me know. 
text(x=-0.1, y=1.7, "Coleoptera                                Diptera                          Hymenoptera   Lepidoptera", col = "black", pos = 4)

# Or we can just color matching the heatmap and the intensity of interactions:
intensity_color2 <- viridis(max(max(rowSums(pollination_num_matrix), max(colSums(pollination_num_matrix)))))
plotweb(pollination_num_matrix
        , method = "normal"
        , sequence = list(seq.high=poll_order, seq.low=plant_order)
        , text.rot = 90
        , col.high = intensity_color2[colSums(pollination_num_matrix)]
        , col.low = intensity_color2[rowSums(pollination_num_matrix)]
        , col.interaction = t(alt_link_cols)
        , bor.col.interaction	= t(alt_link_cols)
        , text.high.col = "white" # the pollinator labels are uninformative
        , empty = FALSE
        #, plot.axes = TRUE # I used this temporarily to figure out where the top labels should go
)
text(x=-0.1, y=1.7, "Coleoptera                                Diptera                          Hymenoptera   Lepidoptera", col = "black", pos = 4)


