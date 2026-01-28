# Header-------------------------
## R script for ECOL496G 596G 
## Complex Systems and Networks
## UA class Spring 2025
## written by Anna Dornhaus

## This copy belongs to:

# Prep --------------------------
# You may not need all these libraries, and I strongly suggest adding notes
# about what you need each for. 
library(viridis) # Good colorblind-proof color set
library(scales) # Easily make color scales

library(sna) # Social network analysis, which I am using for .... ?
library(intergraph) # Convert sna objects into igraph ones and vice versa, and matrix
# conversions
library(igraph) # Network analysis, which I am using for... ?

library(bipartite) # Handling and plotting bipartite networks 
library(networkD3) # Making Sankey diagrams

# Importing data -----------------------

# You really only need to import the data you are planning to use. 

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

# Graphing prep for consistency -------------------

# Set color arrays to use throughout your paper
color_closeness <- viridis(100)
colors_dolphinsex <- c("darkseagreen", "darkmagenta", "orange")
color_distributions <- "plum"

# Margins and other formats may need to be set individually for different figure
# types, but if you are doing for example a series of boxplots, you might want to
# set consistent par() parameters in the beginninng. 

# 'par()' can be used to set inner and outer margins - https://r-graph-gallery.com/74-margin-and-oma-cheatsheet.html
# I always leave a note in my script: bottom, left, top, right
par(mar = c(5,5,1,1) # margin of 5 on left and bottom, 1 otherwise
    , oma = c(0,0,0,0) # no outer margins
    , mgp=c(3, 1, 0) # distance of labels
    , las=1  # orientation of labels
    , xpd=T # allow plotting outside inner margins
)
# You may or may not also want to define some other parameters for your own use,
# e.g. to plot sample sizes in a consistent position, or labels in a consistent 
# place near nodes, etc.
label_xoffset <- 0.002 
label_yoffset <- 0
y_max <- 0.010
# And, you may or may not define a sorting or positioning early on. We won't do it 
# here since we'll be working on different spatial ordering.
#sorted_nodes <- order(V(dolphins)$position)
#dolphins_coordinates <- layout_nicely(dolphins)

# Totally different network diagrams ---------------

# Better heatmaps ------------------------------
# resetting par values
par(oma = c(0,0,0,0), mar = c(1,1,1,1), mgp=c(3, 1, 0), las=1) # bottom, left, top, right
par(mfrow=c(1,1))
# We can add additional information in the built-in color sidebars, for example
# the sex and degree of dolphins. The dendrograms and order of dolphins is identical
# on the x and y axis. 
dolphin_num_matrix <- as.matrix(dolphins_df)
d_map <- heatmap(dolphin_num_matrix
                 , col = c("slateblue4", "slateblue1")
                 , margins=c(8,8)
                 , scale = "none"
                 , cexRow = 0.5
                 , cexCol = 0.5
                 #, Rowv = TRUE
                 #, Colv = TRUE
                 , revC = TRUE
                 , RowSideColors = magma(24)[degree(dolphins)]
                 , ColSideColors = colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))]
)
# This is an undirected, and thus symmetrical matrix. It is also unweighted, so 
# all squares of the heatmap are either 0 or 1. Which makes the most interesting part
# of this graph the clustering algorithm, which grouped the dolphins by who they are 
# most similar to in their interaction pattern. 
# Note that you may have to make the plotting window large to see all dolphin names
# - the names on the x and y axis should be identical. 
# What I don't like: no built-in continuous legend function (not needed in this
# example perhaps); and no possibility of suppressing the dendrogram while keeping
# the clustering (unless one does the clustering manually).




# Chord diagrams ----------------------------------
# https://jokergoo.github.io/circlize/
# https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html
# install.packages("circlize")
# This package gives a lot of functionality for plotting anything on a circle. 
# Careful: it overwrites the function degree() to mean something completely different.
library(circlize)
chordDiagram(as.matrix(flo))
# Here the nodes are sections on a ring, with their length = network degree
# (in + out).
# The edges are plotted as arcs among nodes (color = origin node). 
# Remember that the argument has to be just the interaction matrix (not a variable
# in igraph format).
chordDiagram(dolphin_num_matrix)
# Coloring nodes
colorscale_dolphins <- viridis(1+max(igraph::degree(dolphins, mode = "total")))
V(dolphins)$color <- colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))]
# Coloring the edges by some attribute of the origin node
colormatrix <- matrix(, nrow = vcount(dolphins), ncol = vcount(dolphins))
for (i in 1:vcount(dolphins)) {
  colormatrix[i,] <- colorscale_dolphins[1+igraph::degree(dolphins, mode = "total")[i]] 
}
# Reordering nodes around the circle
#dolphin_sequence <- levels(reorder(V(dolphins)$name, igraph::degree(dolphins, mode = "total")))
dolphin_sequence <- V(dolphins)$name[order(V(dolphins)$sex)]
color_sequence <- V(dolphins)[dolphin_sequence]$color
# 'order' expects the same as union(rownames(dolphin_num_matrix), colnames(dolphin_num_matrix)) 
# just in the order you want it to be. 
# Similar to par() for normal plots, circos.par() can set various graphical 
# parameters. Unlike par(), it really wants to be reset with circos.clear() after.
circos.par(start.degree = 90 # This is starting at 'up' and going clockwise
           , canvas.xlim = c(-1.5,1.5)) 
chordDiagram(dolphin_num_matrix
             , order = dolphin_sequence
             , grid.col = color_sequence #needs to be length nodelist
             #, col = colormatrix #needs to be length of nodesxnodes (i.e. whole matrix)
             #, symmetric = TRUE # For undirected graphs, so each link won't be plotted twice
             , annotationTrack = c("grid")
             , annotationTrackHeight = c(0.1)
)
# Customizing the labeling in this function is somewhat laborious; you 
# have to separately format the 'track' (label circle)
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  #circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
  #            facing = "bending.inside", niceFacing = TRUE, col = "white", cex = 0.5)
  circos.text(mean(xlim) # 'x-axis' is around the circle
              , ylim[1] + 1.2 # 'y-axis' is radial, i.e. away from circle
              , si, sector.index = si
              , track.index = 1 # 'track' is the concentric ring in which labels are
              , facing = "clockwise" # orientation of text
              , niceFacing = TRUE
              , col = "black"
              , cex = 0.8 # text size
              , adj = c(0, 0.5)
  )
}
circos.clear() # required after changing circos.par()

## IMHO, there are too many links for this to be a helpful visualization, plus
## this visualization is particularly suited for weighted, undirected links.
## How about the zebras:

path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Rubenstein etal zebras/"
# The zebra data .csv file is a terrible format for our purposes, so we have to 
# work a little to extract the info we want.
# First, in Excel, I split it into two files:
grevys_raw <- read.table(paste(path, "grevys.csv", sep = ""), header=T, sep = ",")
wildasses_raw <- read.table(paste(path, "wildasses.csv", sep = ""), header=T, sep = ",")
# Now, each 'group' is interpreted as an interaction among all members.
# Edges are undirected, and in the resulting network an edge between two individuals
# is weighted by how often these two individuals occurred together in a group.
# We can ignore 'day' because group numbers are not repeated. 
# First, I am going to generate a dataframe that lists all individuals.
grevys_indiv <- data.frame(ID = unique(grevys_raw$Individual.ID))
wildasses_indiv <- data.frame(ID = unique(wildasses_raw$Individual.ID))
# Now we add a column with lists of all the groups that individual was in.
for (i in 1:nrow(grevys_indiv)) {
  groups_list <- list()
  for (r in 1:nrow(grevys_raw)) {
    if(grevys_raw$Individual.ID[r]==grevys_indiv$ID[i]) {
      groups_list <- append(groups_list, grevys_raw$Group.No[r])
    }
  }
  grevys_indiv$groups[[i]] <- groups_list
}  
# Now we can make an adjacency matrix, by counting the number of groups shared
# between any two individuals.
# We initialize with zeros.
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

grevys <- graph_from_adjacency_matrix(grevys_interactions_m, weighted=TRUE)
plot(grevys) # Not a nice plot, just to test that this worked
is_weighted(grevys)

# For chordDiagram, there is an additional problem: individuals with few links are
# automatically removed. If we don't want this to create problems, we should remove them
# from the dataset prophylactically.
grevys <- delete_vertices(grevys, "05s")
grevys_interactions_m <- grevys_interactions_m[, colnames(grevys_interactions_m) != "05s"]
grevys_interactions_m <- grevys_interactions_m[rownames(grevys_interactions_m) != "05s",]

# We don't know any attributes for the zebras, so we can just color and sort again by degree
grevy_sequence <- reorder(V(grevys)$name, igraph::degree(grevys, mode = "total"))
# And color by degree
colorscale_grevys <- turbo(1+max(igraph::degree(grevys, mode = "total")))
V(grevys)$color <- colorscale_grevys[1+igraph::degree(grevys, mode = "total")]  
color_sequence <- V(grevys)[levels(grevy_sequence)]$color

circos.par(start.degree = 90 # This is starting at 'up' and going clockwise
           , canvas.xlim = c(-1.5,1.5)) 
chordDiagram(grevys_interactions_m
             , order = levels(grevy_sequence)
             , grid.col = color_sequence #needs to be length nodelist
             #, col = "grey" #needs to be length of nodesxnodes (i.e. whole matrix)
             , symmetric = TRUE # For undirected graphs, so each link won't be plotted twice
             , annotationTrack = c("grid")
             , annotationTrackHeight = c(0.1)
             , scale = TRUE # this makes each sector the same size and shows relative strength of connections
)
# Customizing the labeling in this function is somewhat laborious; you 
# have to separately format the 'track' (label circle)
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  #circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
  #            facing = "bending.inside", niceFacing = TRUE, col = "white", cex = 0.5)
  circos.text(mean(xlim) # 'x-axis' is around the circle
              , ylim[1] + 1.2 # 'y-axis' is radial, i.e. away from circle
              , paste(si, "deg:", igraph::degree(grevys, mode = "total")[si])
              , sector.index = si
              , track.index = 1 # 'track' is the concentric ring in which labels are
              , facing = "clockwise" # orientation of text
              , niceFacing = TRUE
              , col = "black"
              , cex = 0.8 # text size
              , adj = c(0, 0.5)
  )
}
circos.clear() # required after changing circos.par()

## If we were also doing community detection, one could label and set a bit apart
## the fully or partially disconnected modules on this network. 

## Problems with chordDiagram still:
# I am not sure about the interaction of in- vs out-degrees and how the function decides
# sector size. In the scaled version this is irrelevant, but otherwise the calculated 
# degree doesn't agree with the sector size as plotted. 
warnings()


# Arc diagrams ---------------------------------------
# https://rdrr.io/github/gastonstat/arcdiagram/f/vignettes/introduction.Rmd
detach("package:circlize", unload=TRUE)

library(devtools)
devtools::install_github('gastonstat/arcdiagram')
library(arcdiagram)
weights <- edge_attr(UScities, "how_close")
arcplot(as_edgelist(UScities), lwd.arcs = weights, show.nodes=TRUE, las=3)
# See https://www.r-bloggers.com/2013/02/arc-diagrams-in-r-les-miserables/ for
# a fully formatted example. 

par(oma=c(1,0,0,0)
    , mar=c(4, 1, 0, 1)
    , xpd=TRUE) # bottom, left, top, right

arcplot(as_edgelist(grevys)
        , lwd.arcs = E(grevys)$weight
        , show.nodes=TRUE
        , las=3
        , cex.node = 1 + degree(grevys)/10
        , col.nodes = V(grevys)$color
        , col.labels = V(grevys)$color
        , line = 1
        , col.arcs = alpha("black", 0.2) #making sure arcs are transparent
        # is important for this kind of plot, so that multiple arcs on 
        # top of each other have different color 
)

warnings()

# Spatial arrangement of nodes ----------------------

## Different algorithm for automated rearrangement
# With SNA
gplot(dolphins_sna, mode="circle")
gplot(dolphins_sna, mode="random")
# With iGraph
plot(dolphins, layout = layout_with_fr)

# And using multidimensional scaling, like PCA, we can 
# arrange them so that distances in the plot correspond to 
# degree of similarity in the network (in terms of set of nodes connecting to).
par(oma=c(0,0,0,0)
    , mar=c(4, 4, 1, 1)
) # bottom, left, top, right
loc <- cmdscale(dolphin_num_matrix)
x <- loc[, 1]
y <- loc[, 2]
plot(x, y
     , xlab = "First dimension"
     , ylab = "Second dimension"
     , axes = TRUE
     , pch = 19
     , col = colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))]
     , cex = 1.5 # point size
)
text(x, y-0.03, rownames(loc), cex = 0.6)

# Or, we can use the positions derived from cmdscale to drive the network plot.
plot(dolphins
     , vertex.color = alpha(colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))], 0.4)
     #, vertex.label.cex = 0.7
     #, vertex.label.color = "black"
     #, vertex.label.family = "Calibri"
     , vertex.size = 1+degree(dolphins, mode = "all")/2
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = loc
     , vertex.label = "" # V(dolphins)$name
     , rescale = FALSE
     , axes = TRUE
     , xlim = c(-0.6, 0.6)
     , ylim = c(-0.6, 0.6)
     , pch = 19
)
text(x, y, rownames(loc), cex = 0.8)
# We can see this is not ideal because a bunch of nodes are too close to see much;
# this is what the automated network plot algorithms help with. But the information here
# is that the males are all mores similar to each other, connection-wise, than
# the females, who end up on the outside of this plot. 
# Larger nodes (with higher degree, more connections) are not necessarily similar to each
# other in who they connect to. 


# In a weighted network, this method means deciding locations of points based on 
# strength of ties in a weighted network.
# Before you run the next command, think about what you expect will happen - what
# will this plot look like?

loc <- cmdscale(UScitiesD)
x <- -loc[, 1]
y <- -loc[, 2] 
plot(x, y
     #, type = "n"
     , xlab = "", ylab = "", axes = TRUE
     , asp = 1
     , pch = 19
     , col="blue"
)
text(x, y-50, rownames(loc), cex = 0.6)
## Note asp = 1, to ensure an aspect ratio of 1:1 between axes.
## And the '-' in front of loc for x and y flips axes so North is at the top
## and West is left. 


# Plotting on actual maps -------------------
library('maps')
library('geosphere')

# I just want a more or less random example dataset, but if you are interested
# in available species occurrence data, check out this description of how to use
# the package and data: https://rspatialdata.github.io/species_occurrence.html
install.packages('spocc', dependencies = TRUE)
library(spocc)
occurrence_data_gbif <- occ(query = "Pheidole"
                       , from = "gbif"  
                       , limit = 1000 # just because right now we don't need to work with a huge or complete dataset
                       , has_coords = TRUE
                       )
occurrence_data_inat <- occ(query = "Pheidole"
                       , from = "inat"  
                       , limit = 1000 
                       , has_coords = TRUE
)
#head(occurrence_data_gbif$gbif$data$Pheidole, 3) 
# The reduced dataset is this:
Pheidole_US <- rbind(occ2df(obj = occurrence_data_gbif), occ2df(obj = occurrence_data_inat))
#head(Pheidole_US)
Pheidole_US$latitude <- as.numeric(Pheidole_US$latitude)
Pheidole_US$longitude <- as.numeric(Pheidole_US$longitude)
Pheidole_US <- subset(Pheidole_US, (latitude >= 20 & latitude <= 50))
Pheidole_US <- subset(Pheidole_US, (longitude >= -125 & longitude <= -65) ) 

# Plot a map of the united states:
map("state", col="grey20", fill=TRUE, bg="black", lwd=0.1)

# Add a point on the map for each occurrence:
points(x=Pheidole_US$longitude, y=Pheidole_US$latitude
       , pch=19
       , col= c("orange", "blue")[as.factor(Pheidole_US$prov)]
)
# The species occurrence data above are real, but now I'm going to make up some fake interaction
# data across the map just to demonstrate how we would plot them:
nodes1 <- unique(Pheidole_US$key)
nodes2 <- sample(unique(Pheidole_US$key))
fake_links <- data.frame(nodes1, nodes2)
fake_links$Weight <- rep(0, times = nrow(fake_links))
colnames(fake_links) <- c("Source", "Target", "Weight")
for(i in 1:nrow(fake_links)) {
  ifelse(runif(1)>0.6
         , fake_links[i,3] <- runif(1)
         , 0)
}
edge.col <- alpha(magma(100), 0.3)

for(i in 1:nrow(fake_links))  {
  node1 <- Pheidole_US[Pheidole_US$key == fake_links[i,]$Source,][1,]
  node2 <- Pheidole_US[Pheidole_US$key == fake_links[i,]$Target,][1,]

  arc <- gcIntermediate(c(node1$longitude, node1$latitude)
                        , c(node2$longitude, node2$latitude)
                        , n=1000
                        , addStartEnd=TRUE)
  
  edge.ind <- round(100*fake_links[i,]$Weight / max(fake_links$Weight))
  
  lines(arc
        , col=edge.col[edge.ind]
        , lwd=edge.ind/20)
}





