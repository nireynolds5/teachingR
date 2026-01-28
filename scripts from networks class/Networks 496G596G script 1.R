# Header-------------------------
## R script for ECOL496G 596G 
## Complex Systems and Networks
## UA class Spring 2025
## written by Anna Dornhaus

## This copy belongs to:

# R -----------------------------

## This is an R script. 
## Any line that starts with "##" is a comment, and directed to you, the user;
## it will be ignored by R.
## To learn how to use R, do each command by itself (rather than the whole script or sections
## at once). See what the outcome is. If you just get a new prompt, no error message etc.,
## it probably worked.

## If you are new to R, make sure you ask for help & tutorials to cover the basics. 

## Before you start, make sure you know where your script is saved on your computer
## (and any datasets you downloaded and are planning to use).

## Find out what working directory your R program is using:
getwd()
## and then set it to the directory/folder that contains your R files.
## E.g.
#setwd("C:/Users/Anna Dornhaus/Dropbox/Learning R")


# Preparation for the class -----------------------------

## There are a lot of packages for network analysis. In particular, network analysis needs 
## data structures that can contain both information about nodes and about links. 
## If you do it yourself, the simplest way is typically to think of the information
## about links as a table (or a matrix); but then if you want node properties, these need to 
## be in a separate vector. 
## Ready-made packages usually define their own data formats that contain both of these 
## types of information. 
## Most commonly used is the 'network' type or class (in package sna - for Social Network Analysis) and
## the 'igraph' type (in package iGraph). The package intergraph allows you to convert
## one into the other. 

## You only have to 'install' packages once on a given computer/RStudio instance - 
## afterwards, i.e. in new scripts or when running this script again, you only 
## need the 'library' command to activate them. 
#install.packages("sna")
#install.packages("intergraph")
#install.packages("igraph")

## I'm also going to use these packages for making color gradients:
#install.packages("viridis")
#install.packages("scales")

# Activate the Network Analysis packages ------------------

# Note that these packages use some of the same commands, and what the command does
# is determined by whichever package was last activated (with 'library()').
library(sna)
library(intergraph)
library(igraph)

# And for color/graphing:
library(viridis)
library(scales)

# Importing data -----------------------

## It is always good practice to import data directly from the original source, then 
## do all converting and cleaning you need in R, where it is documented. 

## It is also good practice to use informative, 'self commenting' variable names. 
## For example, if you are going to use a matrix format and an igraph format of the same 
## dataset, make sure they are stored in variables that tell you from the name which
## they are, to avoid difficult debugging later. 

## This is how you can load .csv files (stands for 'comma separated values')
## into R. The functions below generate two R variables, myMatrix and myMatrix2,
## which both contain tables from the .csv files. The second one includes both
## row and column names.
#dataset1_df <- read.table("data.csv", sep = ",")

## I named this '_df' to indicate it is of the 'dataframe' variable type. You can check:
#class(dataset1_df)

## Or, load like this:
#dataset2_df <- read.table("data_w_headers.csv", header =T, row.names=1, sep = ",")

## If you do some manipulations and then want to save your table again, use e.g.
#write.table(dataset2_df, "./mod_data.csv", col.names=NA, sep=",")


# Where to get data I: built-in ---------------------------

## I'm hoping you will use your own research data, as well as scouring the internet 
## for either fun or relevant datasets by others. R also has some built-in ones that
## are automatically available for you. 
## For example
UScitiesD
eurodist
data(emon)
data(flo)
data(coleman)

## Confusingly, these are all different formats though. To make them all into igraphs:

florentines <- graph_from_adjacency_matrix(flo)
emon1 <- emon[[1]]
organizations <- asIgraph(emon1)
V(organizations)$name <- V(organizations)$vertex.names
# UScities <- graph.adjacency(as.matrix(UScitiesD))
# The above command does make the UScitiesD dataset into an igraph object. Try it out. 
# Maybe plot it. Do you notice what is odd (and really unwieldy) here?
# A better, though a bit convoluted strategy is to do the following:
dist_matrix <- as.matrix(UScitiesD)
dist_matrix[dist_matrix>0] <- 1
UScities <- graph_from_adjacency_matrix(dist_matrix)
dist_list <- as.list(as.matrix(UScitiesD))
omit_zeros_list <- dist_list[dist_list != 0]
dist_to_weight <- function(x) {return(2000/x)}
E(UScities)$how_close <- lapply(omit_zeros_list, dist_to_weight)
# Who can tell me what this does?
# Remember to do the same thing to the eurodist variable if you want to use it
# with igraph. 

# If you are new to R, or if you don't yet want to be bothered with the details of 
# implementation, feel free to simply run the code above and forget all about it. 
# Just use the newly defined igraph objects from here on out:
florentines
organizations
UScities

# Where to get data II: download ---------------------------

# On D2L, I'm also giving you a dataset about interactions between dolphins. You
# will have to download it and put it into your working directory before running
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Lusseau dolphins/"
dolphins_df <- read.table(paste(path, "dolphinnetwork.csv", sep = ""), sep = ",")
dolphins_nodeinfo <- read.table(paste(path, "dolphinvertexinfo.csv", sep = ""), header = TRUE, sep = ",")
dolphins <- graph_from_adjacency_matrix(as.matrix(dolphins_df))
V(dolphins)$sex <- dolphins_nodeinfo$dolphin.sex
V(dolphins)$disappeared <- dolphins_nodeinfo$dolphin.disappeared
# You can now just use the igraph 
dolphins

# On Google Drive (linked on D2L) I'm providing a number of other animal-related
# datasets that are downloadable from the internet. 
# Again, to run this you need to download the file first then enter the correct
# path for your computer:
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Memmott pollinator network/"
pollination_df <- read.table(paste(path, "memmott_1999.csv", sep = ""), header = TRUE, sep = ",")
pollination_df <- pollination_df[-1]
pollination_nodeinfo1 <- read.table(paste(path, "plants.csv", sep = ""), sep = ",")
pollination_nodeinfo2 <- read.table(paste(path, "pollinators.csv", sep = ""), header = TRUE, sep = ",")
colnames(pollination_df) <- pollination_nodeinfo1$V2
rownames(pollination_df) <- pollination_nodeinfo2$species
pollination <- graph_from_biadjacency_matrix(as.matrix(pollination_df), directed = FALSE, weighted = TRUE)


# iGraph for network description and plotting ----------------

# iGraph has a fairly nice and comprehensive tutorial here:
# https://r.igraph.org/articles/igraph.html

# For any of the igraph objects above you already have, you can check a list of
# their nodes (or 'vertices') with
V(florentines)
V(dolphins)
# and a list of the links (or 'edges') with
E(organizations)
E(UScities)

# In iGraph, you plot just like this:
plot(dolphins)
# But just like the base R 'plot', the iGraph one has many options - check the 
# tutorial. For example,
plot(UScities, edge.width = E(UScities)$how_close)

# If you downloaded the pollination network, you need a bipartite graph, like this:
plot(pollination, layout = layout.bipartite)
# Don't worry, we will talk soon about how to make this look better. 

# SNA for network description and plotting ----------------

# To convert any of the above into 'network' objects, use
dolphins_sna <- asNetwork(dolphins)

# Note that this conversion function comes from the intergraph package, which 
# has a tutorial here:
# https://cran.r-project.org/web/packages/intergraph/vignettes/howto.html

# The plotting function in SNA is gplot (not to be confused with ggplot...)
gplot(dolphins_sna)

## Check the file "M Clarkson gplot network graphing manual 1.pdf" (provided on D2L)
## for details on how to modify the appearance of the graph resulting from gplot.
## Remember that you can add any number of extra parameters, so try the function
## below with just one new element at a time, to see what each part does.
## For example
gplot(dolphins_sna, diag=TRUE, vertex.col=1, edge.lwd=1, loop.cex=2)
## adds, in this order, self-referential links (loops), nodes are black, edge thickness is
## scaled to interaction strength * 1, loops are 2x normal size.

# Visualize the matrix instead -----------------------

# I am using the color scales from the package 'viridis' since I find them attractive,
# they are suitable for color-blind people, AND they maximize color distinguishability.
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

# Since 'heatmap' invokes a temperature-related intensity, I am using 'magma'. 
heatmap(as.matrix(UScitiesD)
        , col = magma(2000)
        , margins = c(8,8)
        , scale = "none"
)
# Note that this is just like writing the numbers in the interaction matrix, but using darker
# colors for smaller numbers. 

# Compare the above to this.
heatmap(as.matrix(UScitiesD)
        , Rowv = NA
        , Colv = NA
        , col = magma(2000)
        , scale = "none"
        , margins = c(8,6)
)

# The above does not show a legend. This is the basic version, but note that if
# the values in the matrix are not scaled, the 'middle' value is not the average 
# or median of the matrix values, just the middle color. 
legend("bottomright"
       , legend = c("min", "middle", "max")
       , fill = magma(3)
)

# A fancier version using the package IMIFA
par(oma=c(0,0,1,3))
heatmap(as.matrix(UScitiesD)
        , col = magma(2000)
        , margins = c(8,8)
        , scale = "none"
)
library(IMIFA)
heat_legend(as.matrix(UScitiesD), magma(2000), side = 2)

# As always in R, there are many packages to do the same thing in slightly different ways or styles.
# Feel free to experiment. 
# E.g. for heatmaps there is also the package 'ComplexHeatmaps' https://github.com/jokergoo/ComplexHeatmap



