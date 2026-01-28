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

# Machine learning for network data ------------------------------
# If you want to learn more about these machine learning algorithms, I like the 
# tutorial provided by
# https://lgatto.github.io/IntroMachineLearningWithR/index.html

# Of course you may have multiple dimensions of information about nodes outside
# of the network itself, and these may be relevant to categorize nodes. 
# As far as the network data themselves are concerned though, before we simply 
# interpreted each row in the adjacency matrix as an attribute vector for each
# node, so that the 'weight of link to node A' was one dimension of the attributes
# of each node, and 'weight of link to node B' was another such dimension.

net <- dolphins
net_adjmatrix <- as.matrix(dolphins_df)
#net <- wildassN
#net_adjmatrix <- wildasses_interactions_m[,-5]

# 'Unsupervised learning' - i.e. categorization of multidimensional data ---------------
## Center and scale, or not?  ------------------
# Often, especially if different dimensions with different units are to be
# compared, it makes sense to center and scale values first. With network
# data this is a bit odd, since the 'weight of link to node A' across nodes
# is almost certainly not normally distributed. 
scaled_matrix <- scale(net_adjmatrix, center = TRUE, scale = TRUE)

## k-means clustering  ------------------
colorscale <- viridis(5)
node_coordinates <- layout_nicely(net)
par(mfrow = c(2,2), oma = c(0,0,0,0), mar = c(0,0,0,0))
output <- kmeans(net_adjmatrix, centers = 3, nstart = 10)
outputS <- kmeans(scaled_matrix, centers = 3, nstart = 10)
plot(net
     , vertex.color = colorscale[output$cluster]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , vertex.size = degree(net) + 3
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = node_coordinates
)
plot(net
     , vertex.color = colorscale[outputS$cluster]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , vertex.size = degree(net)+ 3
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = node_coordinates
)
output <- kmeans(net_adjmatrix, centers = 5, nstart = 10)
outputS <- kmeans(scaled_matrix, centers = 5, nstart = 10)
plot(net
     , vertex.color = colorscale[output$cluster]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , vertex.size = degree(net)+ 3
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = node_coordinates
)
plot(net
     , vertex.color = colorscale[outputS$cluster]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , vertex.size = degree(net)+ 3
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = node_coordinates
)
# In this output, the two graphs on the top have 3 clusters, the two on the
# bottom have 5; and the two graphs on the right are on scaled data.
# The most obvious difference is that when scaling, there is one cluster of 
# small to medium sized nodes distributed throughout. 

# How do we choose how many clusters to define? K-means clustering is based on
# finding the 'middle' of each cluster and minimizing the distance of each point
# to the middle of 'its' cluster. The sum of distances of all points is 
output$tot.withinss
# We can check what number of clusters gets us most of the way to a minimum:
tot_dist <- numeric(0)
for (i in 1:(vcount(net)/2)) {
  output <- kmeans(net_adjmatrix, centers = i, nstart = 10)
  tot_dist <- c(tot_dist, output$tot.withinss)
}
par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(4,5,1,1))
plot(tot_dist, xlab = "Number of clusters"
     , ylab = "1-'fit' [total sum of squares of distances\n of nodes from cluster centers]")
# If you did this with the dolphin data, there really ARE NO obvious clusters
# - the total distance of nodes from cluster centers just keeps decreasing with 
# more clusters, without an obvious 'elbow' indicating that we've captured most
# of the variation. 

## hierarchical clustering  ------------------
distance_object <- dist(net_adjmatrix)
hier_clust <- hclust(distance_object)
plot(hier_clust)
# The hier_clust object contains the information of the dendrogram.
# To cut the dendrogram into a predefined number of clusters, we can
# use 'cutree', which outputs a list of the nodes with their cluster IDs.
# This in turn can be used to define the colors in a plot, for example:
plot(net
     , vertex.color = colorscale[cutree(hier_clust, k = 3)]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , vertex.size = degree(net)+ 3
     , edge.arrow.size = 0
     , edge.width = 2
     , layout = node_coordinates
)


## PCA ------------------
# We already used a multidimensional scaling method like this:
loc <- cmdscale(net_adjmatrix)
x <- loc[, 1]
y <- loc[, 2]
plot(x, y
     , xlab = "First dimension"
     , ylab = "Second dimension"
     , axes = TRUE
     , pch = 19
     , col = colorscale[outputS$cluster]
    # , col = colorscale[cutree(hier_clust, k = 3)]
     , cex = 1.5 # point size
)
text(x, y-0.03, rownames(loc), cex = 0.6)
# (Note the colors here are either from the k means or hierarchical clustering,
# but you could use something completely different...)
# Now this is conceptually very similar to a PCA, but technically PCA is a 
# specific version of multidimensional scaling methods, and there are several
# different types - for more details check out for example
# https://stats.stackexchange.com/questions/14002/whats-the-difference-between-principal-component-analysis-and-multidimensional
# https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_07_Functions_PCA.pdf

# For a implementation of a more classic PCA that also gives more details, check out
net_pca <- prcomp(net_adjmatrix)
summary(net_pca) # shows what proportion of the variation is explained by each 
# principal component. 
plot(net_pca$x[, 1:2]
     , pch = 19
     , col = colorscale[outputS$cluster]
)
text(net_pca$x[,1], net_pca$x[,2]-0.03, rownames(loc), cex = 0.6)
# (Note again I colored here using the output from the k means clustering, 
# showing again some similarity but not an obvious correspondence.)

# All the algorithms above aim to classify a set of observations, for which
# you have multidimensional data, by finding their similarity. In some of these
# algorithms you define the number of groups you want to identify, in some
# a 'map' of similarity is created, in others the algorithm finds a number of 
# groups (categories).

# 'Supervised learning' - let R learn YOUR categories ----------
# These methods are not so much about finding structure in your data when you
# don't yet know if there is any; instead, they help you classify data points
# when you have a particular categorization in mind, but don't necessarily know
# what combination of the multidimensional information you have about each point
# is characteristic of each category. 

# To try out this method, we'll import a dataset that has a bit more information
# and structure than the dolphins'. 
# If you haven't used it before, get this package from github:
# https://github.com/schochastics/networkdata/blob/8e4994bc261654d040b8438a6568c6087db9dc60/README.md
# If necessary do 
# install.packages("remotes")
# first, then 
# remotes::install_github("schochastics/networkdata")
library(networkdata)
library("class")
# The following is a dataset on the interactions in the movie between characters
# in 'The Big Lebowski'. 
data("movie_103")
dude_net <- movie_103
# I'm going to manually add information on gender of characters as well as 
# whether these characters have an actual name or just a role description like
# 'waitress'. 
V(dude_net)$gender <- "M"
V(dude_net)$gender[c(5, 20,31,33,36)] <- "W"
V(dude_net)$named_char <- "Y"
V(dude_net)$named_char[c(2,3,7,12,14,19,21,23,25,28,31,33,35,36,37)] <- "N"
check_table <- data.frame(names = V(dude_net)$name, gender = V(dude_net)$gender, named = V(dude_net)$named_char)
dude_matrix <- as_adjacency_matrix(dude_net)

# Now we can try to classify nodes as being 'named characters' or not. 
# We'll first pick a few nodes as our training dataset:
tr <- sample(37, 15) # This picks 10 numbers from the set 1-37 (total number of
#                      characters)
nw <- seq(1:37)[-tr] # Rest of the numbers 1-37 that are not in the training set.

# First we'll try a k nearest neighbors algorithm.
# The actual classification is just one line:
knn_result <- knn(dude_matrix[tr, ], dude_matrix[nw, ], V(dude_net)$named_char[tr])
table(knn_result, V(dude_net)$named_char[nw])
# The diagonal shows the cases the algorithm classified correctly. Note that one 
# error is much more frequent than another: because some named characters are 
# misclassified as 'unnamed', but hardly any unnamed are misclassified as 'named'.
V(dude_net)$knn_result <- NA
V(dude_net)$knn_result[nw] <- knn_result
check_table$knn_result <- as.factor(V(dude_net)$knn_result)
levels(check_table$knn_result) <- c("N", "Y")

# There are many intricacies that we are not exploiting here. One is that there
# are parameters (arguments to the function knn() ) that can be changed; another
# is the choice and amount of training data. 

#install.packages("caret")
#install.packages("rpart")
#install.packages("rpart.plot")
#install.packages("caTools")
# The caret package has different tools for classification learning and 
# detailed outputs.
library(caret)
library(caTools)
# The 'trainControl' function allows us to predefine some settings for how
# the model training should proceed. 
# One particularly relevant setting here is the 'number' one, which controls 
# into how many training data vs test data sets the original data are split.
# The caret package function 'train' then trains a ML model repeatedly on 
# the different subsets of the data, using the parameters derived from this to 
# build the final overall model. 
myControl <- trainControl(
  method = "cv", ## cross validation
  number = 10,   ## 10-fold
  summaryFunction = twoClassSummary, ## NEW
  classProbs = TRUE, # IMPORTANT
  verboseIter = FALSE
)
# Most of the ML functions need the data as a data frame, with one column
# defining the correct classification.
dude_df <- as.data.frame(as.matrix(dude_matrix))
dude_df$Class <- as.factor(V(dude_net)$named_char)
# This is the actual model training. Here, for a test case, we are simply using
# a logistic regression, the parameters of which are 'trained' by the model,
# effectively using a bootstrapping method.
model <- train(Class ~ ., dude_df,
               method = "glm", ## to use glm's logistic regression
               trControl = myControl)
# Now that we've trained the model, we can see what classification the model
# 'predicts' for all the nodes. As a first step, we get a probability for each
# node
predicted_prob <- predict(model, dude_df, type = "prob")
# and in the second step we define a decision rule, like 'classify as "N" if the
# probability of N is >50%'. 
model_classification <- ifelse(predicted_prob[,1] > 0.5, "N", "Y")
# Now we want to see how well our model performed, which we can in this case, 
# since we know the 'true' classification of all the nodes in our dataset.
# The 'confusionmatrix' tabulates the true and false positives and negatives:
# (remember the way we set it up here is that "N" is the 'positive' detection)
# the prediction x reference table tells you how many "N" from the reference 
# (real classification) were detected by the model (prediction). So values in the
# diagonal of this table are correct classifications, off the diagonal are incorrect.
confusionMatrix(factor(model_classification), dude_df$Class)
# There are actually a lot of other outputs describing the model performance
# as well: accuracy, sensitivity (proportion of correct detections of all true positives),
# specificity (proportion of all positive detections that are actually true
# positives), etc.
# Another important tool for evaluating really any decision algorithm is the
# Receiver Operating Characteristic or ROC curve:
caTools::colAUC(predicted_prob, dude_df[["Class"]], plotROC = TRUE)
# In the ROC curve, instead of predefining the cutoff probability as we did in
# 'model_classification' above, the algorithm plots the proportion of correct
# positives (out of true positives) against the proportion of false positives
# (out of positive detections). If we refused to detect anything as positive at all,
# we'd be in the bottom left corner; we we classify everything as positive, we're
# in the top right corner. So no matter how good your model is at actually recognizing
# anything, the ROC curve should start from the bottom left and end at the top right.
# Better models skew to the top left, i.e. quickly increase sensitivity without 
# increasing false alarms. 
# Based on this particular ROC curve, we might decide that the best tradeoff 
# between sensitivity and specificity is at a different probability cutoff.
model_classification <- ifelse(predicted_prob[,1] > 0.3, "N", "Y")
confusionMatrix(factor(model_classification), dude_df$Class)

# In the above we were actually just using a logistic glm to do the classification.
# Below are a couple of examples of more sophisticated model fitting algorithms. 

## Random forest ------------------

library("rpart") ## recursive partitioning
library("rpart.plot")
# The 'trees' in the random 'forest' are decision trees. The algorithm generates
# a bunch of possible decision trees, and then evaluates their performance. 
# This is a versatile, non-linear method; but what I like particularly is that
# the output is easy to understand as a decision tree. 
model <- rpart(Class ~ ., data = dude_df, method = "class")
rpart.plot(model) # What this shows us is that the model essentially just uses
                  # 'Is the character connected to Walter?' as the test rule to 
                  # decide whether to estimate that hte character is a 'named'
                  # character.
plot(dude_net) # Just looking at the network to check 'Walter's' role shows that
# he is a sort of secondary main character, connected to most of the main characters
# (who have names), but not all the unnamed characters. 

predicted_prob <- predict(model, dude_df, type = "class")
table(predicted_prob, dude_df$Class)
# The classification based on this rule is not perfect (in fact considerably 
# worse than the glm-based one above), but given its simplicity, still quite good.

# Using the function from caret to implement a random forest model allows us
# to fine-tune a lot more parameters about how the search for the best decision
# tree proceeds, how many variables are taken into account, etc.; unfortunately
# the model produced can't be plotted as easily as an actual decision tree any more. 
model <- train(Class ~ ., data = dude_df, method = "ranger")
print(model)
# This initial version isn't very satisfactory, so we systematically try out a 
# set of values for tuning parameters:
myGrid <- expand.grid(mtry = c(1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35), # mtry has to be < than number of rows in data
                      splitrule = c("gini", "extratrees"),
                      min.node.size = 1) ## Minimal node size; default 1 for classification
model <- train(Class ~ .,
               data = dude_df,
               method = "ranger",
               tuneGrid = myGrid,
               trControl = trainControl(method = "cv",
                                        number = 5,
                                        verboseIter = FALSE))
print(model)
plot(model)
# Our accuracy here never got better than just using 'has a link with Walter' :-)

# -------------------------

