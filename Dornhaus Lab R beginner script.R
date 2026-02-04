# Anna Dornhaus & Lab
# R script for data analysis and figures for 
# just learning to do it in R

# Paper reference: NONE

### Libraries ------------------
#library(lme4) #needed for GLMMs
#library(lmerTest) #needed for obtaining p-values in lmm
#library(emmeans) #post-hoc comparisons

# General data handling
library(tidyverse) # for bind.rows
# File access
library(readxl)
library(googlesheets4) # for working with Google Sheets
library(googledrive) # for drive_download
# Colors
library(scales)  # for number_format & color scales
library(viridis)
# Making output tables for statistical models
library(sjPlot) # for linear models output tables


# Graphics setup -----------------------------
no_categories <- 6
## Colorpalette
twogroupcolors <- c("#5a9a8f", "#7b6ea8")
threegroupcolors <- viridis(3)
severalcategories_colors <- magma(no_categories)

## Similar to colors, it often makes sense to define some other things
## universally for all your plots, e.g. margins, where the sample sizes
## are plotted, y-axis range. This depends on your figures though.
y_max <- 6
y_min <- 0
N_y_offset <- 0.95 # This puts sample size numbers 5% below max, for example



# IMPORT ----------------------------------

# A common problem when reading any files is that you must make sure
# that R uses the folder you want as 'working directory'. You can
# pick the working directory from the 'Session' menu above; or 
# you can write the path here in the code, e.g.:
setwd("C:/Users/dornh/Dropbox/Github/teachingR")
# But note that this is a little problematic as this detailed path
# would never work on someone else's computer. Instead,
# you can also use
setwd("../teachingR")
# or similar ... the ".." means 'go up one directory level' and then
# it will go into the folder to the right of the "/". 
MyData <- read.table("bb col data.csv"
                     , header=T
                     , row.names=1
                     , sep = ","
                     , dec = "."
                     )

# Or directly from a Google Sheet:
MyDatafromGoogleSheets <- read.csv("https://docs.google.com/spreadsheets/d/1K2D2rH770iDfZzlwakHK89bwvoPqWDLaHbVNsanJFX8/gviz/tq?tqx=out:csv")
# Or an existing .csv file on Google Drive:
MyDatafromGoogleSheets <- read.csv("https://drive.google.com/uc?export=download&id=1woJYnpsfMBPCC3kPSErcydgDCduEKjzy")
# Note that in both cases you are not just pasting the entire link from your browser, you are modifying it to include 
# the file id and some other stuff. Pattern match here!

# An .xlsx file from Google Drive is more complicated:
# You have to first download the file into a local file, then
# import the local file into R. 
googlepath <- "https://docs.google.com/spreadsheets/d/1sF5WnJs0uxeLKApn7_JzGJMozu6wgSPJ/edit"
tempxlsfile <- tempfile(fileext = ".xlsx")
drive_download(as_id(googlepath), path = tempxlsfile, overwrite = TRUE)
MyDatafromGoogleDrive  <- read_excel(tempxlsfile)
unlink(tempxlsfile)

# Or directly from github:
MyDatafromGithub <- read.csv("https://raw.githubusercontent.com/shannonmcwaters/Directed-exploration/refs/heads/main/Maze%20Data%20Raw")

# Import a whole list of files
# Note that this will only work if all the files have the same columns (and otherwise
# it is anyway doubtful that this would make sense)
files <- (Sys.glob("./example_data/similardatasheets/*.csv"))
# Initiate a blank data frame
EnormousData <- data.frame()
# Read content of all files into a list
listOfDataframes <- lapply(files, 
                             function(x) {
                               read.table(x, 
                                          header = T,
                                          sep = ",",
                                          skip = 6
                               )
                             }
)
# Add all the rows from all the files together
EnormousData <- do.call("bind_rows", listOfDataframes)

# Import a whole list of xls files
# Headache! Saving as csv is better...
googledrivefolderfiles <- drive_ls(as_id("https://drive.google.com/drive/folders/1aZWvhJjgTH9QTrD9hV1LiPiOKFRxmH7x"))
files2 <- googledrivefolderfiles$id
names2 <- googledrivefolderfiles$name
files_local <- vector(mode = "character", length = length(files2))
dir.create(file.path("./", "temp"), showWarnings = FALSE)
for (i in 1:length(files2)) {
  # Put name in files_local
  files_local[i] <- paste("./temp/", names2[i])
  # Download file
  drive_download(files2[i], path = files_local[i], overwrite = TRUE)
}
EnormousData <- data.frame()
# Read content of all files into a list
listOfDataframes <- lapply(files_local, 
                           function(x) {
                             read_excel(x)
                           }
)
# Add all the rows from all the files together
EnormousData <- do.call("bind_rows", listOfDataframes)


# GRAPHING THINGS ------------------------------

### BOXPLOTS ### ---------------
# Generally you use a boxplot when plotting a continuous y
# against a categorical x axis (e.g. an outcome against 2 treatments).

# Bare boxplot
boxplot(avgsize ~ treatment
        , data = MyData
        , xlab = "X Concept [unit measured]"
        , ylab = "Y Concept [unit measured]"
        , col = threegroupcolors
        , range = 0 # to get whiskers to extend to range (no outlier points)
)

# Fancy boxplot - but this should be your default
# What makes it fancy:
# Adjust margins
# Add sample sizes
# Add data points

# Margins:
par(oma = c(2,2,2,2), mar = c(4,4,1,1), mgp=c(3, 1, 0), las=1) 
# bottom, left, top, right
par(mfrow=c(1,1))

# How to make a great boxplot -

# Saving the plot into a variable allows us to access plot parameters afterwards.
Nice_Plot <- boxplot(avgsize ~ treatment
                     , data = MyData
                     , xlab = "X Concept [unit measured]"
                     , ylab = "Y Concept [unit measured]"
                     , range = 0
# Always make axis descriptions as clear and comprehensive as possible
                     , names = c("Large bees", "Middling bees", "Small bees")
                     , col = alpha(threegroupcolors, 0.5) # use same colors as elsewhere, 
                     # but slightly transparent so we can see data points
                     , ylim = c(y_min, y_max) # always think about the scale - starting from zero is typically better
)

# Putting sample sizes above bars
nbGroup <- nlevels(as.factor(Nice_Plot$names)) # this is just a way to extract
# category names from the plot - you could get this directly from data
text(x=c(1:nbGroup) 
  , y=N_y_offset*y_max
  , cex = 1
  , col = threegroupcolors
  , paste("N=", Nice_Plot$n, sep="")  # again, the sample size 'n' is directly extracted from the plot
)

mtext("additional margin label", side=2, line=2, las=0)
mtext("additional margin label on outside", side=1, line=5, las=0, xpd = TRUE)

# Represent the raw data as well, especially for mid- to low sample sizes.
stripchart(avgsize ~ treatment
           , data = MyData
           , add = TRUE # this plots this graph on top of the existing one
           , pch = 19
           , col = threegroupcolors
           , method = "jitter"
           , jitter = 0.2
           , vertical = TRUE
)



### SCATTERPLOTS ### ---------------

### NOT DONE !!!!!!!!!!!!!!!! -------------

# working on this 

graph_data <- MyData
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
boxplot(table(graph_data$avgsize)
        , data = graph_data
        , xaxt = 'n'
        , yaxt = 'n'
        , frame = FALSE
        , range = 0
#        , col = color_distributions
)
# Panel 2: The original graph
par(mar = c(4,4,1,2)) # bottom, left, top, right
plot(avgsize ~ stdevsize # we could have stuck with the 
     # subsets above, but this is a little cleaner - note that we still assign
     # colors by the variable $sex in the node list
     , pch = 19 # set point shape
#     , col = colors_dolphinsex[sapply(V(dolphins)$sex, function(x) switch(x, "M"=1, "F"=2, "UNKNOWN"=3))]
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
boxplot(table(stdevsize)
        , data = graph_data
        , xaxt = 'n'
        #, yaxt = 'n'
        , frame = FALSE
        , range = 0
        , col = color_distributions
        , horizontal=TRUE
)

### SCATTERPLOT CODE IS NOT DONE!!

