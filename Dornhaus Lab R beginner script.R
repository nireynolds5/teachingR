# Anna Dornhaus & Lab
# R script for data analysis and figures for 
# just learning to do it in R

# Paper reference: NONE

### Libraries ------------------
#library(lme4) #needed for GLMMs
#library(lmerTest) #needed for obtaining p-values in lmm
#library(emmeans) #post-hoc comparisons
#library(sjPlot) # for linear models output tables
#library(tidyverse)

library(scales)  # for number_format & color scales
library(viridis)

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
# Or you can load a .csv file from Google Drive:
path <- "C:/Users/dornh/Dropbox/TEACHING/ECOL496G 596G Spring2025 Networks/network datasets/Memmott pollinator network/"
MyDatafromGoogleDrive <- read.table(paste(path, "memmott_1999.csv", sep = ""), header = TRUE, sep = ",")

# Or directly from a Google Sheet:
MyDatafromGoogleSheets <- read.csv("https://docs.google.com/spreadsheets/d/1K2D2rH770iDfZzlwakHK89bwvoPqWDLaHbVNsanJFX8/gviz/tq?tqx=out:csv")
# Or directly from github:
MyDatafromGithub <- read.csv("https://raw.githubusercontent.com/shannonmcwaters/Directed-exploration/refs/heads/main/Maze%20Data%20Raw")


# GRAPHING THINGS ------------------------------

# Bare boxplot
boxplot(avgsize ~ treatment
        , data = MyData
        , xlab = "X Concept [unit measured]"
        , ylab = "Y Concept [unit measured]"
        , col = threegroupcolors
        , range = 0
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






