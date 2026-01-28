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

## Colorpalette
twogroupcolors <- c("#5a9a8f", "#7b6ea8")
threegroupcolors <- viridis(3)
severalcategories_colors <- magma(6)

## Similar to colors, it often makes sense to define some other things
## universally for all your plots, e.g. margins, where the sample sizes
## are plotted, y-axis range. This depends on your figures though.
y_max <- 6
y_min <- 0
N_y_offset <- y_max * 0.95 # This puts sample size numbers 5% below max, for example



# IMPORT ----------------------------------
MyData <- read.table("bb col data.csv"
                     , header=T
                     , row.names=1
                     , sep = ","
                     , dec = "."
                     )


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
par(oma = c(0,0,0,0), mar = c(3,5,1,1), mgp=c(4, 1, 0), las=1) # bottom, left, top, right
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
  , y=N_y_offset
  , cex = 1
  , col = threegroupcolors
  , paste(Nice_Plot$n,sep="")  # again, the sample size 'n' is directly extracted from the plot
)

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


