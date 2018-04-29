############################################################
# Useful Libraries
############################################################
library(ggplot2) # For data visualizations
library(parallel) # For parallel computing
library(sqldf) # For manipulating dfs using SQL
library(lubridate)
library(pomp)
require(foreach)
require(doMC)
library(dplyr)
library(tseries)
library(ggfortify) # For ggplot ACF
############################################################
# Useful Functions
############################################################


registerDoMC(cores=detectCores())        
## number of cores 
## usually the number of cores on your machine, or slightly smaller

set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)

##################################################################
# SQL Related
##################################################################
sql = function(str){ 
  #Function for fast querying with SQL
  require(sqldf)
  sqldf()
  ret<-sqldf(str,drv="SQLite")
  sqldf()
  return(ret)
}


##################################################################
# GG-Plot Related
##################################################################
formatGG = function(ggobject, 
                    axis_numbers_size = 10,
                    axis_titles_size = 11,
                    title_size = 17,
                    legend_pos = c(0.9,0.9)){
  require(ggplot2)
  #Function to preformat ggplot objects
  gg = ggobject + 
    theme_minimal()+
    theme(legend.position=legend_pos, legend.justification=c(1,1),
                legend.box = "horizontal",
                legend.background=element_rect(color="lightgrey"))+
    theme(legend.title=element_blank()) +
    theme(plot.title=element_text(face="bold",hjust=-.08,vjust=2,colour="#3C3C3C",size=title_size))+
    theme(axis.title.y=element_text(size=axis_titles_size,colour="#535353",face="bold",vjust=1.5)) +
    theme(axis.title.x=element_text(size=axis_titles_size,colour="#535353",face="bold",vjust=-.5))+
    theme(axis.text.x=element_text(size=axis_numbers_size))+
    theme(axis.text.y=element_text(size=axis_numbers_size))
  return(gg)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
