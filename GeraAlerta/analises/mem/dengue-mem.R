"
Reads historical incidence data and extract epidemic thresholds per AP,
using MEM package https://cran.r-project.org/web/packages/mem/

:apsids: list with APs name.
:epithresholds: list with full epimem report for each ap, in the same order as apsids
:dfthresholds: dataframe with epidemic threhsold per AP, with pre, pos, mid, high and very high
               intensity thresholds
"

library(mem)
library(plyr)
library(logging)

# Initialize logger
basicConfig()
addHandler(writeToFile, file="./dengue-mem.log", level='INFO')
bindseason <- function(df1=data.frame(), df2=data.frame(), baseyear=integer()){
  "
  Function to bind season incidences from df1 onto df2, placing each season in a new colum
  Returns:
  :df3: data.frame df2 with 2 new colums from df1$SE and df1$inc based on df1$SE within season
        of base year, which is from winter of baseyear to winter of baseyear+1 (south hemisphere).
  "
  
  if (missing(df1) | missing(df2) | missing(baseyear)){
    logerror('missing argument on function call', logger='dengue-mem.bindseason')
    return(NULL)
  }

    ti <- (baseyear+1)*100
  tf <- (baseyear+1)*100 + 41
  df3 <- cbind(df2, df1[(df1$SE > max(df1$SE[df1$SE<ti])-12 & df1$SE < tf),
                      c('SE', 'inc')])
  suff <- paste(as.character(baseyear),as.character(baseyear+1), sep='-')
  newse <- paste0('SE',suff)
  newinc <- paste0('inc',suff)
  df3 <- rename(df3, c('SE'=newse, 'inc'=newinc))
  
  loginfo('Function executed and exited with status 0', logger='dengue-mem.bindseason')
  return(df3)
}

applymem <- function(df.data, l.seasons){
  "
  Function to apply epimem algorithm on df.data and generate full reports as well as a data frame
  with summary of relevant thresholds (pre, pos, mid, high, veryhigh)

  Input:
  :df.data: data frame with APS in a column named APS, and incidence seasons in each column
            each row gives the incidence in each season, for each APS, for each week.
  :l.seasons: vector with incidence columns to be used

  Returns:
  :epithresholds: list with full epimem report for each APS, keyed by AP's name.
  :dfthresholds: data frame with thresholds for each AP.
  "

    if (missing(df.data)){
    logerror('missing argument on function call', logger='dengue-mem.applymem')
    return(NULL)
  }
  
  # List of APSs
  apsids <- unique(df.data$APS)
  
  # Apply mem algorithm to obtain epidemic thresholds, using 0.60 confidence interval
  epithresholds <- list()
  dfthresholds <- data.frame(apsids)
  dfthresholds <- rename(dfthresholds, c('apsids'='aps'))
  dfthresholds['pre'] <- NULL
  dfthresholds['pos'] <- NULL
  dfthresholds['mid'] <- NULL
  dfthresholds['high'] <- NULL
  dfthresholds['veryhigh'] <- NULL
  for (aps in apsids){
    # Firstly, use all seasons
    epitmp <- epimem(i.data=subset(df.data[df.data$APS==as.character(aps),], select=l.seasons),
                     i.n.max=10, i.level=0.60, i.level.threshold=0.60)
    
    # Discard seasons that are below threshold and rerun.
    # This is useful for properly defining activity levels during an epidemic
    discard <- NULL
    #for (ss in season){
    #  # Obtain seasons below threshold 
    #  # discard <- 
    #}
    episeasons <- l.seasons[! l.seasons %in% discard]
    epitmp <- epimem(i.data=subset(df.data[df.data$APS==aps,], select=episeasons),
                     i.n.max=20, i.level=0.60, i.level.threshold=0.60)
    
    # Store full report in epithresholds:
    epithresholds[[aps]] <- epitmp
    
    # Store epidemic thresholds
    dfthresholds$pre[dfthresholds$aps==aps] <- epitmp$pre.post.intervals[1,3]
    dfthresholds$pos[dfthresholds$aps==aps] <- epitmp$pre.post.intervals[2,3]
    dfthresholds$mid[dfthresholds$aps==aps] <- epitmp$epi.intervals[1,4]
    dfthresholds$high[dfthresholds$aps==aps] <- epitmp$epi.intervals[2,4]
    dfthresholds$veryhigh[dfthresholds$aps==aps] <- epitmp$epi.intervals[3,4]
    
  }
  loginfo('Function executed and exited with status 0', logger='dengue-mem.applymem')
  return(list("epimemthresholds"=epithresholds, "dfthresholds"=dfthresholds))
}

# "
# Example of usage, based on ./alertaAPS_201539.csv input file.
# Applies method and generate full report and plots for each AP.
# "
# # Read historical data
# dfcomplete <- read.csv('./alertaAPS_201539.csv')
# 
# # Store only necessary data, separating seasons by columns
# dfsimple <- dfcomplete[dfcomplete$SE > max(dfcomplete$SE[dfcomplete$SE<201100])-12 &
#                          dfcomplete$SE < 201141,
#                        c('APS', 'SE', 'inc')]
# dfsimple <- rename(dfsimple, c('SE'='SE2010-2011', 'inc'='inc2010-2011'))
# seasons <- c('inc2010-2011')
# for (i in 2011:2014){
#   if (max(dfcomplete$SE) >= (i+1)*100 + 41){
#     dfsimple <- bindseason(dfcomplete, dfsimple, i)
#     seasons <- cbind(seasons, paste0('inc',i,'-',i+1))
#   }
# }
# 
# thresholds <- applymem(dfsimple)
# # Plotting structure:
# 
# for (aps in unique(dfsimple$APS)){
#   sprintf("APS: %s\n", aps)
#   print(thresholds$epimemthresholds[[aps]])
#   plot(thresholds$epimemthresholds[[aps]])
#   title(main=aps)
# }