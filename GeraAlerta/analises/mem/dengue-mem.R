'''
Reads historical incidence data and extract epidemic thresholds per AP,
using MEM package https://cran.r-project.org/web/packages/mem/

:apsids: list with APs name.
:epithresholds: list with full epimem report for each ap, in the same order as apsids
:dfthresholds: dataframe with epidemic threhsold per AP, with pre, pos, mid, high and very high
               intensity thresholds
'''

library(mem)
library(plyr)

bindseason <- function(df1, df2, baseyear){
  # Function to bind season incidences, placing each season in a colum
  
  ti <- (baseyear+1)*100
  tf <- (baseyear+1)*100 + 27
  df3 <- cbind(df2, df1[(df1$SE > max(df1$SE[df1$SE<ti])-26 & df1$SE < tf),
                      c('SE', 'inc')])
  suff <- paste(as.character(baseyear),as.character(baseyear+1), sep='-')
  newse <- paste0('SE',suff)
  newinc <- paste0('inc',suff)
  df3 <- rename(df3, c('SE'=newse, 'inc'=newinc))
  
  return(df3)
}

# Read historical data
dfcomplete <- read.csv('alertaAPS_201539.csv')

# Store only necessary data, separating seasons by columns
dfsimple <- dfcomplete[dfcomplete$SE > max(dfcomplete$SE[dfcomplete$SE<201100])-26 &
                         dfcomplete$SE < 201127,
                       c('APS', 'SE', 'inc')]
dfsimple <- rename(dfsimple, c('SE'='SE2010-2011', 'inc'='inc2010-2011'))
seasons <- c('inc2010-2011')
for (i in 2011:2014){
  dfsimple <- bindseason(dfcomplete, dfsimple, i)
  seasons <- cbind(seasons, paste0('inc',i,'-',i+1))
}

# List of APSs
apsids <- unique(dfsimple$APS)

# Apply mem algorithm to obtain epidemic thresholds, using 0.60 confidence interval
epithresholds <- NULL
dfthresholds <- data.frame(apsids)
dfthresholds <- rename(dfthresholds, c('apsids'='aps'))
dfthresholds['pre'] <- NULL
dfthresholds['pos'] <- NULL
dfthresholds['mid'] <- NULL
dfthresholds['high'] <- NULL
dfthresholds['veryhigh'] <- NULL
for (aps in apsids){
  # Firstly, use all seasons
  epitmp <- epimem(i.data=subset(dfsimple[dfsimple$APS==as.character(aps),], select=seasons),
                   i.n.max=10, i.level=0.60, i.level.threshold=0.60)
  
  # Plotting structure:
  #plot(epitmp)
  #title(main=aps)
  
  # Discard seasons that are below threshold and rerun.
  # This is useful for properly defining activity levels during an epidemic
  discard <- NULL
  #for (ss in season){
  #  # Obtain seasons below threshold 
  #  # discard <- 
  #}
  episeasons <- seasons[! seasons %in% discard]
  epitmp <- epimem(i.data=subset(dfsimple[dfsimple$APS==aps,], select=episeasons),
                   i.n.max=10, i.level=0.60, i.level.threshold=0.60)
  
  # Store full report in epithresholds:
  epithresholds <- rbind(epithresholds, list(epitmp))
  
  # Store epidemic thresholds
  dfthresholds$pre[dfthresholds$aps==aps] <- epitmp$pre.post.intervals[1,3]
  dfthresholds$pos[dfthresholds$aps==aps] <- epitmp$pre.post.intervals[2,3]
  dfthresholds$mid[dfthresholds$aps==aps] <- epitmp$epi.intervals[1,4]
  dfthresholds$high[dfthresholds$aps==aps] <- epitmp$epi.intervals[2,4]
  dfthresholds$veryhigh[dfthresholds$aps==aps] <- epitmp$epi.intervals[3,4]
}