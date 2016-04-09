#packages
library(reshape)
library(dplyr)
library(ggplot2)

####COLORS####
night <- rgb(26,31,30, max=255)
beach <- rgb(227,223,186, max=255)
red60 <- rgb(1,.4,.4)
tangyblue <- rgb(108,189,181, max=255)
purp <- rgb(181.5,145.5, 141.5, max=255)


green1 <- rgb(178,179,159, max = 255)
lightbrown <- rgb(205,140,82, max = 255)

green2 <- rgb(200,201,181, max = 255)
green3 <- rgb(222,223,197, max = 255)
orange <- rgb(61,66,60, max = 255)
aqua <- rgb(240,236,201, max = 255)

####TIME SERIES####

#tscs.lag <- function(data,unit,var,n=1){
#  library(dplyr,quietly = T)
#  lag.name <- paste('l',var[1],sep='')
#  lagged.data <- data %>% group_by(unit[1]) %>% mutate(lag.name=lag(var,n=n))
#  return(lagged.data)
#}

#lagging time series variables by 1 time unit
#data %>% group_by(unit) %>% mutate(lvar=lag(var))
