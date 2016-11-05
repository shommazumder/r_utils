####GENERAL FUNCTIONS FOR R####
##AUTHOR: SHOM MAZUMDER
##

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

olive <- rgb(164,199,6,max = 255)
orangejuice <- rgb(252,135,16,max = 255)
flymetothemaroon <- rgb(159,25,22,max = 255)

####DATA PREP####

stata_codebook <- function(x) {
  cb <- data.frame(attr(x, "var.labels"))
  rownames(cb) <- names(x)
  cb
}

####AUXILLARY####
got_Internet <- function(){
  if (is.character(RCurl::getURL("www.google.com"))) {
    out <- TRUE
  } else {
    out <- FALSE
  }
}

####REGRESSION####

get_rse <- function(models.list,hc='HC2'){
  ##Purpose: this function takes in a list of regression models and outputs a list of robust SEs
  se.list <- lapply(models.list,FUN = function(x){return(sqrt(diag(vcovHC(x,type = hc))))})
  return(se.list)
}

ch.row <- function(name, yesno) {
  ##Purpose: adds checkmarks to table
  c(name, ifelse(yesno, "$\\checkmark$", ""))
}

####GIS####


####PLOTTING####
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

##Coefficient Plots##

myCoefPlot <- function(models = models.list,names = model.names,coef.name = coef.name,se = NULL){
  #takes in regression models, the name of the coefficient you want to plot, and an optional list of standard errors
  #and returns a ggplot coefficient plot with 95% confidence intervals
  
  #initialize vector of regression coefficients
  reg.coefs <- as.vector(unlist(lapply(models,FUN = function(x){return(coef(x)[coef.name])})))
  
  #initialize vector of regression standard errors
  if(is.null(se)){
    #regular standard errors
    se.vector <- unlist(lapply(models,FUN = function(x){return(sqrt(diag(vcov(x)))[coef.name])}))
  }else{
    #user supplied standard errors
    se.vector <- unlist(lapply(se,FUN = function(x){return(x[coef.name])}))
  }
  
  #generate vector of 95% confidence intervals
  reg.lower.95 <- reg.coefs-qnorm(p = 0.975)*se.vector
  reg.upper.95 <- reg.coefs+qnorm(p = 0.975)*se.vector
  
  #generate vector of 90% confidence intervals
  reg.lower.90 <- reg.coefs-qnorm(p = 0.95)*se.vector
  reg.upper.90 <- reg.coefs+qnorm(p = 0.95)*se.vector
  
  #generate dataframe for plotting in ggplot2
  plot.df <- data.frame(cbind(reg.coefs,reg.lower.95,reg.upper.95,reg.lower.90,reg.upper.90,Model=names),row.names = NULL)
  plot.df$reg.coefs <- as.numeric(as.character(plot.df$reg.coefs))
  plot.df$reg.lower.95 <- as.numeric(as.character(plot.df$reg.lower.95))
  plot.df$reg.upper.95 <- as.numeric(as.character(plot.df$reg.upper.95))
  plot.df$reg.lower.90 <- as.numeric(as.character(plot.df$reg.lower.90))
  plot.df$reg.upper.90 <- as.numeric(as.character(plot.df$reg.upper.90))
  
  #plot coefficients
  reg.coef.plot <- plot.df %>% ggplot(aes(x=Model,y=reg.coefs,group=Model,colour=Model)) + 
    geom_point() + 
    geom_hline(yintercept = 0,size=0.5) + 
    geom_linerange(aes(ymax = reg.upper.95,ymin=reg.lower.95))+
    geom_linerange(aes(ymax = reg.upper.90,ymin=reg.lower.90),size=1)+
    ylab('Estimated Coefficient')+
    theme(text = element_text(size=15))
  return(reg.coef.plot)
}


####TIME SERIES####

#tscs.lag <- function(data,unit,var,n=1){
#  library(dplyr,quietly = T)
#  lag.name <- paste('l',var[1],sep='')
#  lagged.data <- data %>% group_by(unit[1]) %>% mutate(lag.name=lag(var,n=n))
#  return(lagged.data)
#}

#lagging time series variables by 1 time unit
#data %>% group_by(unit) %>% mutate(lvar=lag(var))

tslag <- function(x, d=1) {
  x <- as.vector(x)
  n <- length(x)
  c(rep(NA,d),x)[1:n]
}

pastmin <- function(x) {
  xt <- x
  xt[!is.na(xt)] <- cummin(na.omit(x))
  return(xt)
}
pastmax <- function(x) {
  xt <- x
  xt[!is.na(xt)] <- cummax(na.omit(x))
  return(xt)
}

pastsum <- function(x) {
  xt <- x
  xt[!is.na(xt)] <- cumsum(na.omit(x))
  return(xt)
}

pan.lag <- function(x,ind) {
  unlist(tapply(x,ind, function(x) c(NA,x[-length(x)])))
}

pan.lag <- function(x, ind, lag = 1) {
  unlist(tapply(x,ind, function(x) c(rep(NA, times = lag),x[-((length(x) - lag +1):length(x))])))
}


pan.sum <- function(x, ind) {
  unlist(tapply(x,ind, function(x) {
    xt <- x
    xt[!is.na(xt)] <- cumsum(na.omit(x))
    return(xt)
  }))
}

pan.prod <- function(x,ind) {
  unlist(tapply(x,ind, function(x) {
    xt <- x
    xt[!is.na(xt)] <- cumprod(na.omit(x))
    return(xt)
  }))
}

pan.mean <- function(x,ind) {
  unlist(tapply(x, ind, function(x) rep(mean(x, na.rm=TRUE), length(x))))
}

pan.first <- function(x,ind) {
  unlist(tapply(x, ind, function(x) rep(ifelse(any(x),which(x)[1],NA),length(x)) ))
}


pan.min <- function(x,ind) {
  unlist(tapply(x, ind, function(x) rep(min(x, na.rm=TRUE), length(x))))
}


pan.cummin <- function(x, ind) {
  unlist(tapply(x,ind, function(x) {
    xt <- x
    xt[!is.na(xt)] <- cummin(na.omit(x))
    return(xt)
  }))
}

####INSTRUMENTAL VARIABLES###

#get first-stage f-stat
getFStat <- function(ivobject){
  sumiv <- summary(ivobject,diagnostics=T)
  f.stat <- as.numeric(sumiv$diagnostics['Weak instruments','statistic'])
  return(round(f.stat,digits = 2))
}

#iv exclusion sensitivity analysis (Conley et al 2012, ReStat)
#code built off of Jim Bisbee's function (http://static1.squarespace.com/static/559c0c95e4b0c021321c7435/t/56d4a18307eaa0ae80a870eb/1456775555253/handout-2-29.pdf)
#exclusion_sens <- function(ivobject,dat,g){
#  gamma <- seq(-1,1,0.1)
#  newY <- dat$dv - g*dat$instrument
#}
#gamma <- seq(-1,1,.25)
#ExclSens <- function(g) {
#  newY <- dat$logpgp95 - g*dat$logem4
#  coef(ivreg(newY~avexpr+f_brit+f_french,~logem4+f_brit+f_french,cbind(dat,newY)))[2]
#}
#sens.coefs <- sapply(gamma,ExclSens)
#names(sens.coefs)<-gamma
#round(sens.coefs,3)

####DIFF IN DIFF####
#function to create parallel trend plot
parallel.trend <- function(dv,upper,lower,df,treatment){#need to generalize this
  dv <- substitute(dv)
  y.max <- substitute(upper)
  y.min <- substitute(lower)
  y.lab <- NULL
  title <- NULL
  parallel.plot <- ggplot(df,aes(x=year,y=eval(dv),group=treatment,color=factor(treatment)))+
    geom_point()+
    geom_line()+
    geom_vline(xintercept = 2009,size=1.5)+
    geom_linerange(aes(ymax=eval(y.max),ymin=eval(y.min)))
  return(parallel.plot)
}