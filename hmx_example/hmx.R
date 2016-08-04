
## How Much Should We Trust Estimates from Multiplicative
## Interaction Models? Simple Tools to Improve Empirical Practice

## Authors: Jens Hainmueller; Jonathan Mummolo; Yiqing Xu
## This Version: 2016.2.28

######   Interpreting Interaction Models  #######

## 1. inter.raw: first look at the data: D, X, Y
## 2. inter.binning: estimates by X bins
## 3. inter.kernel: kernel smooth plot
## 4. inter.wald: wald test

## Appendix: vcovCluster: clustered standard error

#################################################

inter.rawplot<-function(Y, D, X, data, weights=NULL,
                        Ylabel=NULL, Dlabel=NULL, Xlabel=NULL, 
                        nbins=3, cuts=NULL, span=1, pos=NULL,
                        na.rm = FALSE){ 
    
    ## Y: outcome
    ## D: "treatment" indicator
    ## X: covariate to be interacted with D
    ## nbins: No of bins
    ## cuts: specified cutoff points
    ## span: bandwidth for loess

    ## check missing values
    if (na.rm == TRUE) {
        data <- na.omit(data[,c(Y, D, X)])
    } else {
        if (sum(is.na(data[,Y]))>0) {
            stop(paste("Missing values in ",Y,". Try na.rm=TRUE\n",sep=""))
        }
        if (sum(is.na(data[,D]))>0) {
            stop(paste("Missing values in ",D,". Try na.rm=TRUE\n",sep=""))
        }
        if (sum(is.na(data[,X]))>0) {
            stop(paste("Missing values in ",X,". Try na.rm=TRUE\n",sep=""))
        }
    }
    
    
    ## load packages
    library(ggplot2)
    
    ## variable labels
    if (is.null(Ylabel)==TRUE) {
        Ylabel=Y
    }
    if (is.null(Xlabel)==TRUE) {
        Xlabel=X
    }
    if (is.null(Dlabel)==TRUE) {
        Dlabel=D
    }
    
    
    ## ploting
    if (length(unique(data[,D]))==2) { ## binary case
        
        ## plotting
        D.values <- names(table(data[,D]))
        treat.lab<-c(paste(Dlabel, "=", D.values[1]),
                     paste(Dlabel, "=", D.values[2]))
        data.aug<-data
        data.aug$treat<-factor(data[,D], labels=treat.lab)

        if (is.null(pos)==TRUE) {
            box.pos.co <- min(data[,Y])
            box.pos.tr <- max(data[,Y])
        } else {
            box.pos.co <- pos[1]
            box.pos.tr <- pos[2]
        }

        D.tr <- max(data[,D])
        D.co <- min(data[,D])
        tr <- which(data[,D]==D.tr)
        co <- which(data[,D]==D.co)
        qt90.tr <- quantile(data[tr,X],c(0.05,0.95))
        qt90.co <- quantile(data[co,X],c(0.05,0.95))
        qt50.tr <- quantile(data[tr,X],c(0.25,0.75))
        qt50.co <- quantile(data[co,X],c(0.25,0.75))
        med.tr <- median(data[tr,X])
        med.co <- median(data[co,X])
        
        data.aug$qt90 <- NA
        data.aug$qt90[which(data.aug[,D]==D.tr &
                            data.aug[,X]>=qt90.tr[1] &
                            data.aug[,X]<=qt90.tr[2])]<- box.pos.tr
        data.aug$qt90[which(data.aug[,D]==D.co &
                            data.aug[,X]>=qt90.co[1] &
                            data.aug[,X]<=qt90.co[2])]<- box.pos.co
        data.aug$qt50 <- NA
        data.aug$qt50[which(data.aug[,D]==D.tr &
                            data.aug[,X]>=qt50.tr[1] &
                            data.aug[,X]<=qt50.tr[2])]<- box.pos.tr
        data.aug$qt50[which(data.aug[,D]==D.co &
                            data.aug[,X]>=qt50.co[1] &
                            data.aug[,X]<=qt50.co[2])]<- box.pos.co
        data.aug$med <- NA
        data.aug$med[tr[which.min(abs(data.aug[tr,X]-med.tr))]]<- box.pos.tr
        data.aug$med[co[which.min(abs(data.aug[co,X]-med.co))]]<- box.pos.co
        
        
        ## plotting
        if (is.null(weights)==TRUE) {
            p1 <- ggplot(data.aug, aes_string(X, Y))
        } else {
            p1 <- ggplot(data.aug, aes_string(X, Y, weight=weights))
                
        } 
        p1 <- p1 + geom_point() +
            geom_smooth(method = "lm", se = F, fullrange = T,
                        colour = "steelblue", size = 1)
        
        if (is.null(span)==TRUE) {
            p1 <- p1+ geom_smooth(method = "loess", formula = y ~ x,
                            se = F, colour="red") 
        } else {
            p1 <- p1 + geom_smooth(method = "loess", formula = y ~ x,
                            se = F, colour="red",span=span)
        }
        p1 <- p1 + xlab(Xlabel) + ylab(Ylabel)

       
        p1 <- p1 + geom_line(aes_string(X,"qt90"),
                             size=1,colour="grey50")
        p1 <- p1 + geom_line(aes_string(X,"qt50"),
                             size=3,colour="grey50")
        p1 <- p1 + geom_point(aes_string(X,"med"),size=3,
                              shape=21,fill="white",colour="red")

        p1 <- p1 + theme(axis.title = element_text(size=15))

        p1 <- p1 + facet_wrap(~treat, ncol=1) 
        #p1 <- p1 + facet_grid(treat ~.)                
        suppressWarnings(print(p1))
        
        
    } else { # continuous case
        
        ## grouping by X
        if (is.null(cuts)==TRUE) {
            cutoff<-quantile(data[,X],probs=seq(0,1,1/nbins))
        } else {
            cutoff<-cuts
        }
        while (length(unique(cutoff))!=nbins+1) {
            nbins<-nbins-1; cutoff<-quantile(data[,X],probs=seq(0,1,1/nbins));
        } 
        groupID<-cut(data[,X],breaks=cutoff,label=F)
        groupID[which(data[,X]==min(data[,X]))]<-1
        
        ## X labels
        if (nbins==2) {
            gp.lab<-c(paste(Xlabel,": low",sep=""),paste(Xlabel,": high",sep="")) 
        } else if (nbins==3) {
            gp.lab<-c(paste(Xlabel,": low",sep=""),paste(Xlabel,": medium",sep=""),
                      paste(Xlabel,": high",sep="")) 
        } else {
            gp.lab<-c();for (i in 1:nbins) gp.lab<-c(gp.lab,paste(lab,": Q",i,sep=""))
        } 
        groupID <- factor(groupID, labels=gp.lab)
        data.aug <- data
        data.aug$groupID<-groupID
        
        ## plotting
        if (is.null(weights)==TRUE) {
            p1 <- ggplot(data.aug, aes_string(D, Y))
        } else {
            p1 <- ggplot(data.aug, aes_string(D, Y,weight=weights))
        }
        p1 <- p1 + geom_point() + 
            geom_smooth(method = "lm", se = F, fullrange = T,
                        colour = "steelblue", size = 1)
        if (is.null(span)==TRUE) {
            p1 <- p1 +
                geom_smooth(method = "loess", formula = y ~ x,
                            se = F, colour="red") 
        } else {
            p1 <- p1 +
                geom_smooth(method = "loess", formula = y ~ x,
                            se = F, colour="red",span=span) 
        }
        
        p1 <- p1 + xlab(Dlabel) + ylab(Ylabel) + facet_grid(.~groupID)
        print(p1)
    }
    
}


inter.gam<-function(Y, D, X, Z=NULL,FE=NULL,data, SE=0, k=10,
                    weights=NULL,angle=c(30,100,-30,-120),
                    na.rm = FALSE){


    ## check missing values
    if (na.rm == TRUE) {
        data <- na.omit(data[,c(Y, D, X, Z)])
    } else {
        if (sum(is.na(data[,c(Y, D, X, Z)]))>0) {
            stop("Missing values. Try na.rm=TRUE\n")
        } 
    }

    ## X more than 5 values
    if (length(unique(data[ ,X]))< 5) {
        warning("Moderator less than 5 values; consider a fully saturated model.")
    }
    
    library(mgcv)
    
    if (is.null(FE)==FALSE) {
        if (is.null(Z)==TRUE) {Z<-c()}
        for (i in 1:length(FE)) {
            Z<-c(Z,paste("as.factor(",FE[i],")",sep=""))
        } 
    }
    
    if (is.null(Z)==TRUE) { # no controls
        formula<-as.formula(paste(Y,"~","s(",D,",",X,",k=",k,")"))
    } else {
        formula<-as.formula(paste(Y,"~","s(",D,",",X,",k=",k,")+",
                                  paste(Z,collapse="+")))
    }
    if (is.null(weights)==TRUE) {
        model<-gam(formula, data=data)
    } else {
        model<-gam(formula, data=data,weights=data[,weights])
    }
    par(mfrow=c(2,2),mar=c(2,2,0,0))
    if (SE==0) {
        for (i in angle) {
            vis.gam(model, view=c(D,X), type="response",cex.lab=1.5,zlab=Y,
                    ticktype="detailed",plot.type="persp",
                    n.grid=40,too.far=0.5,theta=i,phi=20)
        }
    } else {
        for (i in angle) {
            vis.gam(model, view=c(D, X), type="response",cex.lab=1.5,zlab=Y,
                    ticktype="detailed",plot.type="persp",
                    se=2,n.grid=40,too.far=0.5,theta=30,phi=20)
        }
    }
}

####################################################

inter.binning<-function(Y,D,X,Z=NULL,FE=NULL,data,weights=NULL,
                        nbins=3, cuts=NULL,
                        vartype = "homoscedastic", cl=NULL,
                        time=NULL,pairwise=TRUE, figure=TRUE,
                        main=NULL, Xlabel=NULL, Dlabel=NULL,
                        Ylabel=NULL, ylim=NULL,interval=NULL,
                        na.rm = FALSE){
    
    ## Y: outcome
    ## D: "treatment" indicator
    ## X: covariate to be interacted with D
    ## Z: control variables
    ## weights: weighting variable
    ## nbins: # of X bins
    ## vartype: homoscedastic (default); "robust"; "cluster", "pcse"
    ## cl: variable to be clustered on
    ## time: time variable for pcse
    ## pairwaise: pcse option

    ## check missing values
    if (na.rm == TRUE) {
        data <- na.omit(data[,c(Y, D, X, Z, cl, time)])
    } else {
        if (sum(is.na(data[,c(Y, D, X, Z, cl, time)]))>0) {
            stop("Missing values. Try na.rm=TRUE\n")
        } 
    }
    
    
    library(ggplot2)
    library(Lmoments)
    N<-dim(data)[1]
    data[,D]<-as.numeric(data[,D])
    
    ## parsing fixed effects
    if (is.null(FE)==FALSE) {
        if (is.null(Z)==TRUE) {Z<-c()}
        for (i in 1:length(FE)) {
            Z<-c(Z,paste("as.factor(",FE[i],")",sep=""))
        } 
    } 

    ## X more than 5 values
    if (length(unique(data[ ,X]))< 5) {
        warning("Moderator less than 5 values; consider a fully saturated model.")
    }

    ## grouping by X
    if (is.null(cuts)==TRUE) {
        cuts.X<-quantile(data[,X],probs=seq(0,1,1/nbins))
    } else {
        cuts.X<-cuts
        nbins<-length(cuts)-1
    }
    while (length(unique(cuts.X))!=nbins+1) {
        nbins<-nbins-1
        cuts.X<-quantile(data[,X],probs=seq(0,1,1/nbins))
    }      
    groupX<-cut(data[,X],breaks=cuts.X,label=F)
    groupX[which(data[,X]==min(data[,X]))]<-1
    
    ## variance of treatment in each group 
    varD <- c()
    obsD <- c()
    for (i in 1:nbins) {
        varD <- c(varD, var(data[groupX==i,D]/mean(data[groupX==i,D])))
        obsD <- c(obsD, length(data[groupX==i,D])/length(data[,D]))
    }
    if (min(obsD) < 0.01) {
        print(list(D.share.byGroup = obsD ))
        stop("One bin has few observations (< 1%).")
    }

    ##############
    ## a naive fit
    if (is.null(Z)==FALSE) {
        mod.f<-as.formula(paste(Y,"~",D,"+",X,"+",D,"*",X,"+",paste(Z,collapse="+"),sep=""))
    } else {
        mod.f<-as.formula(paste(Y,"~",D,"+",X,"+",D,"*",X,sep=""))
    }
    if (is.null(weights)==TRUE) {
        mod.naive<-lm(mod.f,data=data)
    } else {
        mod.naive<-lm(mod.f,data=data,weights=data[,weights])
    }    
    ## coefficients
    coefs<-summary(mod.naive)$coefficients[,1]
    coef.D<-coefs[D]
    coef.X<-coefs[X]
    coef.DX<-coefs[paste(D,X,sep=":")] #interaction
    
    ## # variance
    if (is.null(vartype)==TRUE) {
        vartype <-"homoscedastic"
    }
    if (vartype == "homoscedastic") {
        v<-vcov(mod.naive)
    } else if (vartype == "robust") {
        require(sandwich)
        v<-vcov(mod.naive,type="HC1") # White with small sample correction
    } else if (vartype == "cluster") {
        if (is.null(cl)==FALSE) {
            v<-vcovCluster(mod.naive,cluster = data[,cl])
        } else {  # variable to be clustered not supply
            stop("Warning: clustering variable not found, set cl=varname")
        }
    } else if (vartype =="pcse") {
        library(pcse)
        if (is.null(cl)==FALSE & is.null(time)==FALSE) {
            v<-pcse(mod.naive,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
        } else {  # variable to be clustered not supply
            stop("Warning: please supply unit and time indicators, set cl=varname1, time=varname2")
        }
    }
    
    if (vartype =="pcse") {
        var.D<-v[D,D]
        var.DX<-v[paste(D,X,sep="."),paste(D,X,sep=".")]
        cov.DX<-v[D,paste(D,X,sep=".")]
    } else {
        var.D<-v[D,D]
        var.DX<-v[paste(D,X,sep=":"),paste(D,X,sep=":")]
        cov.DX<-v[D,paste(D,X,sep=":")]
    }
    
    ###  make a vector of the marginal effect of D on Y as X changes 
    X.lvls<-as.numeric(quantile(data[,X], probs=seq(0,1,0.01)))
    marg<-coef.D + coef.DX*X.lvls
    marg
    
    ## the variance is var(B1_D) + X^2*var(B_3) + 2*inst*cov(D, X)
    se<-sqrt(var.D +  X.lvls^2*var.DX + 2*X.lvls*cov.DX)
    df<-mod.naive$df.residual
    crit<-abs(qt(.025, df=df)) # critical values
    
    ##make 95% confidence bands. 
    lb<-marg-crit*se
    ub<-marg+crit*se
    
##################################################
    

    
############## Discre#tize X #################
    
    ## mid points
    x0<-rep(NA,nbins)
    for (i in 1:nbins) x0[i]<-median(data[which(groupX==i),X], na.rm=TRUE)
    
    ## create dummies for bins and interactions
    ## G -- a matrix of group dummies
    ## DG -- a matrix of interactions 
    
    G<-DG<-GX<-DGX<-matrix(0,N,nbins)
    for (i in 1:nbins) {
        G[which(groupX==i),i]<-1
        DG[,i]<-data[,D]*G[,i]
        GX[,i]<-G[,i]*(data[,X]-x0[i])
        DGX[,i]<-DG[,i]*(data[,X]-x0[i])
    }
    
    ## formula and esitmation
    Gs<-GXs<-DGs<-DGXs<-c()
    for (i in 1:nbins)  {
        Gs<-c(Gs,paste("G[,",i,"]",sep=""))
        GXs<-c(GXs,paste("GX[,",i,"]",sep=""))
        DGs<-c(DGs,paste("DG[,",i,"]",sep=""))
        DGXs<-c(DGXs,paste("DGX[,",i,"]",sep=""))
    }
    Xf<-paste(Y,"~ -1+",paste(DGs,collapse="+"),"+",paste(DGXs,collapse="+"),
              "+",paste(Gs,collapse="+"),"+",paste(GXs,collapse="+"),sep="")
    if (is.null(Z)==FALSE) {
        Xf<-paste(Xf,"+",paste(Z,collapse="+"),sep="")
    }
    mod.Xf<-as.formula(Xf)
    if (is.null(weights)==TRUE) {
        mod.X<-lm(mod.Xf,data=data) 
    } else {
        mod.X<-lm(mod.Xf,data=data,weights=data[,weights])
    }    
    
    ## coefficients and CIs
    Xcoefs<-mod.X$coefficients[1:nbins]
    if (vartype == "homoscedastic") {
        X.v<-vcov(mod.X)
    } else if (vartype == "robust") {
        X.v<-vcov(mod.X,type="HC1") ## White with small sample correction
    } else if (vartype == "cluster") {
        X.v<-vcovCluster(mod.X,cluster=data[,cl])
    } else if (vartype == "pcse") {
        if (is.null(Z)==FALSE) {
            exclude<-names(which(is.na(mod.X$coefficients)==TRUE))  ## drop colinear variables
            Z.ex<-setdiff(Z,exclude)
            Xf<-paste(Y,"~ -1+",paste(DGs,collapse="+"),"+",paste(DGXs,collapse="+"),
                      "+",paste(Gs,collapse="+"),"+",paste(GXs,collapse="+"),"+",paste(Z.ex,collapse="+"),sep="")
            mod.X<-lm(as.formula(Xf),data=data)
        }
        X.v<-pcse(mod.X,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
    }
    X.v<-X.v[1:nbins,1:nbins] 
    X.se<-sqrt(diag(X.v)) 
    df.X<-mod.X$df.residual
    crit.X<-abs(qt(.025, df=df.X))
    lb.X<-Xcoefs-crit.X*X.se
    ub.X<-Xcoefs+crit.X*X.se

   
################### plotting ###############################
    
    ## margin and label adjustment
    
    if (figure==TRUE) {
        
        if(is.null(Xlabel)==FALSE){
            x.label<-c(paste("Moderator: ", Xlabel, sep=""))
            y.label<-c(paste("Marginal effect of ",Dlabel," on ",Ylabel,sep=""))
        } else {
            x.label<-c(paste("Moderator: ", X, sep=""))
            y.label<-c(paste("Marginal effect of ",D," on ",Y,sep=""))
        }
        out<-data.frame(X.lvls,marg,lb,ub)
        out.bin<-data.frame(x0,Xcoefs,lb.X,ub.X)
        out.bin2<-out.bin[which(is.na(Xcoefs)==FALSE),] ## non missing part
        out.bin3<-out.bin[which(is.na(Xcoefs)==TRUE),]  ## missing part
        errorbar.width<-(max(X.lvls)-min(X.lvls))/20
        
        p1 <- ggplot()
        yrange<-na.omit(c(marg,lb,ub,Xcoefs,lb.X,ub.X))
        if (is.null(ylim)==FALSE) {yrange<-c(ylim[2],ylim[1]+(ylim[2]-ylim[1])*1/6)}
        maxdiff<-(max(yrange)-min(yrange))
        pos<-max(yrange)-maxdiff/20

        ## mark zero
        p1 <- p1 + geom_hline(yintercept=0,colour="white",size=2)

        ## mark the original interval
        if (is.null(interval)==FALSE) {
            p1<- p1 + geom_vline(xintercept=interval,colour="steelblue",linetype=2)
        }
        
        ## histogram
        if (length(unique(data[,D]))==2) { ## binary D

            if (is.null(weights)==TRUE) {
                hist.out<-hist(data[,X],breaks=80,plot=FALSE)
            } else {
                suppressWarnings(library(plotrix))
                suppressWarnings(hist.out<-hist(data[,X],data[,weights],
                                               breaks=80,plot=FALSE))
            } 
            n.hist<-length(hist.out$mids)
            dist<-hist.out$mids[2]-hist.out$mids[1]
            hist.max<-max(hist.out$counts)
            ## count the number of treated
            count1<-rep(0,n.hist)
            treat<-which(data[,D]==max(data[,D]))
            for (i in 1:n.hist) {
                count1[i]<-sum(data[treat,X]>=hist.out$breaks[i] &
                               data[treat,X]<hist.out$breaks[(i+1)])
            }
            count1[n.hist]<-sum(data[treat,X]>=hist.out$breaks[n.hist] &
                                data[treat,X]<=hist.out$breaks[n.hist+1])
            ## put in a data frame
            histX<-data.frame(ymin=rep(min(yrange)-maxdiff/5,n.hist),
                              ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                              xmin=hist.out$mids-dist/2,
                              xmax=hist.out$mids+dist/2,
                              count1=count1/hist.max*maxdiff/5+min(yrange)-maxdiff/5)
            p1 <- p1 + geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                                 colour="gray50",alpha=0,size=0.5) + # control
                geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=count1),
                          fill="red",colour="grey50",alpha=0.3,size=0.5) # treated
            
        } else { ## continuous D
            hist.out<-hist(data[,X],breaks=80,plot=FALSE)
            n.hist<-length(hist.out$mids)
            dist<-hist.out$mids[2]-hist.out$mids[1]
            hist.max<-max(hist.out$counts)
            
            histX<-data.frame(ymin=rep(min(yrange)-maxdiff/5,n.hist),
                              ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                              xmin=hist.out$mids-dist/2,
                              xmax=hist.out$mids+dist/2)
            p1 <- p1 + geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                                 colour="gray50",alpha=0,size=0.5)
        }
        
        if (is.null(main)==FALSE) {
            p1<-p1+ggtitle(main)+
                theme(plot.title = element_text(size=35,lineheight=.8, face="bold"))
        } 

      
        ## labels: L, M, H
        if (nbins==3) {
            p1<-p1 + annotate(geom="text", x=out.bin[1,1], y=pos, label="L",colour="gray50",size=10) +
                annotate(geom="text", x=out.bin[2,1], y=pos, label="M",colour="gray50",size=10) +
                    annotate(geom="text", x=out.bin[3,1], y=pos, label="H",colour="gray50",size=10)
        } else if (nbins==4) {
            p1<-p1 + annotate(geom="text", x=out.bin[1,1], y=pos, label="L",colour="gray50",size=10) +
                annotate(geom="text", x=out.bin[2,1], y=pos, label="M1",colour="gray50",size=10) +
                    annotate(geom="text", x=out.bin[3,1], y=pos, label="M2",colour="gray50",size=10) +
                        annotate(geom="text", x=out.bin[4,1], y=pos, label="H",colour="gray50",size=10)
        }
        
        ## linear plot 
        p1<-p1 + geom_line(data=out,aes(X.lvls,marg))+
            geom_ribbon(data=out, aes(x=X.lvls,ymin=lb,ymax=ub),alpha=0.2)+
                xlab(x.label) + ylab(y.label)
        ## bin estimates
        p1<-p1+ geom_errorbar(data=out.bin2, aes(x=x0, ymin=lb.X, ymax=ub.X),colour="red",
                              width= errorbar.width)+
                                  geom_point(data=out.bin2,aes(x0,Xcoefs),size=4,shape=21,fill="white",colour="red") +
                                      theme(axis.title = element_text(size=15))
        ## in case there's non-overlap
        p1<-p1+annotate(geom="text", x=out.bin3[,1], y=rep(0,dim(out.bin3)[1]),
                        label="NaN",colour="red") 

        if (is.null(ylim)==FALSE) {p1<-p1+coord_cartesian(
                                              ylim = c(ylim[1]-(ylim[2]-ylim[1])*0.25/6,
                                                       ylim[2]+(ylim[2]-ylim[1])*0.4/6))}
        
        plot(p1)

        
    }  # end of plotting

    ## ################# testing  ###############################

    ## binary treatment
    btreat<-length(unique(data[,D]))==2
    
    
    
    
    ## if the signs are correct
    correctOrder<-ifelse(as.numeric((Xcoefs[1]-Xcoefs[2])*(Xcoefs[2]-Xcoefs[3]))>0,TRUE,FALSE) 
    
    ## p values
    pvalue<-function(i,j){
        stat<-(Xcoefs[i]-Xcoefs[j])/sqrt(X.v[i,i]+X.v[j,j]-2*X.v[i,j])
        p<-(1-pt(abs(stat),df.X))*2
        return(p)
    }
    if (nbins==3) {
        p.twosided<-round(c(pvalue(1,2),pvalue(2,3),pvalue(1,3)),digit=4)
        names(p.twosided)<-c("p.1v2","p.2v3","p.1v3")
        names(Xcoefs)<-c("X_low","X_med","X_high")
    } else if (nbins==2) {
        p.twosided<-round(pvalue(1,2),digit=4)
        names(p.twosided)<-c("p.LvH")
        names(Xcoefs)<-c("X_low","X_high")
    } else if (nbins==4) {
        names(Xcoefs)<-c("X_low","X_med1","X_med2","X_high")
    }

    ## L-kurtosis of the moderator
    lkurtosis<-Lmoments(data[,X],returnobject=TRUE)$ratios[4]
 
    ## storage
    out<-list("D.binary" = btreat,
              "D.share.byGroup" = round(obsD, 3),
              "D.variation.byGroup" = round(varD,3),
              "lkurtosis" = round(lkurtosis,3)
              )
    if (nbins==3) {
        out<-c(out,list(correctOrder=correctOrder))
    }
    out<-c(out,list(coefs=round(Xcoefs,4)))
    if (nbins%in%c(2,3)) {
        out<-c(out,list(p.twosided=p.twosided)) 
    }
    return(out)
}

###########################

inter.kernel<-function(Y,D,X,Z=NULL,FE=NULL,data,nboot=500,grid=30,h=NULL,cl=NULL,
                       Dlabel=NULL,Xlabel=NULL,Ylabel=NULL,
                       cores=4,parallel=TRUE,ylim=NULL,file=NULL,
                       na.rm = FALSE){
    
    ## nboot: number of bootstraps
    ## grid: points at which kernel regressions are estimated
    ## h: bandwidth
    ## cl: clustering variable

    ## check missing values
    if (na.rm == TRUE) {
        data <- na.omit(data[,c(Y, D, X, Z, cl)])
    } else {
        if (sum(is.na(data[,c(Y, D, X, Z, cl)]))>0) {
            stop("Missing values. Try na.rm=TRUE\n")
        }
    }

    ## X more than 5 values
    if (length(unique(data[ ,X]))< 5) {
        warning("Moderator less than 5 values; consider a fully saturated model.")
    }
    
    library(ggplot2)
    library(np)
    
    ## parpare
    n<-dim(data)[1]
    newx<-seq(min(data[,X]),max(data[,X]),length=grid) ## grid: at which regressions are run
    newy<-newd<-rep(0,grid)
    newz<-matrix(0,grid,length(Z))
    
    ## partial out fixed effects
    if (is.null(FE)==FALSE) {
        fe<-c()
        for (i in 1:length(FE)) {
            fe<-c(fe,paste("as.factor(",FE[i],")",sep=""))
        }
        Y.tilde<-lm(paste(Y,"~",paste(fe,collapse="+")),data=data)$residuals
        D.tilde<-lm(paste(D,"~",paste(fe,collapse="+")),data=data)$residuals        
    } else {
        Y.tilde<-data[,Y]
        D.tilde<-data[,D]
    }
    DX.tilde<-D.tilde*data[,X]
    
    ## bandwidth selection
    if (is.null(h)==TRUE) {
        bw<-npscoefbw(ydat=Y.tilde, xdat=cbind(D.tilde,DX.tilde,
                                               data[,c(X,Z)]),zdat=data[,X])
        h.max<-(max(data[,X])-min(data[,X]))/5
        if (h.max<bw$bw) {
            bw<-npscoefbw(ydat=Y.tilde, xdat=cbind(D.tilde,DX.tilde, data[,c(X,Z)]),zdat=data[,X],bws=h.max,
                          bandwidth.compute=FALSE)
        } 
    } else { ## set bandwidth mannually
        bw<-npscoefbw(ydat=Y.tilde, xdat=cbind(D.tilde,DX.tilde,data[,c(X,Z)]),zdat=data[,X],
                      bws=h,bandwidth.compute=FALSE)
        
    }
    evaluate.data<-cbind(newd,newd,newx,newz)
    mod.new<-npscoef(bw,betas=TRUE,exdat=evaluate.data,ezdat=newx)
    coefs<-coef(mod.new)
    coef<-coefs[,2]+coefs[,3]*newx
    cat("Bandwidth =", round(mod.new$bw,3),"\n")
    
    ## bootstrap standard errors and CI
    if (is.null(cl)==FALSE) { ## find clusters
        clusters<-unique(data[,cl])
        id.list<-split(1:n,data[,cl])
    }
    
    
    ## the function
    getCoef<-function(i) {
        if (is.null(cl)==TRUE) {
            smp<-sample(1:n,n,replace=TRUE)
        } else { ## block bootstrap
            cluster.boot<-sample(clusters,length(clusters),replace=TRUE)
            smp<-unlist(id.list[match(cluster.boot,clusters)])
        }   
        s<-data[smp,]
        if (is.null(FE)==FALSE) {
            s.Y<-lm(paste(Y,"~",paste(fe,collapse="+")),data=s)$residuals
            s.D<-lm(paste(D,"~",paste(fe,collapse="+")),data=s)$residuals            
        } else {
            s.Y<-s[,Y]
            s.D<-s[,D]
        }
        s.DX<-s.D*s[,X]
        mod.boot<-npscoef(bw,tydat=s.Y, txdat=cbind(s.D,s.DX,s[,c(X,Z)]), tzdat=s[,X],
                          exdat=evaluate.data, ezdat=newx,betas=TRUE)
        coefs<-coef(mod.boot)
        out<-coefs[,2]+coefs[,3]*newx
        return(out)
    }
    
    ## start loop
    cat("Bootstrapping:\n")
    if (parallel==TRUE) {
        library(foreach)
        library(doParallel)
        pcl<-makeCluster(cores)  
        registerDoParallel(pcl)
        cat("parallel computing with", cores,"cores...")
        suppressWarnings(
            coefs<-foreach (i=1:nboot, .combine=cbind, .packages="np",
                            .export="getCoef",.inorder=FALSE) %dopar% {getCoef(i)}
        )
        stopCluster(pcl)
        cat("\n") 
    } else {
        coefs<-matrix(NA,grid,nboot)
        for (i in 1:nboot) {
            coefs[,i]<-getCoef(i)
            if (i%%50==0) cat(i) else cat(".")
        }
        cat("\n")        
    }
    
    ## summary
    CI<-t(apply(coefs,1,quantile,c(0.025,0.975)))
    est<-data.frame(cbind("X"=newx,"Coef"=coef,"SE"=apply(coefs,1,sd),
                          "CI_lower"=CI[,1], "CI_upper"=CI[,2]))
    
    
    ## plotting
    if(is.null(Xlabel)==FALSE){
        x.label<-c(paste("Moderator: ", Xlabel, sep=""))
        y.label<-c(paste("Marginal effect of ",Dlabel," on ",Ylabel,sep=""))
    } else {
        x.label<-c(paste("Moderator: ", X, sep=""))
        y.label<-c(paste("Marginal effect of ",D," on ",Y,sep=""))
    }
    ##main<-paste("Gaussian kernel, bandwidth =",signif(bw$bw,3))

    ## mark zero
    p1 <- ggplot() + geom_hline(yintercept=0,colour="white",size=2)
    
    p1 <-  p1 + geom_line(data=est,aes(X,Coef))+
        geom_ribbon(data=est, aes(x=X,ymin=CI_lower,ymax=CI_upper),alpha=0.2)+
        xlab(x.label) + ylab(y.label) + theme(axis.title = element_text(size=15))

    
    ## histogram
    yrange<-na.omit(c(CI))
    if (is.null(ylim)==FALSE) {yrange<-c(ylim[2],ylim[1]+(ylim[2]-ylim[1])*1/6)}
    maxdiff<-(max(yrange)-min(yrange))

    if (length(unique(data[,D]))==2) { ## binary D
        
        hist.out<-hist(data[,X],breaks=80,plot=FALSE)
        n.hist<-length(hist.out$mids)
        dist<-hist.out$mids[2]-hist.out$mids[1]
        hist.max<-max(hist.out$counts)
        ## count the number of treated
        count1<-rep(0,n.hist)
        treat<-which(data[,D]==max(data[,D]))
        for (i in 1:n.hist) {
            count1[i]<-sum(data[treat,X]>=hist.out$breaks[i] &
                           data[treat,X]<hist.out$breaks[(i+1)])
        }
        count1[n.hist]<-sum(data[treat,X]>=hist.out$breaks[n.hist] &
                            data[treat,X]<=hist.out$breaks[n.hist+1])
        ## put in a data frame
        histX<-data.frame(ymin=rep(min(yrange)-maxdiff/5,n.hist),
                          ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                          xmin=hist.out$mids-dist/2,
                          xmax=hist.out$mids+dist/2,
                          count1=count1/hist.max*maxdiff/5+min(yrange)-maxdiff/5)
        p1 <- p1 + geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                             colour="gray50",alpha=0,size=0.5) + # control
            geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=count1),
                      fill="red",colour="grey50",alpha=0.3,size=0.5) # treated
        
    } else { ## continuous D
        hist.out<-hist(data[,X],breaks=80,plot=FALSE)
        n.hist<-length(hist.out$mids)
        dist<-hist.out$mids[2]-hist.out$mids[1]
        hist.max<-max(hist.out$counts)
        
        histX<-data.frame(ymin=rep(min(yrange)-maxdiff/5,n.hist),
                          ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                          xmin=hist.out$mids-dist/2,
                          xmax=hist.out$mids+dist/2)
        p1 <- p1 + geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                             colour="gray50",alpha=0,size=0.5)
    }
    
    if (is.null(ylim)==FALSE) {p1<-p1+coord_cartesian(
                                          ylim = c(ylim[1]-(ylim[2]-ylim[1])*0.25/6,
                                                   ylim[2]+(ylim[2]-ylim[1])*0.4/6))}
    if (is.null(file)==FALSE) {
        pdf(file)
        plot(p1)
        graphics.off() 
    } else {
        plot(p1)
    }
    
    output<-list(bw=bw$bw,kernel=bw$ckertype,est=est)
    return(output)
}



#############
inter.wald<-function(Y,D,X,Z=NULL,FE=NULL, data,
                     vartype = NULL, cl = NULL, 
                     pairwise = TRUE,  time=NULL,
                     weights=NULL, nbins=3, cuts=NULL,
                     na.rm = FALSE){
    
    ## Y: outcome
    ## D: "treatment" indicator
    ## X: covariate to be interacted with D
    ## Z: control variables
    ## weights: weighting variable
    ## nbins: # of X bins
    ## var: homoscedastic (default); "robust"; "cluster", "pcse"
    ## cl: variable to be clustered on
    ## time: time variable for pcse
    ## pairwaise: pcse option

    ## check missing values
    if (na.rm == TRUE) {
        data <- na.omit(data[,c(Y, D, X, Z, cl, time)])
    } else {
        if (sum(is.na(data[,c(Y, D, X, Z, cl, time)]))>0) {
            stop("Missing values. Try na.rm=TRUE\n")
        }
    }
    
    N<-dim(data)[1]
    data[,D]<-as.numeric(data[,D])
    
    ## parsing fixed effects
    if (is.null(FE)==FALSE) {
        if (is.null(Z)==TRUE) {Z<-c()}
        for (i in 1:length(FE)) {
            Z<-c(Z,paste("as.factor(",FE[i],")",sep=""))
        } 
    }
    
    ## formula
    formula0 <- paste(Y,"~",D,"+",X,"+",D,"*",X)
     
    
    ## ############ Binning #########################
    
    ## grouping by X
    if (is.null(cuts)==TRUE) {
        cuts.X<-quantile(data[,X],probs=seq(0,1,1/nbins))
    } else {
        cuts.X<-cuts
        nbins<-length(cuts)-1
    }
    while (length(unique(cuts.X))!=nbins+1) {
        nbins<-nbins-1; cuts.X<-quantile(data[,X],probs=seq(0,1,1/nbins));
    }      
    groupX<-cut(data[,X],breaks=cuts.X,label=F)
    groupX[which(data[,X]==min(data[,X]))]<-1

    ## create dummies for bins and interactions
    ## G -- a matrix of group dummies
    ## DG -- a matrix of interactions
    
    ## mid points
    x0<-rep(NA,nbins)
    for (i in 1:nbins) x0[i]<-median(data[which(groupX==i),X], na.rm=TRUE)
     
    
    G<-DG<-GX<-DGX<-matrix(0,N,(nbins-1))
    for (i in 1:(nbins-1)) {
        G[which(groupX==(i+1)),i]<-1
        DG[,i]<-data[,D]*G[,i]
        GX[,i]<-data[,X]*G[,i]
        DGX[,i]<-data[,D]*data[,X]*G[,i]
    }
    
    ## formula and esitmation
    Gs<-GXs<-DGs<-DGXs<-c()
    for (i in 2:nbins)  {
        Gs<-c(Gs,paste("G",i,sep=""))
        GXs<-c(GXs,paste("GX",i,sep=""))
        DGs<-c(DGs,paste("DG",i,sep=""))
        DGXs<-c(DGXs,paste("DGX",i,sep=""))
    }
    colnames(G) <- Gs
    colnames(DG) <- DGs
    colnames(GX) <- GXs
    colnames(DGX) <- DGXs
    
    data.aug <- cbind.data.frame(data, G, DG, GX, DGX)
    formula1<-paste(formula0,
              "+",paste(Gs,collapse=" + "),
              "+",paste(GXs,collapse=" + "),
              "+",paste(DGs,collapse=" + "),
              "+",paste(DGXs,collapse=" + "))

    if (is.null(Z)==FALSE) {
        formula0 <- paste(formula0, "+",paste(Z,collapse=" + "))
        formula1 <- paste(formula1, "+",paste(Z,collapse=" + "))
    } 
    
    ##### fit the two models
 
    if (is.null(weights)==TRUE) {
        mod.re<-lm(as.formula(formula0),data=data.aug) 
        mod.un<-lm(as.formula(formula1), data=data.aug) 
    } else {
        mod.re<-lm(as.formula(formula0), data=data.aug, weights=data.aug[,weights])
        mod.un<-lm(as.formula(formula1), data=data.aug, weights=data.aug[,weights])
    }

    ## vcov
    if (is.null(vartype)==TRUE) {
        vartype <- "homoscedastic"
    }
    if (vartype=="homoscedastic") {
        v<-vcov(mod.un)
    } else if (vartype=="robust") {
        require(sandwich)
        v<-vcov(mod.un,type="HC1") # White with small sample correction
    } else if (vartype=="cluster") {
        if (is.null(cl)==FALSE) {
            v<-vcovCluster(mod.un,cluster = data.aug[,cl])
        } else {  # variable to be clustered not supply
            stop("Warning: clustering variable not found, set cl=varname")
        }
    } else if (vartype=="pcse") {
        library(pcse)
        if (is.null(cl)==FALSE & is.null(time)==FALSE) {
            v<-pcse(mod.un,
                    groupN=data.aug[,cl],
                    groupT=data.aug[,time],
                    pairwise=pairwise)$vcov
        } else {  # variable to be clustered not supply
            stop("Warning: please supply unit and time indicators, set cl=varname1, time=varname2")
        }
    }
    library(lmtest)
    wald <- waldtest(mod.re, mod.un, test="Chisq", vcov=v)
    p.wald <- wald[[4]][2]
    print(p.wald)
     
    out <- list(wald = wald)

    return(out)
    
    
}

#############

## vcovCluster.r 
## function to compute var-cov matrix using clustered robust standard errors
## inputs:
## model = model object from call to lm or glm
## cluster = vector with cluster ID indicators
## output:
## cluster robust var-cov matrix
## to call this for a model directly use:
## coeftest(model,vcov = vcovCluster(model, cluster))
## formula is similar to Stata's , cluster command

vcovCluster <- function(
    model,
    cluster
)
{
    require(sandwich)
    ## require(lmtest)
    
    if(nrow(model.matrix(model))!=length(cluster)){
        stop("check your data: cluster variable has different N than model")
    }
    M <- length(unique(cluster))
    N <- length(cluster)           
    K <- model$rank   
    if(M<50){
        warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
    }
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
    rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
    ##colnames(rcse.cov)<-rownames(rcse.cov)<-names(model$coefficients)
    return(rcse.cov)
}



