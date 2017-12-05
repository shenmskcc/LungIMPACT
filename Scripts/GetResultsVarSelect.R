##########################################
###### AUTOMATIZE READING RESULTS ########
##########################################

## INPUT : 
#         1) studyType : name of the IMPACT directory where the results/data are saved
#         2) method : Method implemented to be analysed (LASSO,RF,GBM,CF)
#         3) LT : is the data left truncated : TRUE v FALSE

## OUTPUT :
#         For each of the methods to be analyzed will give 
#         1) CI distribution bar plot
#         2) An influence plot of top hit genes
#         3) Kaplan meier plot based on 4 risk groups
#         4) Proportion of mutated genes in each of the 4 risk groups, per top hit gene

getResults <- function(studyType,method,CNV,OnlyCNV,geneList){
  
  if(CNV){studyType2 <- paste0(studyType,"CNV")}
  if(OnlyCNV){studyType2 <- paste0(studyType,"OnlyCNV")}
  if(!CNV && !OnlyCNV){studyType2 <- studyType}
  
  if("LASSO" %in% method) try(load(paste0("./Study/",studyType,"/results/",studyType2,"_LASSO.Rdata")))
  if("CF" %in% method) try(load(paste0("./Study/",studyType,"/results/",studyType2,"_CF.Rdata")))
  data <- read.csv(paste0("./Study/",studyType,"/data/",studyType2,"ReadyStudy.csv"),header = TRUE, row.names = 1)
  
  ## make clinical names for LUNG
  
  if(!CNV && !OnlyCNV) data.clin <- read.csv(paste0("./Study/",studyType,"/data/",studyType,"ClinReadyStudy.csv"),header = TRUE, row.names = 1)
  if(CNV) data.clin <- read.csv(paste0("./Study/",studyType,"/data/",studyType,"ClinCNVReadyStudy.csv"),header = TRUE, row.names = 1)
  if(OnlyCNV) data.clin <- read.csv(paste0("./Study/",studyType,"/data/",studyType,"ClinOnlyCNVReadyStudy.csv"),header = TRUE, row.names = 1)
  clin.names <- colnames(data.clin[,which(is.na(match(colnames(data.clin),colnames(data))))]) 
  data.clin[,match(clin.names,colnames(data.clin))] <- data.clin[,match(clin.names,colnames(data.clin))] -1
  
  ## determine if left truncated
  if(colnames(data)[1]  == "time") {LT = FALSE}
  if(colnames(data)[1] == "time1") {LT = TRUE}
  
  if( LT && max(data$time2) > 1000){MD = 365
  time.type = "Days"}
  if( LT && max(data$time2) < 1000){MD = 12
  time.type = "Months"}
  if( !LT && max(data$time) > 1000){MD = 365
  time.type = "Days"}
  if( !LT && max(data$time) < 1000){MD = 12
  time.type = "Months"}
  ### LASSO ANALYSIS ###
  
  #try(setwd("./results/"))
  
  if("LASSO" %in% method) {
    
    if(CNV && studyType == "Lung"){LASSO <- ENET2}
    # ConcordanceIndex <- as.data.frame(as.vector(unlist(sapply(LASSO, "[[", "CI"))))
    # summary.CI <- round(as.data.frame(c(quantile(ConcordanceIndex[,1],c(0.1,0.25,0.5,0.75,0.9)))),digits = 2)
    # colnames(summary.CI) <- "Concordance Index"
    # rownames(summary.CI) <- c("Lower 10%","1st Quarter","Median","3rd Quarter","Upper 10%")
    # CI.BP <- as.data.frame(t(summary.CI))
    # ### FOR APP
    # # CI.BP <- ggplot(ConcordanceIndex, aes(x="",y = ConcordanceIndex)) +
    # #   geom_boxplot(outlier.colour = "red", outlier.shape = 1,colour = "#3366FF") +
    # #   geom_jitter(width = 0.2) +
    # #   theme_grey() +
    # #   labs(title = "Concordance Index distribution", subtitle = "Using genetic data only")
    # ###
    # 
    # if(LT){genes <- colnames(data)[-c(1:3)]}
    # if(!LT){genes <- colnames(data)[-c(1:2)]}
    # selected.genes.lasso <- as.data.frame(matrix(0L,nrow=1,ncol = length(genes)))
    # colnames(selected.genes.lasso) <- genes
    # 
    # for(i in 1:length(LASSO)){
    #   temp <-LASSO[[i]]$fit
    #   if(!is.null(temp)){
    #     selected.temp <- rownames(temp)
    #     selected.genes.lasso[,match(selected.temp,colnames(selected.genes.lasso))] <-
    #       selected.genes.lasso[,match(selected.temp,colnames(selected.genes.lasso))] +1
    #   }
    # }
    # 
    # topHits <- names(sort(selected.genes.lasso,decreasing = TRUE))[1:10]
    # if(length(selected.genes.lasso)  >= 20){ranked.lasso <- sort(selected.genes.lasso,decreasing = TRUE)[1:20]}
    # if(length(selected.genes.lasso)  < 20){ranked.lasso <- sort(selected.genes.lasso,decreasing = TRUE)[1:length(selected.genes.lasso)]}
    # melt.rank.lasso <- melt(ranked.lasso)
    # melt.rank.lasso$value <- melt.rank.lasso$value/200
    # colnames(melt.rank.lasso) <- c("Gene","Frequency")
    # 
    # influencePlot <- ggplot(melt.rank.lasso,aes(x=Gene,y=Frequency,fill=Gene))+geom_col()+
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #   #scale_fill_discrete(guide = guide_legend(title = "Category")) +
    #   #scale_fill_manual(values=rainbow(length(unique(melt.rank.lasso$Gene)))) +
    #   labs(title = "Most recurrent selected genes out of 200 runs") +
    #   theme(legend.position="none")

    # save(influencePlot,file="./Study/Lung/results/influencePlot.Rdata")
    # 
    if(LT) Variables <- colnames(data)[-c(1:3)]
    if(!LT) Variables <- colnames(data)[-c(1:2)]
    allCoefs <- as.data.frame(matrix(nrow=length(LASSO),ncol=length(Variables)))
    colnames(allCoefs) <- Variables
    
    for(x in 1:length(LASSO)){
      coefsValues <- LASSO[[x]]$fit[,1]
      allCoefs[x,match(names(coefsValues),colnames(allCoefs))] <- as.numeric(coefsValues)
    }
    allCoefs[is.na(allCoefs)] <- 0
    
    meanCoefs <- apply(allCoefs,2,function(x){mean(x,na.rm = TRUE)})
    selectFreq <- apply(allCoefs,2,function(x){
      length(which(x!=0))/length(x)
    })
    
    ## get mu freq
    if(LT) data.temp <- data[,-c(1:3)]
    if(!LT) data.temp <- data[,-c(1:2)]
    MutationFrequency <- apply(data.temp,2,function(x){
      sum(x)/length(x)
    })
    
    resultsAll <- as.data.frame(cbind(meanCoefs,selectFreq,MutationFrequency))
    colnames(resultsAll) <- c("MeanCoefficient","SelectionFrequency","MutationFrequency")
    rownames(resultsAll) <- names(meanCoefs)
    resultsAll <- resultsAll[complete.cases(resultsAll),]
    resultsAll$GeneName <- rownames(resultsAll)
    resultsAll$MutationFrequency2 <- cut(resultsAll$MutationFrequency, c(0,0.10,0.20,0.40))
    
    if(length(geneList)!=0){
      m <- resultsAll[match(geneList,rownames(resultsAll)), ]
      
      a <- list(
        x = m$MeanCoefficient,
        y = m$SelectionFrequency,
        text = rownames(m),
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 7,
        ax = 20,
        ay = -40
      )
      
      selectInflPlot <- plot_ly(data = resultsAll, x = ~MeanCoefficient, y = ~SelectionFrequency,
                                text = ~paste('Gene :',GeneName,
                                              '</br> Hazard Ratio :',round(exp(MeanCoefficient),digits=2)),
                                mode = "markers",size = ~MutationFrequency,color = ~MutationFrequency) %>% 
        layout(title ="Volcano Plot",annotations = a)
    }
    
    else{
      selectInflPlot <- plot_ly(data = resultsAll, x = ~MeanCoefficient, y = ~SelectionFrequency,
                                text = ~paste('Gene :',GeneName,
                                              '</br> Hazard Ratio :',round(exp(MeanCoefficient),digits=2)),
                                mode = "markers",size = ~MutationFrequency,color = ~MutationFrequency) %>% 
        layout(title ="Volcano Plot")
    }
    
    # final.pred <- as.data.frame(matrix(nrow= nrow(data),ncol = length(LASSO)))
    # rownames(final.pred) <- rownames(data)
    # for(i in 1:length(LASSO)){
    #   temp <- LASSO[[i]]$predicted
    #   final.pred[match(names(temp),rownames(final.pred)),i] <- as.numeric(temp)
    # }
    # 
    # average.risk <- apply(final.pred,1,function(x){
    #   mean(as.numeric(x),na.rm = TRUE)
    # })
    # RiskScore <- rescale(average.risk, to = c(0, 10), from = range(average.risk, na.rm = TRUE, finite = TRUE))
    # summary.RiskScore <- round(as.data.frame(c(quantile(RiskScore,c(0.1,0.25,0.5,0.75,0.9)))),digits = 2)
    # colnames(summary.RiskScore) <- "Risk Score"
    # rownames(summary.RiskScore) <- c("Lower 10%","1st Quarter","Median","3rd Quarter","Upper 10%")
    # ## refit coxph model with average risk as covariate
    # meanRS <- mean(RiskScore)
    # #RiskScore <- average.risk #- meanRS
    # if(LT) refit.risk <- coxph(Surv(data$time1,data$time2,data$status)~RiskScore)
    # if(!LT) refit.risk <- coxph(Surv(data$time,data$status)~RiskScore)
    # Risk <- as.data.frame(RiskScore)
    # 
    # RiskHistogram <- ggplot(Risk, aes(x = RiskScore, y = ..density..)) +
    #   geom_histogram(show.legend = FALSE, aes(fill=..x..),
    #                  breaks=seq(min(Risk$RiskScore), max(Risk$RiskScore), by=0.05)) +
    #   geom_density(show.legend = FALSE) +
    #   theme_minimal() +
    #   labs(x = "Average risk score", y = "Density") +
    #   scale_fill_gradient(high = "red", low = "green")
    # 
    # # save(RiskHistogram,file="./Study/Lung/results/RiskHistogram.Rdata")
    # 
    # ## Refit with clinical
    # #1) Match patients
    # average.risk.clin <- RiskScore[match(rownames(data.clin),names(RiskScore))]
    # if(LT) data.refit.clin <- cbind(data.clin[,match(c("time1","time2","status",clin.names),colnames(data.clin))],average.risk.clin)
    # if(!LT) data.refit.clin <- cbind(data.clin[,match(c("time","status",clin.names),colnames(data.clin))],average.risk.clin)
    # colnames(data.refit.clin)[ncol(data.refit.clin)] <- "RiskScore"
    # # refit
    # if(LT) refit.risk.clin <- coxph(Surv(time1,time2,status)~.,data=data.refit.clin)
    # if(!LT) refit.risk.clin <- coxph(Surv(data$time,data$status)~.,data=data.refit.clin)
    # refit.risk <- as.data.frame(summary(refit.risk)$coefficients)
    # refit.risk$CI <- paste0("(",round(exp(refit.risk$coef-1.96*refit.risk$`se(coef)`),digits = 2),
    #                         ",",round(exp(refit.risk$coef+1.96*refit.risk$`se(coef)`),digits = 2),")")
    # refit.risk[,(ncol(refit.risk)-1)] <- formatC(refit.risk[,(ncol(refit.risk)-1)], format = "e", digits = 3) #format(round(refit.risk[,ncol(refit.risk)],digits = 5),scientific = TRUE)
    # 
    # refit.risk <- refit.risk[,c(1,2,6,5)]
    # colnames(refit.risk) <- c("Coefficients","Hazard Ratio","CI","P-value")
    # 
    # ClinRefitTable <- as.data.frame(summary(refit.risk.clin)$coefficients)
    # ClinRefitTable$CI <- paste0("(",round(exp(ClinRefitTable$coef-ClinRefitTable$`se(coef)`),digits = 2),
    #                             ",",round(exp(ClinRefitTable$coef+ClinRefitTable$`se(coef)`),digits = 2),")")
    # ClinRefitTable[,(ncol(ClinRefitTable)-1)] <- formatC(ClinRefitTable[,(ncol(ClinRefitTable)-1)], format = "e", digits = 3) #format(ClinRefitTable[,ncol(ClinRefitTable)],scientific = TRUE)
    # ClinRefitTable <- ClinRefitTable[,c(1,2,6,5)]
    # colnames(ClinRefitTable) <- c("Coefficients","Hazard Ratio","CI","P-value")
  }
  
  
  ### ENET ANALYSIS ###
  
  if("ENET" %in% method) {
    
    ConcordanceIndex <- as.data.frame(as.vector(sapply(ENET, "[[", "CI")))
    summary.CI <- round(as.data.frame(c(quantile(ConcordanceIndex[,1],c(0,0.25,0.5,0.75,1)),mean(ConcordanceIndex[,1]))),digits = 2)
    colnames(summary.CI) <- "Concordance Index"
    rownames(summary.CI) <- c("Minimum","1st Quarter","Median","3rd Quarter","Maximum","Mean")
    CI.BP <- as.data.frame(t(summary.CI))
    ### FOR APP
    # CI.BP <- ggplot(ConcordanceIndex, aes(x="",y = ConcordanceIndex)) +
    #   geom_boxplot(outlier.colour = "red", outlier.shape = 1,colour = "#3366FF") +
    #   geom_jitter(width = 0.2) +
    #   theme_grey() +
    #   labs(title = "Concordance Index distribution", subtitle = "Using genetic data only")
    ###
    
    if(LT){genes <- colnames(data)[-c(1:3)]}
    if(!LT){genes <- colnames(data)[-c(1:2)]}
    selected.genes.enet <- as.data.frame(matrix(0L,nrow=1,ncol = length(genes)))
    colnames(selected.genes.enet) <- genes
    for(i in 1:length(ENET)){
      temp <-ENET[[i]]$fit
      if(!is.null(temp)){
        selected.temp <- rownames(temp)
        selected.genes.enet[,match(selected.temp,colnames(selected.genes.enet))] <- 
          selected.genes.enet[,match(selected.temp,colnames(selected.genes.enet))] +1
      }
    }
    
    topHits <- names(sort(selected.genes.enet,decreasing = TRUE))[1:10]
    if(length(selected.genes.enet)  >= 20){ranked.enet <- sort(selected.genes.enet,decreasing = TRUE)[1:20]}
    if(length(selected.genes.enet)  < 20){ranked.enet <- sort(selected.genes.enet,decreasing = TRUE)[1:length(selected.genes.enet)]}
    melt.rank.enet <- melt(ranked.enet)
    colnames(melt.rank.enet) <- c("Gene","Count")
    
    influencePlot <- ggplot(melt.rank.enet,aes(x=Gene,y=Count,fill=Gene))+geom_col()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      #scale_fill_discrete(guide = guide_legend(title = "Category")) +
      #scale_fill_manual(values=rainbow(length(unique(melt.rank.lasso$Gene)))) +
      labs(title = "Most recurrent selected genes out of 200 runs") +
      theme(legend.position="none")
    
    
    if(LT) Variables <- colnames(data)[-c(1:3)]
    if(!LT) Variables <- colnames(data)[-c(1:2)]
    allCoefs <- as.data.frame(matrix(nrow=length(ENET),ncol=length(Variables)))
    colnames(allCoefs) <- Variables
    
    for(x in 1:length(ENET)){
      coefsValues <- ENET[[x]]$fit[,1]
      allCoefs[x,match(names(coefsValues),colnames(allCoefs))] <- as.numeric(coefsValues)
    }
    allCoefs[is.na(allCoefs)] <- 0
    
    meanCoefs <- apply(allCoefs,2,function(x){mean(x,na.rm = TRUE)})
    selectFreq <- apply(allCoefs,2,function(x){
      length(which(x!=0))/length(x)
    })
    
    ## get mu freq
    if(LT) data.temp <- data[,-c(1:3)]
    if(!LT) data.temp <- data[,-c(1:2)]
    MutationFrequency <- apply(data.temp,2,function(x){
      sum(x)/length(x)
    })
    
    resultsAll <- as.data.frame(cbind(meanCoefs,selectFreq,MutationFrequency))
    colnames(resultsAll) <- c("MeanCoefficient","SelectionFrequency","MutationFrequency")
    rownames(resultsAll) <- names(meanCoefs)
    resultsAll <- resultsAll[complete.cases(resultsAll),]
    resultsAll$GeneName <- rownames(resultsAll)
    resultsAll$MutationFrequency2 <- cut(resultsAll$MutationFrequency, c(0,0.10,0.20,0.40))
    selectInflPlot <- plot_ly(data = resultsAll, x = ~MeanCoefficient, y = ~SelectionFrequency,
                              text = paste('Gene :'~GeneName,
                                           'Hazard Ratio :' ~round(exp(MeanCoefficient)),digits=2), mode = "markers",size = ~MutationFrequency,color = ~MutationFrequency) %>% 
      layout(title ="Volcano Plot")
    
    # selectInflPlot <- ggplot(resultsAll, aes(x=MeanCoefficient, y=SelectionFrequency)) + 
    #   geom_point(aes(size=MutationFrequency,color=MutationFrequency)) +
    #   geom_vline(xintercept=0) + 
    #   geom_text(data=subset(resultsAll,SelectionFrequency > 0.4),label=rownames(subset(resultsAll,SelectionFrequency > 0.4))) +
    #   theme(text = element_text(size=12)) +
    #   scale_x_continuous(limits = c(-max(abs(resultsAll$MeanCoefficient)),max(abs(resultsAll$MeanCoefficient)))) + 
    #   scale_colour_gradient(low = "cyan3",high = "blue")
    
    final.pred <- as.data.frame(matrix(nrow= nrow(data),ncol = length(ENET)))
    rownames(final.pred) <- rownames(data)
    for(i in 1:length(ENET)){
      temp <- ENET[[i]]$predicted
      final.pred[match(names(temp),rownames(final.pred)),i] <- as.numeric(temp)
      
    }
    average.risk <- apply(final.pred,1,function(x){
      mean(as.numeric(x),na.rm = TRUE)
    })
    #TopHits<- ranked.enet[1:10]
    if(LT) refit.risk <- coxph(Surv(data$time1,data$time2,data$status)~average.risk)
    if(!LT) refit.risk <- coxph(Surv(data$time,data$status)~average.risk)
    
    RiskHistogram <- ggplot(Risk, aes(x = average.risk, y = ..density..)) +
      geom_histogram(show.legend = FALSE, aes(fill=..x..),
                     breaks=seq(min(Risk$average.risk), max(Risk$average.risk), by=0.05)) +
      geom_density(show.legend = FALSE) +
      theme_minimal() +
      labs(x = "Average risk score", y = "Proportion") +
      scale_fill_gradient(high = "red", low = "green")
    
    
    
    # RiskHistogram <- ggplot(Risk, aes(Risk$average.risk,..density..)) +
    #   geom_histogram(show.legend = FALSE,aes(fill=..x..),
    #                  breaks=seq(min(Risk$average.risk), max(Risk$average.risk), by=0.05)) +
    #   geom_density(show.legend = FALSE) +
    #   theme_minimal() +
    #   labs(x = "Average risk score", y = "Proportion") +
    #   scale_fill_gradient("Count", low="green", high="red")
  }
  
  
  ### CF ANALYSIS ###
  
  if("CF" %in% method){
    
    ConcordanceIndex <- as.data.frame(as.vector(sapply(CF, "[[", "CI")))
    summary.CI <- round(as.data.frame(c(quantile(ConcordanceIndex[,1],c(0,0.25,0.5,0.75,1)),mean(ConcordanceIndex[,1]))),digits = 2)
    colnames(summary.CI) <- "Concordance Index"
    rownames(summary.CI) <- c("Minimum","1st Quarter","Median","3rd Quarter","Maximum","Mean")
    CI.BP <- as.data.frame(t(summary.CI))
    CI.BP <- ggplot(ConcordanceIndex, aes(x="",y = ConcordanceIndex)) +
      geom_boxplot(outlier.colour = "red", outlier.shape = 1,colour = "#3366FF") +
      geom_jitter(width = 0.2) +
      theme_grey() +
      labs(title = "Concordance Index distribution", subtitle = "Using genetic data only")
    
    if(LT){genes <- colnames(data)[-c(1:3)]}
    if(!LT){genes <- colnames(data)[-c(1:2)]}
    
    influence <- as.data.frame(matrix(ncol=length(genes),nrow = length(CF) ))
    colnames(influence) <- genes
    
    for ( i in 1:length(CF)){
      temp <- CF[[i]]$Vars
      for(j in 1 : length(temp)){
        influence[i,match(names(temp)[j],colnames(influence))] <- temp[j]
      }
    }
    
    
    average.influence  <- as.data.frame(matrix(nrow= 2, ncol = length(genes)))
    colnames(average.influence) <- genes
    for(i in 1:ncol(influence)){
      average.influence[1,i] <- mean(influence[,i],na.rm = TRUE)
      average.influence[2,i] <- median(influence[,i],na.rm = TRUE)
    }
    
    topHits <- names(sort(average.influence[2,],decreasing = TRUE))[1:10]
    
    topHitsMat <- influence[,match(topHits,colnames(influence))]
    topHitsMat.melt <- melt(topHitsMat)
    colnames(topHitsMat.melt) <- c("Gene","Influence")
    
    influencePlot <- ggplot(topHitsMat.melt, aes(x=Gene,y = Influence)) +
      geom_boxplot(outlier.colour = "red", outlier.shape = 1,colour = "#3366FF") +
      #geom_jitter(width = 0.05) +
      theme_grey() +
      labs(title = "Relative Variable Importance", subtitle = "Using genetic data only")
    
    all.selections <- vector(length = length(CF[[1]]$Selection),mode = "numeric")
    numTrees <- 0 
    for(i in 1:length(CF)){
      all.selections <- all.selections + CF[[i]]$Selection
      numTrees <- numTrees + CF[[i]]$count
    }
    all.selections <- all.selections / (numTrees)
    
    if(LT) data.temp <- data[,-c(1:3)]
    if(!LT) data.temp <- data[,-c(1:2)]
    MutationFrequency <- apply(data.temp,2,function(x){
      sum(x)/length(x)
    })
    
    SelectInfl <- as.data.frame(t(rbind(average.influence[2,],all.selections,MutationFrequency)))
    colnames(SelectInfl) <- c("Influence","SelectionFrequency","MutationFrequency")
    SelectInfl$GeneName <- rownames(SelectInfl)
    SelectInfl$MutationFrequency2 <- cut(SelectInfl$MutationFrequency, c(0,0.10,0.20,0.40))
    selectInflPlot <- plot_ly(data = SelectInfl, x = ~Influence, y = ~SelectionFrequency,
                              text = ~GeneName, mode = "markers",size = ~MutationFrequency,color = ~MutationFrequency) %>% 
      layout(title ="Volcano Plot")
    
    final.pred <- as.data.frame(matrix(nrow= nrow(data),ncol = length(CF)))
    rownames(final.pred) <- rownames(data)
    for(i in 1:length(CF)){
      temp <- CF[[i]]$predicted
      final.pred[match(temp[,1],rownames(final.pred)),i] <- as.numeric(temp[,2])
      
    }
    
    average.risk <- apply(final.pred,1,mean,na.rm=TRUE)
    if(LT) refit.risk <- coxph(Surv(data$time1,data$time2,data$status)~average.risk)
    if(!LT) refit.risk <- coxph(Surv(data$time,data$status)~average.risk)
    
    RiskHistogram <- ggplot(Risk, aes(Risk$average.risk,..density..)) +
      geom_histogram(show.legend = FALSE,aes(fill=..density..),
                     breaks=seq(min(Risk$average.risk), max(Risk$average.risk), by=0.05)) +
      geom_density(show.legend = FALSE) +
      theme_minimal() +
      labs(x = "Average risk score", y = "Proportion") +
      scale_fill_gradient("Count", low="green", high="red")
  }
  
  ### RANDOMFOREST SRC ###
  if("RF" %in% method){
    
    ConcordanceIndex <- as.data.frame(as.vector(sapply(RF, "[[", "CI")))
    summary.CI <- round(as.data.frame(c(quantile(ConcordanceIndex[,1],c(0,0.25,0.5,0.75,1)),mean(ConcordanceIndex[,1]))),digits = 2)
    colnames(summary.CI) <- "Concordance Index"
    rownames(summary.CI) <- c("Minimum","1st Quarter","Median","3rd Quarter","Maximum","Mean")
    CI.BP <- as.data.frame(t(summary.CI))
    CI.BP <- ggplot(ConcordanceIndex, aes(x="",y = ConcordanceIndex)) +
      geom_boxplot(outlier.colour = "red", outlier.shape = 1,colour = "#3366FF") +
      geom_jitter(width = 0.2) +
      theme_grey() +
      labs(title = "Concordance Index distribution", subtitle = "Using genetic data only")
    
    if(LT){genes <- colnames(data)[-c(1:3)]}
    if(!LT){genes <- colnames(data)[-c(1:2)]}
    
    ranking <- as.data.frame(matrix(ncol=length(genes),nrow = length(RF) ))
    colnames(ranking) <- genes
    
    for ( i in 1:length(RF)){
      temp <- RF[[i]]$Vars
      for(j in 1 : nrow(temp)){
        ranking[i,match(rownames(temp)[j],colnames(ranking))] <- j
      }
    }
    
    average.ranking  <- as.data.frame(matrix(nrow= 2, ncol = length(genes)))
    colnames(average.ranking) <- genes
    for(i in 1:ncol(ranking)){
      average.ranking[1,i] <- mean(ranking[,i],na.rm = TRUE)
      average.ranking[2,i] <- median(ranking[,i],na.rm = TRUE)
    }
    ranked.averages <- sort(average.ranking[1,],decreasing = FALSE)
    ranked.medians <- sort(average.ranking[2,],decreasing = FALSE)
    
    # for 10 top medians make boxplot
    topHits <- names(ranked.medians[1,1:10])
    topHitsMat <- ranking[,match(topHits,colnames(ranking))]
    topHitsMat.melt <- melt(topHitsMat)
    colnames(topHitsMat.melt) <- c("Gene","Rank")
    
    influencePlot <- ggplot(topHitsMat.melt, aes(x=Gene,y = Rank)) +
      geom_boxplot(outlier.colour = "red", outlier.shape = 1,colour = "#3366FF") +
      #geom_jitter(width = 0.05) +
      theme_grey() +
      labs(title = "Concordance Index distribution", subtitle = "Using genetic data only")
    
    final.pred <- as.data.frame(matrix(nrow= nrow(data),ncol = length(RF)))
    rownames(final.pred) <- rownames(data)
    for(i in 1:length(RF)){
      temp <- RF[[i]]$predicted
      final.pred[match(temp[,1],rownames(final.pred)),i] <- as.numeric(temp[,2])
    }
    average.risk <- c()
    for(j in 1:nrow(final.pred)){
      average.risk[j] <- mean(as.numeric(final.pred[j,]),na.rm = TRUE)
    }
    if(LT) refit.risk <- coxph(Surv(data$time1,data$time2,data$status)~average.risk)
    if(!LT) refit.risk <- coxph(Surv(data$time,data$status)~average.risk)
    
    RiskHistogram <- ggplot(Risk, aes(Risk$average.risk,..density..)) +
      geom_histogram(show.legend = FALSE,aes(fill=..density..),
                     breaks=seq(min(Risk$average.risk), max(Risk$average.risk), by=0.05)) +
      geom_density(show.legend = FALSE) +
      theme_minimal() +
      labs(x = "Average risk score", y = "Proportion") +
      scale_fill_gradient("Count", low="green", high="red")
  }
  
  # return(list("ciSummary" = CI.BP,"inflPlot" = influencePlot,"topHits" = topHits,"average.risk"=average.risk,"data.out"= data,
  #             "selectInflPlot" = selectInflPlot,"MethodUsed" = method,"RiskRefit"=refit.risk,"ClinRefitTable"=ClinRefitTable,
  #             "RiskHistogram"=RiskHistogram,"LassoFits"=allCoefs,"ClinRefit"=refit.risk.clin,"time.type"=time.type,"MD"=MD,
  #             "RiskScoreSummary"=as.data.frame(t(summary.RiskScore))))
 return(list("selectInflPlot" = selectInflPlot))
}
#test <- getResults(studyType="Lung",method="LASSO",CNV=FALSE,OnlyCNV= FALSE)


#FirstRun <- getResults(studyType="Lung",method="LASSO",CNV=FALSE,OnlyCNV= FALSE,geneList = NULL)
#save(FirstRun,file="./FirstRun.Rdata")
