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

getResults <- function(studyType,method,geneList){

  
  load("./Study/Lung/results/Lung_LASSO.Rdata")
  data <- read.csv("./Study/Lung/data/LungReadyStudy.csv",header = TRUE, row.names = 1)
  
  ## determine if left truncated
  LT = T
  MD = 12
  time.type = "Months"
  
  ### LASSO ANALYSIS ###
  
  #try(setwd("./results/"))
  
  if("LASSO" %in% method) {
    
    Variables <- colnames(data)
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
    data.temp <- data
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
    
  }
  
  # return(list("ciSummary" = CI.BP,"inflPlot" = influencePlot,"topHits" = topHits,"average.risk"=average.risk,"data.out"= data,
  #             "selectInflPlot" = selectInflPlot,"MethodUsed" = method,"RiskRefit"=refit.risk,"ClinRefitTable"=ClinRefitTable,
  #             "RiskHistogram"=RiskHistogram,"LassoFits"=allCoefs,"ClinRefit"=refit.risk.clin,"time.type"=time.type,"MD"=MD,
  #             "RiskScoreSummary"=as.data.frame(t(summary.RiskScore))))
 return(list("selectInflPlot" = selectInflPlot))
}
#test <- getResults(studyType="Lung",method="LASSO",CNV=FALSE,OnlyCNV= FALSE)


# FirstRun <- getResults(studyType="Lung",method="LASSO",CNV=FALSE,OnlyCNV= FALSE,geneList = NULL)
# save(FirstRun,file="FirstRun.Rdata")
