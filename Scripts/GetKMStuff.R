KMStuff <-function(data,average.risk,topHits,numGroups,cuts,geneList) {
  
  LT = T
  cuts <- as.numeric(cuts)
  numGroups <- as.numeric(numGroups)
  OUT <- MakeKM(data,average.risk,topHits,LT,numGroups,cuts,geneList)
  #KM <- OUT$KM_Plot
  #Mut <- OUT$mut_Plot
  #survSum <- OUT$SurvSum
  #MajorCasesKM <- OUT$MajorCasesKM
  PieChart <- OUT$PieChart
  GenesUsed <- OUT$GenesUsed
  #coMutPlot <- OUT$coMutation
  
  
  return(list(#"KMPlot" = KM,"MutPlot" =Mut,"SurvSumTable"=survSum,
              "PieChart" = PieChart,"GenesUsed"=GenesUsed))#,"MajorCasesKM" =MajorCasesKM,"coMutPlot" = coMutPlot))
  
}

#test2 <- KMStuff(test$data.out,test$average.risk,test$topHits,2,0.5)
