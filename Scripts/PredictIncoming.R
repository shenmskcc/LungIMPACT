############################################
### PREDICTING INCOMING PATIENT SURVIVAL ###
############################################

### Method

## INPUT:
# 1) Set of mutated genes : All genes not inputed in string will be considered as wild type (from App).
# 2) LASSO runs : All refited (forced beta parameters) cox models
# 3) Set of covariates : All covariates available can be edited and added to the model
# 4) The character vector list of genes available (return colnames of data in GetResults function)
# 5) Refitted model with clinical covariates

## METHOD :
# Iterate through all the runs, using the coefficients found in the LASSO models generate a predicted 
# risk score for the incoming patient. Use the predicted risk score in a refitted coxph model with the 
# given clinical covariates.

## OUTPUT :
# 1) Survival curve for the patients inputed
# 2) Summary Table

### TEST

predictIncomingPatient <- function(mutGenes,clinical,ClinRefit,time.type,MD,LassoFits,RiskScore){
  
  #if(length(mutGenes) == 0){stop("INPUT ERROR : Enter at least one gene name")}
  #library(plotly)
  ### Step 1 : Build the dataset ###
  ## Gene mutation data
  # create frame 
  mut <- as.data.frame(matrix(0L,nrow=1,ncol=ncol(LassoFits)))
  colnames(mut) <- colnames(LassoFits)
  rownames(mut) <- c("Patient")
  
  if(length(mutGenes) > 0){
    mutGenes <- gsub(" ","",mutGenes)}
  
  # Check available genes
  if(anyNA(match(mutGenes,colnames(mut)))){stop(paste("The mutated gene you have inputed :",
                                                      mutGenes[which(is.na((match(mutGenes,colnames(mut)))))],
                                                      "is not available at the moment") )}
  mut[,match(mutGenes,colnames(mut))] <- 1
  
  ## Mutation data complete
  
  ## Clinical 
  ### NEED A WHOLE PART FROM THE APP HERE TRANSFORMING TO VALUES IN THE DATA ###
  
  ## Create frame
  avClinNames <- names(ClinRefit$coefficients)[-length(names(ClinRefit$coefficients))]
  clin <- as.data.frame(matrix(0L,nrow=1,ncol=length(avClinNames)))
  colnames(clin) <- avClinNames 
  rownames(clin) <- "Patient"
  
  # match input from the app
  clin[,match(clinical,colnames(clin))] <- 1
  
  ### Part 2 : Get the predicted risk score
  LassoFits <- as.matrix(LassoFits)
  mut <- as.vector(mut)
  RiskScore.new <- mean(LassoFits%*%t(mut)) 
  RiskScoreRange <- range(RiskScore)
  RiskScore <- rescale(RiskScore.new, to = c(0, 10), from = RiskScoreRange)
  RiskScore <- RiskScore[length(RiskScore)]
  ### part 3 : Refit with the given refitted clinical variable
  clin <- as.data.frame(cbind(clin,RiskScore))
  
  ### FIRST RUN
  
  
  ### GET CI ###
  survival.probs <- as.data.frame(matrix(nrow=(nrow(clin)+3),ncol=15))
  rownames(survival.probs) <- c(rownames(clin),"Lower","Upper","Time")
  surv.temp <- survfit(ClinRefit, newdata = clin)
  for(i in 1:ncol(survival.probs)){
    survival.probs[,i] <- c(as.numeric(summary(surv.temp, times = (i*3-3))$surv),
                            round(summary(surv.temp, times = (i*3-3))$lower,digits=2),
                            round(summary(surv.temp, times = (i*3-3))$upper,digits=2),
                            i*3-3)
  }
  
  a <- list(
    autotick = FALSE,
    dtick = 6,
    tickcolor = toRGB("black")
  )
  
  t.survival.probs <- as.data.frame(t(survival.probs))
  y <- list(
    title = "Survival Probability"
  )
  IndSurvKM <- plot_ly(t.survival.probs, x = ~Time, y = ~Patient, name = 'Estimated Survival', type = 'scatter',
                       mode = 'lines+markers',hoverinfo="hovertext",#color = ~Patient
                       hovertext = ~paste("Genetic Risk Score :",round(RiskScore,digits=3))
  ) %>% layout(yaxis = y,xaxis = ~a) %>%
    layout(shapes=list(type='line', x0= 3*MD, x1= 3*MD, y0=0,
                       y1=1, line=list(dash='dot', width=1,color = "red")),
           xaxis = list(title = paste0("Time (",time.type,")"), showgrid = TRUE)) %>%
    add_ribbons(data = t.survival.probs,
                ymin = ~Lower,
                ymax = ~Upper,
                line = list(color = 'rgba(7, 164, 181, 0.05)'),
                fillcolor = 'rgba(7, 164, 181, 0.2)',
                name = "Confidence Interval")
  
  
  ### MAKE TABLE SUMMARY
  survivalSummary <- as.data.frame(matrix(nrow=1,ncol=5))
  rownames(survivalSummary) <- "New patient"
  colnames(survivalSummary) <- c("MedianOS","95%CI","1Ysurvival","3Ysurvival","Predicted genetic risk")
  # for each group find closest value to median 
  if(time.type == "Months"){YR1 <- 1*12;YR3 <- 3*12}
  if(time.type == "Days"){YR1 <- 1*365;YR3 <- 3*365}
  Fit <- surv.temp 
  #Fit <- survfit(ClinRefit, newdata = clin)#, se.fit = T, conf.int = "log-log")
  med.index <- which.min(abs(Fit$surv-0.5))
  YR1.index <- which.min(abs(Fit$time-YR1))
  YR3.index <- which.min(abs(Fit$time-YR3))
  survivalSummary[1,] <- c(round(Fit$time[med.index],digits=2),paste0("(",round(Fit$time[which.min(abs(Fit$lower-0.5))],digits=2),",",
                                                                      round(Fit$time[which.min(abs(Fit$upper-0.5))],digits=2),")"),
                           round(Fit$surv[YR1.index],digits=2),round(Fit$surv[YR3.index],digits=2),
                           round(RiskScore,digits=3))
  return(list("IndSurvKM"=IndSurvKM,"IndPredTable"=survivalSummary,"RiskScore"=RiskScore))
  #}
}



