# Make first pass
# load("FirstRun.Rdata")
# data <- FirstRun$data.out
# average.risk <- FirstRun$average.risk
# topHits <- FirstRun$topHits
# LT <- TRUE
# numGroups <- 4
# cuts <- c(0.25,0.75,0.9)

MakeKM <- function(data,average.risk,topHits,LT,numGroups,cuts,geneList){
  
  data$lvl4Groups <- rep(NA,nrow(data))
  qts <- as.numeric(quantile(average.risk,cuts))
  
  if(numGroups ==2){
    #qts <- as.numeric(summary(average.risk)[c(2,3,5)])
    #qts <- as.numeric(quantile(average.risk,cuts))
    for(i in 1:nrow(data)){
      
      # for 2 level groups
      if(average.risk[i] < qts) {
        data$lvl4Groups[i] <- "Low"
      }
      
      if(average.risk[i] >= qts) {
        data$lvl4Groups[i] <- "High"
      }
    }
  }
  
  if(numGroups == 3){
    
    for(i in 1:nrow(data)){
      
      # for 3 level groups
      if(average.risk[i] < qts[1]) {
        data$lvl4Groups[i] <- "Low"
      }
      if(average.risk[i] >= qts[1] && average.risk[i] < qts[2]){
        data$lvl4Groups[i] <- "Intermediate"
      }
      if(average.risk[i] >= qts[2]) {
        data$lvl4Groups[i] <- "High"
      }
    }
  }
  
  if(numGroups == 4){
    #qts <- quantile(average.risk,probs = c(0.25,0.75,0.9))
    for(i in 1:nrow(data)){
      
      # for 4 level groups
      if(average.risk[i] < qts[1]) {
        data$lvl4Groups[i] <- "Low"
      }
      if(average.risk[i] >= qts[1] && average.risk[i] < qts[2]){
        data$lvl4Groups[i] <- "Low-intermediate"
      }
      if(average.risk[i] >= qts[2] && average.risk[i] < qts[3]){
        data$lvl4Groups[i] <- "High-intermediate"
      }
      if(average.risk[i] >= qts[3]) {
        data$lvl4Groups[i] <- "High"
      }
    }
  }
  
  
  ##################
  
  ### Is this left truncated ?
  # Will be LT if second column is binary 
  if(LT == TRUE) {
    colnames(data)[1:3] <- c("time1","time2","status")
    survObj <<- with(data, Surv(time=time1, time2=time2, event=status))
    if(max(data$time2) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }
  if(LT == FALSE) {
    colnames(data)[1:2] <- c("time","status")
    survObj <<- with(data, Surv(time=time, event=status))
    if(max(data$time) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }
  ## for 4 levels
  # make KM for Lasso 1 :
  if(numGroups == 2){data$lvl4Groups <- factor(data$lvl4Groups, levels =c("Low","High") )}
  if(numGroups == 3){data$lvl4Groups <- factor(data$lvl4Groups, levels =c("Low","Intermediate","High") )}
  if(numGroups == 4){data$lvl4Groups <- factor(data$lvl4Groups, levels =c("Low","Low-intermediate","High-intermediate","High") )}
  data$lvl4Groups <- as.factor(data$lvl4Groups)
  RiskGroup <<- as.factor(data$lvl4Groups)
  #km.lvl4 <<- survfit(survObj ~ RiskGroup,data=data, conf.type = "log-log")
  ### Get pvalue 
  fit0 <- coxph(survObj ~ RiskGroup,data=data, 
                na.action=na.exclude) 
  log.test.pval <- as.vector(summary(fit0)[10][[1]])[3]
  CI <- as.numeric(as.vector(summary(fit0)[14])[[1]][1])
  limit <- as.numeric(quantile(data$time2,0.95))
  KM_2LVLS <- ggsurvplot(survfit(survObj ~ RiskGroup,data=data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = "hv",
                         data = data,xlim=c(0,limit),break.time.by = 6) + xlab("Time (Months)") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =4)," and CI : ",round(CI,digits=4), ")",sep=""))+ 
    geom_vline(xintercept=intercept,col="red", lty = 2)
  
  # png(paste(file,"_KM_2LVLS.png",sep=""),height = 10,width = 14,unit="in",res = 300)
  # print(KM_2LVLS)
  # dev.off()
  
  
  ### MAKE TABLE 
  if(numGroups == 2) {Groups <- c("Low","High")}
  if(numGroups == 3) {Groups <- c("Low","Intermediate","High")}
  if(numGroups == 4) {Groups <- c("Low","Low-intermediate","High-intermediate","High")}
  # survivalGroup <- as.data.frame(matrix(nrow=length(Groups),ncol=4))
  # rownames(survivalGroup) <- Groups
  # colnames(survivalGroup) <- c("MedianOS","95%CI","1Ysurvival","3Ysurvival")
  # # for each group find closest value to median 
  # if(timeType == "Months"){YR1 <- 1*12;YR3 <- 3*12}
  # if(timeType == "Days"){YR1 <- 1*365;YR3 <- 3*365}
  # for(i in 1:length(Groups)){
  #   if(LT == TRUE){NewObject <- with(data[data$lvl4Groups == Groups[i],],Surv(time1,time2,status))}
  #   if(LT == FALSE){NewObject <- with(data[data$lvl4Groups == Groups[i],],Surv(time,status))}
  #   Fit <- survfit(NewObject ~ 1,data=data[data$lvl4Groups == Groups[i],], conf.type = "log-log")
  #   med.index <- which.min(abs(Fit$surv-0.5))
  #   YR3.index <- which.min(abs(Fit$time-YR1))
  #   YR5.index <- which.min(abs(Fit$time-YR3))
  #   survivalGroup[i,] <- c(round(Fit$time[med.index],digits=2),paste0("(",round(Fit$time[which.min(abs(Fit$lower-0.5))],digits=2),",",
  #                                                                     round(Fit$time[which.min(abs(Fit$upper-0.5))],digits=2),")"),
  #                          round(Fit$surv[YR3.index],digits=2),round(Fit$surv[YR5.index],digits=2))
  # }
  
  
  ### MAKE CORRESPONDING MUTATION PER GROUP PLOT ###
  
  # mutDistrib <- as.data.frame(matrix(nrow = numGroups,ncol = length(topHits)))
  # rownames(mutDistrib) <- Groups
  # colnames(mutDistrib) <- topHits
  # for( gene in 1:length(topHits)){
  #   if(length(data[data$lvl4Groups == "Low",match(topHits[gene],colnames(data))]) != 0 ){
  #     mutDistrib[match("Low",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "Low",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "Low",match(topHits[gene],colnames(data))])
  #   }
  #   if(length(data[data$lvl4Groups == "Low-intermediate",match(topHits[gene],colnames(data))]) != 0 ){
  #     mutDistrib[match("Low-intermediate",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "Low-intermediate",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "Low-intermediate",match(topHits[gene],colnames(data))])
  #   }
  #   if(length(data[data$lvl4Groups == "Intermediate",match(topHits[gene],colnames(data))]) != 0 ){
  #     mutDistrib[match("Intermediate",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "Intermediate",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "Intermediate",match(topHits[gene],colnames(data))])
  #   }
  #   if(length(data[data$lvl4Groups == "High-intermediate",match(topHits[gene],colnames(data))]) != 0  ){
  #     mutDistrib[match("High-intermediate",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "High-intermediate",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "High-intermediate",match(topHits[gene],colnames(data))])
  #   }
  #   if(length(data[data$lvl4Groups == "High",match(topHits[gene],colnames(data))]) != 0 ){
  #     mutDistrib[match("High",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "High",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "High",match(topHits[gene],colnames(data))])
  #   }
  # }
  # mutDistrib$Risk <- Groups
  # melted.prop <- melt(mutDistrib)
  # colnames(melted.prop) <- c("Risk","Gene","Proportion")
  # #melted.prop$Risk <- as.factor(melted.prop$Risk)#c("Low","High")
  # melted.prop$Risk <- factor(melted.prop$Risk, levels =Groups )
  # mut_2LVLS <- ggplot(melted.prop, aes(x = Gene, y = Proportion, fill = Risk)) +
  #   geom_bar(stat = "identity", position=position_dodge()) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   labs(title = "Proportion of mutations for each genes per risk and method", subtitle = "Using genetic data only")
  # 
  ### PIE CHARTS ###
  count.dups <- function(DF){

    DT <- data.table(DF)
    DT[,.N, by = names(DT)]
  }

  ### FIT topHits ###
  fit.data <- data[,match(c("time1","time2","status",topHits[1:5]),colnames(data))]
  fit.topHits <- coxph(Surv(time1,time2,status)~.,data= fit.data)

  if( LT && max(data$time2) > 1000){MD = 365
  time.type = "Days"}
  if( LT && max(data$time2) < 1000){MD = 12
  time.type = "Months"}
  if( !LT && max(data$time) > 1000){MD = 365
  time.type = "Days"}
  if( !LT && max(data$time) < 1000){MD = 12
  time.type = "Months"}

  ####
  if(length(geneList) == 0) {useGenes <- topHits[1:5]}
  else{
    geneList <- gsub(" ","",geneList)
    useGenes <- geneList}
  pie.data <- data[,match(useGenes,colnames(data))]
  Groups <- as.character(data[,match("lvl4Groups",colnames(data))])
  profile <- apply(pie.data,1,function(x){
    return(paste0(paste0(names(x),"=",as.numeric(x),collapse =",")))
  })
  make.pie.data <- as.data.frame(cbind(Groups,profile))
  
  ## Save top 2 scenarios for KM
  top.scenarios <- as.data.frame(matrix(nrow = 8,ncol = 5))
  colnames(top.scenarios) <- useGenes
  start <-1
  
  # high
  high.pie <- try(filter(make.pie.data,Groups == "High"))
  profiles.high <- cbind(count.dups(high.pie)[order(-count.dups(high.pie)$N),],
                         paste("Profile",1:nrow(count.dups(high.pie))))
  colnames(profiles.high) <- c("Groups","profile","N","Profile")
  scenarios <- lapply(1:nrow(profiles.high),function(x){
    geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.high$profile[x]),split=","))))
  })
  profiles.high$profile <- unlist(lapply(1:nrow(profiles.high),function(x){
    temp <- unlist(strsplit(as.character(profiles.high$profile[x]),split = ","))
    prof.temp <- temp[grep("=1",temp)]
    prof <- toString(gsub("=1","",prof.temp))
  }) )
  scenarios <- as.data.frame(do.call(rbind, scenarios))
  colnames(scenarios) <- useGenes
  #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
  #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
  #profiles.high <- profiles.high #cbind(profiles.high,Prob3YSurvival,Prob5YSurvival)
  profiles.high$profile[which(profiles.high$profile =="")] <- "No mutants"
  top.scenarios[start:(start+2),] <- scenarios[1:3,]
  rownames(top.scenarios)[start:(start+2)] <- c("High : Profile 1","High : Profile 2","High : Profile 3")
  start <- start+3
  ## Inter high
  if(numGroups==4){
    high.inter.pie <- try(filter(make.pie.data,Groups == "High-intermediate"))
    profiles.inter.high <- cbind(count.dups(high.inter.pie)[order(-count.dups(high.inter.pie)$N),],
                                 paste("Profile",1:nrow(count.dups(high.inter.pie))))
    colnames(profiles.inter.high) <- c("Groups","profile","N","Profile")
    scenarios <- lapply(1:nrow(profiles.inter.high),function(x){
      geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.inter.high$profile[x]),split=","))))
    })
    profiles.inter.high$profile <- unlist(lapply(1:nrow(profiles.inter.high),function(x){
      temp <- unlist(strsplit(as.character(profiles.inter.high$profile[x]),split = ","))
      prof.temp <- temp[grep("=1",temp)]
      prof <- toString(gsub("=1","",prof.temp))
    }) )
    scenarios <- as.data.frame(do.call(rbind, scenarios))
    colnames(scenarios) <- useGenes
    #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
    #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
    #profiles.inter.high <- profiles.inter.high #cbind(profiles.inter.high,Prob3YSurvival,Prob5YSurvival)
    profiles.inter.high$profile[which(profiles.inter.high$profile =="")] <- "No mutants"
    top.scenarios[start:(start+2),] <- scenarios[1:3,]
    rownames(top.scenarios)[start:(start+2)] <- c("High-Intermediate : Profile 1","High-Intermediate : Profile 2","High-Intermediate : Profile 3")
    start <- start+3}
  
  ## INTERMEDIATE
  if(numGroups==3){
    inter.pie <- try(filter(make.pie.data,Groups == "Intermediate"))
    profiles.inter <- cbind(count.dups(inter.pie)[order(-count.dups(inter.pie)$N),],
                            paste("Profile",1:nrow(count.dups(inter.pie))))
    colnames(profiles.inter) <- c("Groups","profile","N","Profile")
    scenarios <- lapply(1:nrow(profiles.inter),function(x){
      geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.inter$profile[x]),split=","))))
    })
    profiles.inter$profile <- unlist(lapply(1:nrow(profiles.inter),function(x){
      temp <- unlist(strsplit(as.character(profiles.inter$profile[x]),split = ","))
      prof.temp <- temp[grep("=1",temp)]
      prof <- toString(gsub("=1","",prof.temp))
    }) )
    scenarios <- as.data.frame(do.call(rbind, scenarios))
    colnames(scenarios) <- useGenes
    #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
    #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
    #profiles.inter <- profiles.inter #cbind(profiles.inter,Prob3YSurvival,Prob5YSurvival)
    profiles.inter$profile[which(profiles.inter$profile =="")] <- "No mutants"
    top.scenarios[start:(start+2),] <- scenarios[1:3,]
    rownames(top.scenarios)[start:(start+2)] <- c("Intermediate : Profile 1","Intermediate : Profile 2","Intermediate : Profile 3")
    start <- start+3}
  
  
  ## LOW INTER
  if(numGroups==4){
    low.inter.pie <- try(filter(make.pie.data,Groups == "Low-intermediate"))
    profiles.inter.low <- cbind(count.dups(low.inter.pie)[order(-count.dups(low.inter.pie)$N),],
                                paste("Profile",1:nrow(count.dups(low.inter.pie))))
    colnames(profiles.inter.low) <- c("Groups","profile","N","Profile")
    scenarios <- lapply(1:nrow(profiles.inter.low),function(x){
      geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.inter.low$profile[x]),split=","))))
    })
    profiles.inter.low$profile <- unlist(lapply(1:nrow(profiles.inter.low),function(x){
      temp <- unlist(strsplit(as.character(profiles.inter.low$profile[x]),split = ","))
      prof.temp <- temp[grep("=1",temp)]
      prof <- toString(gsub("=1","",prof.temp))
    }) )
    scenarios <- as.data.frame(do.call(rbind, scenarios))
    colnames(scenarios) <- useGenes
    #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
    #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
    profiles.inter.low <- profiles.inter.low #cbind(profiles.inter.low,Prob3YSurvival,Prob5YSurvival)
    profiles.inter.low$profile[which(profiles.inter.low$profile =="")] <- "No mutants"
    top.scenarios[start:(start+2),] <- scenarios[1:3,]
    rownames(top.scenarios)[start:(start+2)] <- c("Low-Intermediate : Profile 1","Low-Intermediate : Profile 2","Low-Intermediate : Profile 3")
    start <- start+3
    }
  
  ## LOW
  low.pie <- try(filter(make.pie.data,Groups == "Low"))
  profiles.low <- cbind(count.dups(low.pie)[order(-count.dups(low.pie)$N),],
                        paste("Profile",1:nrow(count.dups(low.pie))))
  colnames(profiles.low) <- c("Groups","profile","N","Profile")
  scenarios <- lapply(1:nrow(profiles.low),function(x){
    geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.low$profile[x]),split=","))))
  })
  profiles.low$profile <- unlist(lapply(1:nrow(profiles.low),function(x){
    temp <- unlist(strsplit(as.character(profiles.low$profile[x]),split = ","))
    prof.temp <- temp[grep("=1",temp)]
    prof <- toString(gsub("=1","",prof.temp))
  }) )
  scenarios <- as.data.frame(do.call(rbind, scenarios))
  colnames(scenarios) <- useGenes
  #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
  #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
  #profiles.low <- cbind(profiles.low,Prob3YSurvival,Prob5YSurvival)
  profiles.low$profile[which(profiles.low$profile =="")] <- "No mutants"
  top.scenarios[start:(start+2),] <- scenarios[1:3,]
  rownames(top.scenarios)[start:(start+2)] <- c("Low : Profile 1","Low : Profile 2","Low : Profile 3")
  start <- start+3
  

  
  if(numGroups == 2){
    pie.chart <- plot_ly() %>%
      add_pie(data = profiles.low, labels = ~profile, values = ~N,
              name = "Low Group", domain = list(x = c(0, 0.4), y = c(0.3, 0.8)),
              hoverinfo="hovertext",
              hovertext = ~paste('',Profile
                                 #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                 #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                                 ),
              type="pie") %>%
      add_pie(data = profiles.high, labels = ~profile, values = ~N,
              name = "High Group", domain = list(x = c(0.45, 0.95), y = c(0.3, 0.8)),hoverinfo="hovertext",
              hovertext = ~paste('',Profile
                                 #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                 #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                                 ),
              type="pie") %>%
      layout(title = "Genetic Profile Distribution", showlegend = F,
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  }
  
  if(numGroups == 3){
    pie.chart <- plot_ly() %>%
      add_pie(data = profiles.low, labels = ~profile, values = ~N,
              name = "Low Group", domain = list(x = c(0, 0.4), y = c(0.5, 1)),hoverinfo="hovertext",
              hovertext = ~paste('',Profile
                                 #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                 #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                                 ),
              type="pie") %>%
      add_pie(data = profiles.inter, labels = ~profile, values = ~N,
              name = "Intermediate Group", domain = list(x = c(0.45, 0.95), y = c(0.5, 1)),hoverinfo="hovertext",
              hovertext = ~paste('',Profile
                                 #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                 #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                                 ),
              type="pie") %>%
      add_pie(data = profiles.high, labels = ~profile, values = ~N,
              name = "High Group", domain = list(x = c(0.2, 0.6), y = c(0, 0.5)),hoverinfo="hovertext",
              hovertext = ~paste('',Profile
                                 #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                 #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                                 ),
              type="pie") %>%
      layout(title = "Genetic Profile Distribution", showlegend = F,
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  }
  
  if(numGroups == 4){
    pie.chart <- plot_ly() %>%
      add_pie(data = profiles.low, labels = ~profile, values = ~N,
              name = "Lo", domain = list(x = c(0, 0.22), y = c(0, 0.9)),hoverinfo="hovertext",
              hovertext = ~paste('',Profile
                                 #'</br> Mutant genes',profile
                                 #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                 #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                                 ),
              type="pie") %>%
      add_pie(data = profiles.inter.low, labels = ~profile, values = ~N,
              name = "ILo", domain = list(x = c(0.26, 0.48), y = c(0, 0.9)),hoverinfo="hovertext",
              hovertext = ~paste('',Profile
                                 #'</br> Mutant genes',profile
                                 #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                 #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                                 ),
              type="pie") %>%
      add_pie(data = profiles.inter.high, labels = ~profile, values = ~N,
              name = "IHi", domain = list(x = c(0.52, 0.74), y = c(0, 0.9)),hoverinfo="hovertext",
              hovertext = ~paste('',Profile
                                 #'</br> Mutant genes',profile
                                 #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                 #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                                 ),
              type="pie") %>%
      add_pie(data = profiles.high, labels = ~profile, values = ~N,
              name = "Hi", domain = list(x = c(0.78, 1), y = c(0, 0.9)),hoverinfo="hovertext",
              hovertext = ~paste('',Profile
                                 #'</br> Mutant genes',profile
                                 #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                 #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                                 ),
              type="pie") %>%
      add_annotations(x= 0.05, y= 1, xref = "paper", yref = "paper", text = "<b>Low risk</b>", showarrow = F) %>%
      add_annotations(x= 0.3, y= 1, xref = "paper", yref = "paper", text = "<b>Intermediate low</b>", showarrow = F) %>%
      add_annotations(x= 0.61, y= 1, xref = "paper", yref = "paper", text = "<b>Intermediate high</b>", showarrow = F) %>%
      add_annotations(x= 0.9, y= 1, xref = "paper", yref = "paper", text = "<b>High risk</b>", showarrow = F) %>%
      layout(title = "Genetic Profile Distribution", showlegend = F,
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  }
  
  return(list(#"KM_Plot" = KM_2LVLS , "mut_Plot" = mut_2LVLS,"SurvSum" = survivalGroup,#,"Profile" = major.profile
             "PieChart"=pie.chart,"GenesUsed"=paste0(useGenes,collapse = ","))) #"MajorCasesKM" = MajorCasesKM,,"coMutation" = coMutPlot))
}

 #RiskGroupsResults <- MakeKM(data,average.risk,topHits,LT,numGroups,cuts,geneList = NULL)
# 
# # Make2KM(data = data.gen,average.risk = averaged.risk.gen.boot.bag,
# #         file = "GenOnly_Bagging",type.data = "GenOnly")
# save(RiskGroupsResults,file="RiskGroupsResults.Rdata")
