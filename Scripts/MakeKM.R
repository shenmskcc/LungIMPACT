MakeKM <- function(data,average.risk,topHits,LT,numGroups,cuts,geneList){
  
  data$lvl4Groups <- rep(NA,nrow(data))
  qts <- as.numeric(quantile(average.risk,cuts))
  
  if(numGroups ==2){
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
  
  
  count.dups <- function(DF){

    DT <- data.table(DF)
    DT[,.N, by = names(DT)]
  }

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
  
  return(list("PieChart"=pie.chart,"GenesUsed"=paste0(useGenes,collapse = ",")))
}

# RiskGroupsResults <- MakeKM(data,average.risk,topHits,LT,numGroups,cuts,geneList = NULL)
# 
# # Make2KM(data = data.gen,average.risk = averaged.risk.gen.boot.bag,
# #         file = "GenOnly_Bagging",type.data = "GenOnly")
# save(RiskGroupsResults,file="RiskGroupsResults.Rdata")
