library(shiny)
library(shinydashboard)
library(plotly)
library(reshape2)
library(survminer)
library(survival)
library(data.table)
library(rpart.plot)
library(dplyr)
library(gplots)
library(scales)

load("LungDemoData.Rdata")

ui <- dashboardPage(
  dashboardHeader(title = "Prognostic models in metastatic lung adenocarcinoma (BETA)",titleWidth = 400),
  
  dashboardSidebar(width = 200,
    sidebarMenu(
    menuItem("Dashboard", tabName = "fixed", icon = icon("bar-chart")),
    menuItem("Risk Group Stratification", tabName = "strat", icon = icon("scissors")),
    menuItem("Gene View", tabName = "gene", icon = icon("gears")),
    menuItem("Patient View", tabName = "patient", icon = icon("user-circle-o"))
  )
  ),
  
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "fixed",
              sidebarLayout(
                sidebarPanel(
                  width=12,
                  h2("Prognostic Performance")
                  # selectInput("StudyName", "Choose a Study :", 
                  #             choices = c("Lung"
                  #             )),
                  # checkboxInput("AddCNV", "Include copy number data", FALSE),
                  # checkboxInput("OnlyCNV", "Consider only copy number data", FALSE),
                  # 
                  # selectInput("Method", "Choose a Method :", 
                  #             choices = c("Penalized regression"="LASSO")),
                  # submitButton("Submit")
                ),
                mainPanel(width=12,
                          #htmlOutput("VariableHeader"),
                          htmlOutput("RiskHeader"),
                          plotOutput("RiskHistogram"),
                          tableOutput("RiskSummary"),
                          htmlOutput("CIHeader"),
                          tableOutput("CI"),
                          htmlOutput("RefitHeader"),
                          tableOutput("RefitRisk"),
                          htmlOutput("CommentPval"),
                          tableOutput("ClinRefit")#,
                          #htmlOutput("EffectHeader"),
                          #htmlOutput("FreqHeader")#,
                          #plotOutput("influencePlot"),
                )
              )
      ),
      # First tab content
      tabItem(tabName = "strat",
              #
              #   sidebarPanel( width = 12,
              #     checkboxInput("ShowVolvano", "Show Volcano Plot", TRUE),
              #                 checkboxInput("ShowPies", "Show Pie Charts", TRUE),
              #                 submitButton("Submit")
              # ),
              sidebarLayout(
                sidebarPanel(width = 12,
                             h2("Risk Group Stratification")
                             ),
                mainPanel(width = 12,
                  htmlOutput("predRiskText"),
                  htmlOutput("KMText"),
                  plotOutput("KM"),
                  tableOutput("SurvSum"),
                  htmlOutput("MutGroupText"),
                  htmlOutput("MutBPText"),
                  plotOutput("Mut")
                )
              )
      ),
      
      tabItem(tabName = "gene",
              h2("Exploratory interactive gene plots"),
              #   sidebarPanel( width = 12,
              #     checkboxInput("ShowVolvano", "Show Volcano Plot", TRUE),
              #                 checkboxInput("ShowPies", "Show Pie Charts", TRUE),
              #                 submitButton("Submit")
              # ),
              sidebarLayout(
                sidebarPanel(width = 12,
                             textInput("GeneListRisk", 
                                       "Find gene(s) : ",
                                       value = ""),
                             submitButton("Submit")),
                mainPanel(
                  htmlOutput("VolcanoHeader"),width = 12,
                  plotlyOutput("effectPlot"),
                  htmlOutput("profiletext"),
                  plotlyOutput("ProfilePie")#, width="800px", height="400px")
                )
              )
      ),
      
      tabItem(tabName = "patient",
              h2("Predicting survival for an incoming patient"),
              sidebarLayout(
                sidebarPanel(width = 12,
                             textInput("GeneList", 
                                       "Names of the genes you wish to use to create predictive risk : example STK11,KEAP1,KRAS",
                                       value = ""),
                             ##### FOR LUNG ONLY #####
                             checkboxGroupInput(inputId = "Demographics", label = "Choose demographic variables :",
                                                choices = c("Male" = "Sex",
                                                            "Smoker" = "Smoker",
                                                            "Age > 65" = "Age"),
                                                inline = TRUE),
                             submitButton("Submit")
                ),
                mainPanel(width = 12,
                          htmlOutput("predtext"),
                          plotlyOutput("IndSurvKM"),
                          tableOutput("IndPredTable")
                )
              )
      )
      
    )
  )
)




# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # load functions of interest
  # source("./Scripts/GetResultsVarSelect.R")
  source("./Scripts/MakeKM.R")
  # source("./Scripts/GetKMStuff.R")
  # source("./Scripts/PlotTree.R")
  # source("./Scripts/PredictIncoming.R")
  
  # load("FirstRun.Rdata")
  # load("RiskGroupsResults.Rdata")
  
  output$VariableHeader <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Prognostic performance", "</b></font> </u> </h3>") })
  output$RiskHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Histogram of predicted risk score", "</b></font> </u> </h4>")})
  output$CIHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Concordance index in predicting overall survival", "</b></font> </u> </h4>")})
  output$RefitHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Cox regression estimates and significance P-values", "</b></font> </u> </h4>")})
  output$EffectHeader <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Individual gene effect size and relative importance", "</b></font> </u> </h3>") })
  output$FreqHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Gene selection frequency", "</b></font> </u> </h4>")})
  output$VolcanoHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Interactive Volcano plot", "</b></font> </u> </h4>")})
  output$CommentPval <- renderText({paste("<font color=\"black\">","*Note that the pvalue is 0 here because it goes beyond the precision of the machine.", "</font>")})
  
  output$RiskHistogram <- renderPlot({FirstRun$RiskHistogram})
  output$RiskSummary <- renderTable({FirstRun$RiskScoreSummary}, rownames = TRUE)
  output$CI <- renderTable({FirstRun$ciSummary}, rownames = TRUE)
  output$RefitRisk <- renderTable({FirstRun$RiskRefit}, rownames = TRUE)
  output$ClinRefit <- renderTable({FirstRun$ClinRefitTable}, rownames = TRUE)
  output$influencePlot <- renderPlot({FirstRun$inflPlot})
  
  output$predRiskText <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Risk group stratification", "</b></font> </u> </h3>") })
  GetResultsReactive <- reactive({getResults(studyType = "Lung",
                                             method="LASSO",
                                             #clinical = input$AddClin,
                                             CNV = FALSE,
                                             OnlyCNV = FALSE,
                                             geneList = unlist(strsplit(input$GeneListRisk, split ="," )))})
  
  output$effectPlot <- renderPlotly({GetResultsReactive()$selectInflPlot})
  
  KMStuffReactive <- reactive({
    KMStuff(FirstRun$data.out,FirstRun$average.risk,
            FirstRun$topHits,4,
            c(0.25,0.75,0.9),
            geneList=unlist(strsplit(input$GeneListRisk, split ="," )))
  })
  output$KMText <- renderText({ paste("<h4> <u> <font color=\"black\"><b>","Kaplan-Meier plot of overall survival", "</b></font> </u> </h4>") })
  
  output$KM <- renderPlot(RiskGroupsResults$KM_Plot)
  output$SurvSum <- renderTable(RiskGroupsResults$SurvSum,rownames = TRUE)
  
  output$MutGroupText <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Mutation profiles by risk groups", "</b></font> </u> </h3>") })
  output$MutBPText <- renderText({ paste("<h4> <u> <font color=\"black\"><b>","Barplot of mutation frequency", "</b></font> </u> </h4>") })
  
  output$Mut <- renderPlot({
    print(RiskGroupsResults$mut_Plot)})
  
  ## Tab 2
  output$profiletext <- renderText({ paste("<h4> <u> <font color=\"black\"><b>","Piechart of most representative mutation profiles : ",
                                           KMStuffReactive()$GenesUsed, "</b></font> </u> </h4>") })
  output$ProfilePie <- renderPlotly({
    #if(input$ShowPies) 
    print(KMStuffReactive()$PieChart)
  })
  
  
  makePredictionsReactive <- reactive({
    predictIncomingPatient(mutGenes = unlist(strsplit(input$GeneList, split =",")),
                           clinical=c(input$Demographics),#,input$MetSite),
                           ClinRefit=FirstRun$ClinRefit,
                           time.type=FirstRun$time.type,
                           MD=FirstRun$MD,
                           LassoFits=FirstRun$LassoFits,
                           RiskScore=FirstRun$average.risk)
  })
  
  output$IndSurvKM <- renderPlotly({makePredictionsReactive()$IndSurvKM})
  output$IndPredTable <- renderTable({makePredictionsReactive()$IndPredTable},rownames = TRUE) 
  
}

shinyApp(ui, server)
