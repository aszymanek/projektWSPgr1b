shinyUI(fluidPage(
  
  titlePanel("Microarray data analysis"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(id = "sideset",
                  tabPanel("Data upload and preprocessing",
                           helpText("Input your Affymetrix microarray data"),
                           fileInput('celfiles', 'Choose file to upload', multiple=TRUE,
                                     accept = c('text/csv',
                                                '.CEL')),
                           fileInput('scanfile', 'Choose decription file', multiple=F,
                                     accept = c('text/csv',
                                                '.txt')),
                           conditionalPanel(
                             condition="output.uploadedFile",
                             helpText("choose a method")),
                           selectInput("bgc","Background correction method:",
                                       choices=c("Choose method","Mas","Rma")),
                           selectInput("norm", "Normalization method:",
                                       choices=c("Choose method","Quantiles","Loess")),
                           actionButton("normBtn", "Normalize")),
                  
                  ## panel selekcji cech ##
                  tabPanel("Analysis pipeline",
                           helpText("Choose groups to compare"),
                           conditionalPanel(
                             condition="output.uploadedFile",
                             helpText("choose first group")),
                           selectInput("group1","First group",
                                       choices=c("Choose group","Lung Adenocarcinoma","Squamous Cell Lung Carcinoma",
                                                 "Pulmonary Carcinoids", "Small Cell Lung Cancer", "Healthy")),
                           selectInput("group2", "Second group:",
                                       choices=c("Choose group","Lung Adenocarcinoma","Squamous Cell Lung Carcinoma",
                                                 "Pulmonary Carcinoids", "Small Cell Lung Cancer", "Healthy")),
                           checkboxInput("difReg", "Look for differentialy regulated genes"),
                           checkboxInput("downReg", "Look for down-regulated genes"),
                           checkboxInput("upReg", "Look for up-regulated genes"),
                           sliderInput("pvaltresh", "P-value treshold", 
                                       min = 0.01, max = 0.1, step = 0.01, value = 0.05),
                           checkboxInput("normality", "Check distribution normality", value = FALSE),
                           actionButton("runstat", "Run analysis"))
      )
                  ),
    
    mainPanel(
      tabsetPanel(id ="maintabset",
                  tabPanel("About app", 
                           div(img(src="team.png", height = 403, width = 588),style="text-align: center;"),
                           helpText(p("This is a web application which allows Affymetrix expression microarray data analysis.
                                      On the left panel you can load your data, and select a method for background correction and normalization and perform preprocessing.
                                      Next, in separate tab, you have to choose analysis pipeline - define groups, methods and treshold for analysis. You can also save you result plots and gene lists, which you can find in main panels tabs.", align="center"))),
                  
                  tabPanel("Preprocessing and quality control",
                           tabsetPanel(id="nestedpanel",
                                       tabPanel("Raw data",
                                                tabsetPanel(id="rawplots",
                                                            tabPanel("Histogram",
                                                                     conditionalPanel(
                                                                       condition="!output.histo_przed",
                                                                       tags$p("")),
                                                                     conditionalPanel(
                                                                       condition="output.histo_przed",
                                                                       tags$p("")),
                                                                     plotOutput("histo_przed"),
                                                                     br(),
                                                                     downloadButton("downPlot2", label = "Download Plot")),
                                                            tabPanel("Boxplot before normalization",
                                                                     conditionalPanel(
                                                                       condition="!output.box_przed",
                                                                       tags$p("")),
                                                                     conditionalPanel(
                                                                       condition="output.box_przed",
                                                                       tags$p("")),
                                                                     plotOutput("box_przed"),
                                                                     br(),
                                                                     downloadButton("downPlot3", label = "Download Plot")))),
                                       
                                       tabPanel("Normalized data",
                                                conditionalPanel(
                                                  condition="!output.histo",
                                                  tags$p("")),
                                                conditionalPanel(
                                                  condition="output.histo",
                                                  tags$p("")),
                                                plotOutput("histo"),
                                                br(),
                                                downloadButton("downPlot1", label = "Download Plot")),
                                       
                                       tabPanel("RNA Degradation",
                                                conditionalPanel(
                                                  condition="!output.degr",
                                                  tags$p("")),
                                                conditionalPanel(
                                                  condition="output.degr",
                                                  tags$p("")),
                                                  plotOutput("degr"),
                                                  br(),
                                                  downloadButton("downPlot4", label = "Download Plot")))),
                  
                  tabPanel("Differential analysis",
                           tabsetPanel(id = "diffan",
                                       tabPanel("Differentially regulated genes",
                                                DT::dataTableOutput("diffgenes"),
                                                downloadButton('downldiff', 'Save results')),
                                       tabPanel("Down-regulated genes",
                                                DT::dataTableOutput("downgenes"),
                                                downloadButton('downldown', 'Save results')),
                                       tabPanel("Up-regulated genes",
                                                DT::dataTableOutput("upgenes"),
                                                downloadButton('downlup', 'Save results'))))
                  )
      )
  )
  ))
