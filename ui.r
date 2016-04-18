shinyUI(fluidPage(
  
  titlePanel("Nasz projekt"),
  
  sidebarLayout(
    sidebarPanel(     
      helpText("Input your microarray data"),
      fileInput('file', 'Choose file to upload', multiple=TRUE,
                accept = c('text/csv',
                  '.CEL')),
      conditionalPanel(
        condition="output.uploadedFile",
        helpText("choose a method")),
      selectInput("bgc","Background correction method:",
                  choices=c("Choose method","method1","method2")),
      selectInput("norm", "Normalization method:",
                  choices=c("Choose method","method_1","method_2"))),
    
    
    
    
    mainPanel(
      tabsetPanel(id ="maintabset",
                  tabPanel("Describe", 
                           helpText(p("To jest nasza piekna aplikacja i tu bedzie jej opis <3", align="center"))),
                            
                           
                  
                  tabPanel("Zakladka1", 
                           conditionalPanel(
                             condition="!output.histo",
                             tags$p("")),
                           conditionalPanel(
                             condition="output.histo",
                             tags$p("")), 
                           plotOutput("histo")
                        ),
                  
                  tabPanel("Zakladka2", 
                           conditionalPanel(
                             condition="!output.histo_przed",
                             tags$p("")),
                           conditionalPanel(
                             condition="output.histo_przed",
                             tags$p("")), 
                           plotOutput("histo_przed")
                  ),
                  
                  tabPanel("Zakladka3", 
                           conditionalPanel(
                             condition="!output.box_przed",
                             tags$p("")),
                           conditionalPanel(
                             condition="output.box_przed",
                             tags$p("")), 
                           plotOutput("box_przed")
                  ),
                  
                  tabPanel("Zakladka4", 
                           conditionalPanel(
                             condition="!output.degr",
                             tags$p("")),
                           conditionalPanel(
                             condition="output.degr",
                             tags$p("")), 
                           plotOutput("degr")
                  )
             )
      
            )
    )))