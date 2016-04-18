

#limit przyjmowanego pliku
options(shiny.maxRequestSize = 500*1024^2)

## Tutaj ładujemy pliki - prosze nic nie zmieniać
#targets= readTargets("targets.txt")

shinyServer(function(input, output, session) {
  
  
  normalizacja=function(){
    print(input$file$datapath)
    print(input$file)
    print(getwd())
    inputdata=ReadAffy(filenames=unlist(input$file$datapath))
    
    return(inputdata)
  }
  
  
  
  ###########################################################################
  
  ################################korekcja tła
  bgCorrection <- reactive ({
    switch (input$bgc, 
            "Choose method" = NULL, 
            "wybierz"= "metode",
            "wybierz"= "metode"
    )
  })
  
  normal<- reactive({
    switch (input$norm, 
            "Choose method" = NULL, 
            "wybierz" = "metode",
            "wybierz" = "metode"
    )})
  
  

  
  

  

    

  
  

  
  #############################Wykresiki
  output$wykres1<- renderPlot({
    
  })
  
  output$wykres2<- renderPlot({
    
  })
  
  output$wykres3<- renderPlot({
   
  })
  
  output$wykres4<- renderPlot({
    
  })
  
  
})