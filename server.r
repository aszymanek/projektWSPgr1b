library("shiny")
library("affy")
library("Biobase")
library("affydata")
library("DT")
library("gahgu95av2.db")

#limit przyjmowanego pliku
options(shiny.maxRequestSize = 500*1024^2)

#targets= readTargets("targets.txt")

shinyServer(function(input, output, session) {
  
########################################################## Wczytanie danych, przerobka i readAffy
  
  dataRaw=reactive({
    print(input$celfiles$datapath)
    print(getwd())
    #inputdata=ReadAffy(filenames=unlist(input$file$datapath))

    pData <- read.table(input$scanfile$datapath, header=TRUE, sep="\t")
    metadata <- data.frame(labelDescription=c("Anotacja","Tkanka","Probka","Nazwa pliku CEL"),row.names=colnames(pData))
    phenoData=new("AnnotatedDataFrame",data=pData, varMetadata=metadata)
    
############################# WCZYTANIE DANYCH MIKROMACIERZOWYCH #######
    data_raw=ReadAffy(filenames=input$celfiles$datapath,phenoData=phenoData) 
    data_raw@annotation=paste("ga",data_raw@annotation,sep="")
    data_raw@cdfName=paste("GA",data_raw@cdfName,sep="") 
    
    
    return(data_raw)
  })
  
  
  
  ###########################################################################
  
  ################################ korekcja tla i normalizacja metody
  bgCorrection <- reactive ({
    switch (input$bgc, 
            "Choose method" = NULL, 
            "Rma"= "rma",
            "Mas"= "mas"
    )
  })
  
  normal<- reactive({
    switch (input$norm, 
            "Choose method" = NULL, 
            "Quantiles" = "quantiles",
            "Loess" = "loess"
    )})
  
  
###################################################### Korekcja tla i normalizacja
  observeEvent(input$normBtn, {
    progress <- shiny::Progress$new()
    progress$set(message = "Computing data", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    expresionset()
    degradacja()})
  
  skorygowane=reactive({
    
    if (!is.null(bgCorrection()) && !is.null(normal()) ) {
      
      znormalizowane=expresso(dataRaw(),
                              bgcorrect.method=bgCorrection(),
                              normalize.method=normal(),
                              pmcorrect.method="pmonly",
                              summary.method="medianpolish")
      
    }
    return(znormalizowane)
    
  })
  

 skorygowane2 = reactive({
                     skorygowane3=exprs(skorygowane())
    return(skorygowane3)
  })
 
###################### Stworzenie obiektu ExpressionSet
 
 expresionset = reactive({
   
   expset=ExpressionSet(assayData=skorygowane2())
   dane = expset@assayData$exprs
   return(dane)
   
 })
 
######################################################  Degradacja RNA
 degradacja = reactive({
   deg = AffyRNAdeg(dataRaw(),log.it=TRUE)
   return(deg)
 })

  observe({
    if(is.null(input$celfiles)){
      updateSelectInput(session, "bgc", 
                        selected="Choose a Method")
      updateSelectInput(session, "norm", 
                        selected="Choose a Method")
    }
  })
  
  ####################################################### Wykresiki
  
  output$histo_przed<- renderPlot({
    
    hist(dataRaw(),main="Histogram przed normalizacjÄ…")
    
    
    dev.off()
  })
  
  
  
  output$box_przed<- renderPlot({
#     #pasek ladowania wykesu
#     progress <- shiny::Progress$new()
#     progress$set(message = "Computing data", value = 0)
#     # Close the progress when this reactive exits (even if there's an error)
#     on.exit(progress$close())
#     
    boxplot(dataRaw(),main="Wykres pudelkowy przed normalizacjÄ…")
    
    dev.off()
  })
 
  
  
 output$histo<- renderPlot({
#    #pasek ladowania wykesu
#    progress <- shiny::Progress$new()
#    progress$set(message = "Computing data", value = 0)
#    # Close the progress when this reactive exits (even if there's an error)
#    on.exit(progress$close())
#    
   hist(expresionset(),main="Normalized data distribution")
   
   dev.off()
 })
  
 output$degr<- renderPlot({
   #pasek ladowania wykesu
   progress <- shiny::Progress$new()
   progress$set(message = "Computing data", value = 0)
   # Close the progress when this reactive exits (even if there's an error)
   on.exit(progress$close())
   
   plotAffyRNAdeg(degradacja(), transform = "shift.scale", cols=NULL)
   
   dev.off()
 })
 
 
################################################ DownloadHandlers
 
 output$downPlot1=  downloadHandler(
   #specify the file name
   filename = function(){
     paste("Histogram_after_normalization.pdf")
   },
   #open the device -> creat the plot -> close the device
   content = function(file){
     
     pdf(file)
     hist(expresionset(),main="Normalized data distribution")
     dev.off()
     
     
   }
 )
 
 
 output$downPlot2=  downloadHandler(
   #specify the file name
   filename = function(){
     paste("Histogram_before_normalization.pdf")
   },
   #open the device -> creat the plot -> close the device
   content = function(file){
   
     pdf(file)
     hist(dataRaw(),main="Raw data distribution")
     dev.off()
     
     
   }
 )
 
 output$downPlot3=  downloadHandler(
   #specify the file name
   filename = function(){
     paste("Boxplot_before_normalization.pdf")
   },
   #open the device -> creat the plot -> close the device
   content = function(file){
     
     pdf(file)
     boxplot(dataRaw(),main="Raw data distribution")
     dev.off()
     
     
   }
 )
 
 ################ Wykresy ###################
 output$downPlot4=  downloadHandler(
   #specify the file name
   filename = function(){
     paste("DegradationRNA.pdf")
   },
   #open the device -> creat the plot -> close the device
   content = function(file){
     
     pdf(file)
     plotAffyRNAdeg(degradacja(), transform = "shift.scale", cols=NULL)
     dev.off()
     
   })
 
 
 ######################################### ANALIZA RÓ¯NICUJ¥CA 
 observeEvent(input$runstat,{
   
   calculationsdfg()
   calculationsdwng()
   calculationsupg()
   
   output$diffgenes = DT::renderDataTable(calculationsdfg())
   output$downldiff = downloadHandler(
     filename = "differentiating_genes.txt",
     
     content = function(file) {
       write.table(x = calculationsdfg(), file)
     })
   
   
   output$downgenes = DT::renderDataTable(datatable(calculationsdwng()))
   output$downldown = downloadHandler(
     filename = "downregulated_genes.txt",
     
     content = function(file) {
       write.table(x = calculationsdwng(), file)
     })
   
   
   output$upgenes = DT::renderDataTable(datatable(calculationsupg()))
   output$downlup = downloadHandler(
     filename = "upregulated_genes.txt",
     
     content = function(file) {
       write.table(x = calculationsupg(), file)
     })
   
   
   })
 
   calculationsdfg = reactive({
     if(input$difReg == T){
       if(input$normality == F){
         dfg = diffcase()
         }
       else
         {
           dfg = diffnormcase()
         }
       return(dfg)
     }
       
     })
   calculationsdwng = reactive({
     
     if(input$downReg == T){
       if(input$normality == F){
         
         dwng = downcase()
       }
       else
       {
         dwng = downnormcase()
       }
       return(dwng)
     }
     
     
   })
   calculationsupg = reactive({
     if(input$upReg == T){
       if(input$normality == F){
         upg = upcase()
       }
       else
       {
         upg = upnormcase()
       }
       return(upg)
     }
     
     
   })
   
   ############# Podzia³ grup 
   gr1class = reactive ({
   switch (input$group1, 
           "Choose group" = NULL,
           "Lung Adenocarcinoma" = "ADENO",
           "Squamous Cell Lung Carcinoma" = "SQUAMOUS",
           "Pulmonary Carcinoids" = "CARCINOID",
           "Small Cell Lung Cancer" = "SMALL CELL",
           "Healthy" = "NORMAL"
   )
 })
   gr1 = reactive({
   grupa1 = expresionset()[,dataRaw()@phenoData@data$CLASS == gr1class()]
   rownames(grupa1) = rownames(expresionset())
   colnames(grupa1) = dataRaw()@phenoData@data$scan[dataRaw()@phenoData@data$CLASS == gr1class()]
   return(grupa1)
 })
   gr2class = reactive ({
   switch (input$group2, 
           "Choose group" = NULL,
           "Lung Adenocarcinoma" = "ADENO",
           "Squamous Cell Lung Carcinoma" = "SQUAMOUS",
           "Pulmonary Carcinoids" = "CARCINOID",
           "Small Cell Lung Cancer" = "SMALL CELL",
           "Healthy" = "NORMAL"
   )
 })
   gr2 = reactive({
   grupa2 = expresionset()[,dataRaw()@phenoData@data$CLASS == gr2class()]
   rownames(grupa2) = rownames(expresionset())
   colnames(grupa2) = dataRaw()@phenoData@data$scan[dataRaw()@phenoData@data$CLASS == gr2class()]
   return(grupa2)
 })
   
   ################################# SKRYPTY UWZGLÊDNIAJ¥CE ROZK£AD
   normchck = reactive({
     pvalnorm = matrix(NaN, length(rownames(expresionset())), 2,
                       dimnames = list(c(rownames(expresionset())), 
                                       c("Group1", "Group 2")))
     g1 = gr1()
     g2 = gr2()
     
     for (i in 1:nrow(g1)){
       pvalnorm[i,1] = shapiro.test(g1[i,])$p.value
       pvalnorm[i,2] = shapiro.test(g2[i,])$p.value
     }
     
     return(pvalnorm)
   })
   ## geny ró¿nicuj¹ce
   diffnormcase = eventReactive(input$difReg == T,{
     finddiffnorm = reactive({
       progress <- shiny::Progress$new()
       progress$set(message = "Computing data", value = 0)
       # Close the progress when this reactive exits (even if there's an error)
       on.exit(progress$close())
       
       pvaldiff = matrix(NaN, length(rownames(expresionset())), 1,
                         dimnames = list(c(rownames(expresionset())), 
                                         c("P-value")))
       pnorm = normchck()
       g1 = gr1()
       g2 = gr2()
       
       for (i in 1:nrow(pnorm)){
         if (pnorm[i,1] < input$pvaltresh && pnorm[i,2] < input$pvaltresh){
           pvaldiff[i,1] = t.test(x=g1[i,], y=g2[i,], alternative="two.sided")$p.value
         }
         else{
           pvaldiff[i,1] = wilcox.test(x=g1[i,], y=g2[i,], alternative="two.sided")$p.value
         }
         
       }
       
       return(pvaldiff)
     })
     infodiffnorm = reactive({
       
       
       probesdiff = finddiffnorm()
       nanmat = matrix(NaN, nrow(probesdiff), 1)
       
       pvals = nanmat
       gsymb = nanmat
       entrez = nanmat
       gname = nanmat
       
       pvals = unlist(probesdiff[probesdiff < input$pvaltresh])
       
       progress <- shiny::Progress$new()
       progress$set(message = "Summarizing informations", value = 0)
       # Close the progress when this reactive exits (even if there's an error)
       on.exit(progress$close())
       
       ### Getting gene symbols ##
       x_sym <- gahgu95av2SYMBOL
       # Get the probe identifiers that are mapped to gene alias
       mapped_probes_sym <- mappedkeys(x_sym)
       # Convert to a list
       xx_sym <- as.list(x_sym[mapped_probes_sym])
       if(length(xx_sym) > 0) {
         gsymb = unlist(xx_sym[probesdiff < input$pvaltresh])
       }
       
       ### Getting Entrez IDs ##
       x_ent <- gahgu95av2ENTREZID
       # Get the probe identifiers that are mapped to gene alias
       mapped_probes_ent <- mappedkeys(x_ent)
       # Convert to a list
       xx_ent <- as.list(x_ent[mapped_probes_ent])
       if(length(xx_ent) > 0) {
         entrez = unlist(xx_ent[probesdiff < input$pvaltresh])
       }
       
       ### Getting genes names ##
       x_nam <- gahgu95av2GENENAME
       # Get the probe identifiers that are mapped to gene alias
       mapped_probes_nam <- mappedkeys(x_nam)
       # Convert to a list
       xx_nam <- as.list(x_nam[mapped_probes_nam])
       if(length(xx_nam) > 0) {
         gname = unlist(xx_nam[probesdiff < input$pvaltresh])
       }
       
       diffsuminf = data.frame(data = cbind(pvals, gsymb, entrez, gname),
                               stringsAsFactors = F, row.names = NULL)
       rownames(diffsuminf) = rownames(probesdiff)[probesdiff<input$pvaltresh]
       colnames(diffsuminf) = c("P-value", "Gene Symbol", "Entrez ID", "Gene name")
       
       diffsuminf = diffsuminf[ order(diffsuminf[,1]), ]
       return(diffsuminf)
       
     })
     
     return(infodiffnorm())
     
   })
   ## down-regulated
   downnormcase = eventReactive(input$downReg == T,{
     finddownnorm = reactive({
       progress <- shiny::Progress$new()
       progress$set(message = "Computing data", value = 0)
       # Close the progress when this reactive exits (even if there's an error)
       on.exit(progress$close())
       
       pvaldiff = matrix(NaN, length(rownames(expresionset())), 1,
                         dimnames = list(c(rownames(expresionset())), 
                                         c("P-value")))
       pnorm = normchck()
       g1 = gr1()
       g2 = gr2()
       
       for (i in 1:nrow(pnorm)){
         if (pnorm[i,1] < input$pvaltresh && pnorm[i,2] < input$pvaltresh){
           pvaldiff[i,1] = t.test(x=g1[i,], y=g2[i,], alternative="less")$p.value
         }
         else{
           pvaldiff[i,1] = wilcox.test(x=g1[i,], y=g2[i,], alternative="less")$p.value
         }
         
       }
       
       return(pvaldiff)
     })
     infodownnorm = reactive({
       
       
       probesdiff = finddownnorm()
       nanmat = matrix(NaN, nrow(probesdiff), 1)
       
       pvals = nanmat
       gsymb = nanmat
       entrez = nanmat
       gname = nanmat
       
       pvals = unlist(probesdiff[probesdiff < input$pvaltresh])
       
       progress <- shiny::Progress$new()
       progress$set(message = "Summarizing informations", value = 0)
       # Close the progress when this reactive exits (even if there's an error)
       on.exit(progress$close())
       
       ### Getting gene symbols ##
       x_sym <- gahgu95av2SYMBOL
       # Get the probe identifiers that are mapped to gene alias
       mapped_probes_sym <- mappedkeys(x_sym)
       # Convert to a list
       xx_sym <- as.list(x_sym[mapped_probes_sym])
       if(length(xx_sym) > 0) {
         gsymb = unlist(xx_sym[probesdiff < input$pvaltresh])
       }
       
       ### Getting Entrez IDs ##
       x_ent <- gahgu95av2ENTREZID
       # Get the probe identifiers that are mapped to gene alias
       mapped_probes_ent <- mappedkeys(x_ent)
       # Convert to a list
       xx_ent <- as.list(x_ent[mapped_probes_ent])
       if(length(xx_ent) > 0) {
         entrez = unlist(xx_ent[probesdiff < input$pvaltresh])
       }
       
       ### Getting genes names ##
       x_nam <- gahgu95av2GENENAME
       # Get the probe identifiers that are mapped to gene alias
       mapped_probes_nam <- mappedkeys(x_nam)
       # Convert to a list
       xx_nam <- as.list(x_nam[mapped_probes_nam])
       if(length(xx_nam) > 0) {
         gname = unlist(xx_nam[probesdiff < input$pvaltresh])
       }
       
       diffsuminf = data.frame(data = cbind(pvals, gsymb, entrez, gname),
                               stringsAsFactors = F, row.names = NULL)
       rownames(diffsuminf) = rownames(probesdiff)[probesdiff<input$pvaltresh]
       colnames(diffsuminf) = c("P-value", "Gene Symbol", "Entrez ID", "Gene name")
       
       diffsuminf = diffsuminf[ order(diffsuminf[,1]), ]
       return(diffsuminf)
       
     })
     
     return(infodownnorm())
   })
   ## up-regulated
   upnormcase = eventReactive(input$upReg == T,{
     findupnorm = reactive({
       progress <- shiny::Progress$new()
       progress$set(message = "Computing data", value = 0)
       # Close the progress when this reactive exits (even if there's an error)
       on.exit(progress$close())
       
       pvaldiff = matrix(NaN, length(rownames(expresionset())), 1,
                         dimnames = list(c(rownames(expresionset())), 
                                         c("P-value")))
       pnorm = normchck()
       g1 = gr1()
       g2 = gr2()
       
       for (i in 1:nrow(pnorm)){
         if (pnorm[i,1] < input$pvaltresh && pnorm[i,2] < input$pvaltresh){
           pvaldiff[i,1] = t.test(x=g1[i,], y=g2[i,], alternative="greater")$p.value
         }
         else{
           pvaldiff[i,1] = wilcox.test(x=g1[i,], y=g2[i,], alternative="greater")$p.value
         }
         
       }
       
       return(pvaldiff)
     })
     infoupnorm = reactive({
       
       
       probesdiff = findupnorm()
       nanmat = matrix(NaN, nrow(probesdiff), 1)
       
       pvals = nanmat
       gsymb = nanmat
       entrez = nanmat
       gname = nanmat
       
       pvals = unlist(probesdiff[probesdiff < input$pvaltresh])
       
       progress <- shiny::Progress$new()
       progress$set(message = "Summarizing informations", value = 0)
       # Close the progress when this reactive exits (even if there's an error)
       on.exit(progress$close())
       
       ### Getting gene symbols ##
       x_sym <- gahgu95av2SYMBOL
       # Get the probe identifiers that are mapped to gene alias
       mapped_probes_sym <- mappedkeys(x_sym)
       # Convert to a list
       xx_sym <- as.list(x_sym[mapped_probes_sym])
       if(length(xx_sym) > 0) {
         gsymb = unlist(xx_sym[probesdiff < input$pvaltresh])
       }
       
       ### Getting Entrez IDs ##
       x_ent <- gahgu95av2ENTREZID
       # Get the probe identifiers that are mapped to gene alias
       mapped_probes_ent <- mappedkeys(x_ent)
       # Convert to a list
       xx_ent <- as.list(x_ent[mapped_probes_ent])
       if(length(xx_ent) > 0) {
         entrez = unlist(xx_ent[probesdiff < input$pvaltresh])
       }
       
       ### Getting genes names ##
       x_nam <- gahgu95av2GENENAME
       # Get the probe identifiers that are mapped to gene alias
       mapped_probes_nam <- mappedkeys(x_nam)
       # Convert to a list
       xx_nam <- as.list(x_nam[mapped_probes_nam])
       if(length(xx_nam) > 0) {
         gname = unlist(xx_nam[probesdiff < input$pvaltresh])
       }
       
       diffsuminf = data.frame(data = cbind(pvals, gsymb, entrez, gname),
                               stringsAsFactors = F, row.names = NULL)
       rownames(diffsuminf) = rownames(probesdiff)[probesdiff<input$pvaltresh]
       colnames(diffsuminf) = c("P-value", "Gene Symbol", "Entrez ID", "Gene name")
       
       diffsuminf = diffsuminf[ order(diffsuminf[,1]), ]
       return(diffsuminf)
       
     })
     
     return(infoupnorm())
   })
     
   ############################# SKRYPTY NIEPARAMETRYCZNE   
     
   ## geny ró¿nicuj¹ce
   diffcase = eventReactive(input$difReg,{
         finddiff = reactive({
           progress <- shiny::Progress$new()
           progress$set(message = "Computing data", value = 0)
           # Close the progress when this reactive exits (even if there's an error)
           on.exit(progress$close())
           
           pvaldiff = matrix(NaN, length(rownames(expresionset())), 1,
                             dimnames = list(c(rownames(expresionset())), 
                                             c("P-value")))
           g1 = gr1()
           g2 = gr2()
           
           for (i in 1:nrow(g1)){
             pvaldiff[i,1] = wilcox.test(x=g1[i,], y=g2[i,], alternative="two.sided")$p.value
           }
           
           
           return(pvaldiff)
         })
         infodiff = reactive({
           
           probesdiff = finddiff()
           nanmat = matrix(NaN, nrow(probesdiff), 1)
           
           pvals = nanmat
           gsymb = nanmat
           entrez = nanmat
           gname = nanmat
           
           pvals = unlist(probesdiff[probesdiff < input$pvaltresh])
           
           progress <- shiny::Progress$new()
           progress$set(message = "Summarizing informations", value = 0)
           # Close the progress when this reactive exits (even if there's an error)
           on.exit(progress$close())
           
           ### Getting gene symbols ##
           x_sym <- gahgu95av2SYMBOL
           # Get the probe identifiers that are mapped to gene alias
           mapped_probes_sym <- mappedkeys(x_sym)
           # Convert to a list
           xx_sym <- as.list(x_sym[mapped_probes_sym])
           if(length(xx_sym) > 0) {
             gsymb = unlist(xx_sym[probesdiff < input$pvaltresh])
           }
           
           ### Getting Entrez IDs ##
           x_ent <- gahgu95av2ENTREZID
           # Get the probe identifiers that are mapped to gene alias
           mapped_probes_ent <- mappedkeys(x_ent)
           # Convert to a list
           xx_ent <- as.list(x_ent[mapped_probes_ent])
           if(length(xx_ent) > 0) {
             entrez = unlist(xx_ent[probesdiff < input$pvaltresh])
           }
           
           ### Getting genes names ##
           x_nam <- gahgu95av2GENENAME
           # Get the probe identifiers that are mapped to gene alias
           mapped_probes_nam <- mappedkeys(x_nam)
           # Convert to a list
           xx_nam <- as.list(x_nam[mapped_probes_nam])
           if(length(xx_nam) > 0) {
             gname = unlist(xx_nam[probesdiff < input$pvaltresh])
           }
           
           diffsuminf = data.frame(data = cbind(pvals, gsymb, entrez, gname),
                                   stringsAsFactors = F, row.names = NULL)
           rownames(diffsuminf) = rownames(probesdiff)[probesdiff<input$pvaltresh]
           colnames(diffsuminf) = c("P-value", "Gene Symbol", "Entrez ID", "Gene name")
           
           diffsuminf = diffsuminf[ order(diffsuminf[,1]), ]
           return(diffsuminf)
           
         })
         
         return(infodiff())
       })
   ## down-regulated
   downcase = eventReactive(input$downReg,{
         finddown = reactive({
           progress <- shiny::Progress$new()
           progress$set(message = "Computing data", value = 0)
           # Close the progress when this reactive exits (even if there's an error)
           on.exit(progress$close())
           
           pvaldiff = matrix(NaN, length(rownames(expresionset())), 1,
                             dimnames = list(c(rownames(expresionset())), 
                                             c("P-value")))
           g1 = gr1()
           g2 = gr2()
           
           for (i in 1:nrow(g1)){
             pvaldiff[i,1] = wilcox.test(x=g1[i,], y=g2[i,], alternative="less")$p.value
           }
           
           
           return(pvaldiff)
         })
         infodown = reactive({
           
           
           probesdiff = finddown()
           nanmat = matrix(NaN, nrow(probesdiff), 1)
           
           pvals = nanmat
           gsymb = nanmat
           entrez = nanmat
           gname = nanmat
           
           pvals = unlist(probesdiff[probesdiff < input$pvaltresh])
           
           progress <- shiny::Progress$new()
           progress$set(message = "Summarizing informations", value = 0)
           # Close the progress when this reactive exits (even if there's an error)
           on.exit(progress$close())
           
           ### Getting gene symbols ##
           x_sym <- gahgu95av2SYMBOL
           # Get the probe identifiers that are mapped to gene alias
           mapped_probes_sym <- mappedkeys(x_sym)
           # Convert to a list
           xx_sym <- as.list(x_sym[mapped_probes_sym])
           if(length(xx_sym) > 0) {
             gsymb = unlist(xx_sym[probesdiff < input$pvaltresh])
           }
           
           ### Getting Entrez IDs ##
           x_ent <- gahgu95av2ENTREZID
           # Get the probe identifiers that are mapped to gene alias
           mapped_probes_ent <- mappedkeys(x_ent)
           # Convert to a list
           xx_ent <- as.list(x_ent[mapped_probes_ent])
           if(length(xx_ent) > 0) {
             entrez = unlist(xx_ent[probesdiff < input$pvaltresh])
           }
           
           ### Getting genes names ##
           x_nam <- gahgu95av2GENENAME
           # Get the probe identifiers that are mapped to gene alias
           mapped_probes_nam <- mappedkeys(x_nam)
           # Convert to a list
           xx_nam <- as.list(x_nam[mapped_probes_nam])
           if(length(xx_nam) > 0) {
             gname = unlist(xx_nam[probesdiff < input$pvaltresh])
           }
           
           diffsuminf = data.frame(data = cbind(pvals, gsymb, entrez, gname),
                                   stringsAsFactors = F, row.names = NULL)
           rownames(diffsuminf) = rownames(probesdiff)[probesdiff<input$pvaltresh]
           colnames(diffsuminf) = c("P-value", "Gene Symbol", "Entrez ID", "Gene name")
           
           diffsuminf = diffsuminf[ order(diffsuminf[,1]), ]
           return(diffsuminf)
           
         })
         
         return(infodown())
       })
   ## up-regulated
   upcase = eventReactive(input$upReg,{
         findup = reactive({
           progress <- shiny::Progress$new()
           progress$set(message = "Computing data", value = 0)
           # Close the progress when this reactive exits (even if there's an error)
           on.exit(progress$close())
           
           pvaldiff = matrix(NaN, length(rownames(expresionset())), 1,
                             dimnames = list(c(rownames(expresionset())), 
                                             c("P-value")))
           g1 = gr1()
           g2 = gr2()
           
           for (i in 1:nrow(g1)){
             pvaldiff[i,1] = wilcox.test(x=g1[i,], y=g2[i,], alternative="greater")$p.value
           }
           
           
           return(pvaldiff)
         })
         infoup = reactive({
           
           
           probesdiff = findup()
           nanmat = matrix(NaN, nrow(probesdiff), 1)
           
           pvals = nanmat
           gsymb = nanmat
           entrez = nanmat
           gname = nanmat
           
           pvals = unlist(probesdiff[probesdiff < input$pvaltresh])
           
           progress <- shiny::Progress$new()
           progress$set(message = "Summarizing informations for", value = 0)
           # Close the progress when this reactive exits (even if there's an error)
           on.exit(progress$close())
           
           ### Getting gene symbols ##
           x_sym <- gahgu95av2SYMBOL
           # Get the probe identifiers that are mapped to gene alias
           mapped_probes_sym <- mappedkeys(x_sym)
           # Convert to a list
           xx_sym <- as.list(x_sym[mapped_probes_sym])
           if(length(xx_sym) > 0) {
             gsymb = unlist(xx_sym[probesdiff < input$pvaltresh])
           }
           
           ### Getting Entrez IDs ##
           x_ent <- gahgu95av2ENTREZID
           # Get the probe identifiers that are mapped to gene alias
           mapped_probes_ent <- mappedkeys(x_ent)
           # Convert to a list
           xx_ent <- as.list(x_ent[mapped_probes_ent])
           if(length(xx_ent) > 0) {
             entrez = unlist(xx_ent[probesdiff < input$pvaltresh])
           }
           
           ### Getting genes names ##
           x_nam <- gahgu95av2GENENAME
           # Get the probe identifiers that are mapped to gene alias
           mapped_probes_nam <- mappedkeys(x_nam)
           # Convert to a list
           xx_nam <- as.list(x_nam[mapped_probes_nam])
           if(length(xx_nam) > 0) {
             gname = unlist(xx_nam[probesdiff < input$pvaltresh])
           }
           
           diffsuminf = data.frame(data = cbind(pvals, gsymb, entrez, gname),
                                   stringsAsFactors = F, row.names = NULL)
           rownames(diffsuminf) = rownames(probesdiff)[probesdiff<input$pvaltresh]
           colnames(diffsuminf) = c("P-value", "Gene Symbol", "Entrez ID", "Gene name")
           
           diffsuminf = diffsuminf[ order(diffsuminf[,1]), ]
           return(diffsuminf)
           
         })
         
         return(infoup())
       })
       
       
 
 })


 
