#   ESKAPE Act PLUS - Activation Analysis for ESKAPE Pathogens 
#   and other Prokaryotes Labs Usually Study
#   Copyright (C) <2022>  <Katja Koeppen>

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#    General Public License for more details <http://www.gnu.org/licenses/>.

#    If this code is helpful to you, please cite our related publication:
#    https://pubmed.ncbi.nlm.nih.gov/36259735/

library(shiny)
library(shinyjs)

ui <- shinyUI(fluidPage(
    useShinyjs(),
    ## HSA genome UPDATED ++ ##
    #titlePanel(div("ESKAPE Act PLUS - Activation Analysis for ESKAPE Pathogens 
      #and other Prokaryotes Labs Usually Study", style="color:blueviolet")),
    titlePanel(div("HSAPath - Activation Analysis for Human Cellular Pathways", 
                   style="color:blueviolet")),
    ## HSA genome UPDATED -- ##

    
    fluidRow(
        column(3,
               br(),
               wellPanel(# Input: Select a file ----
                         fileInput("file1", "Upload CSV file:", accept = ".csv"),

                         tags$small("Required file format: csv file needs to contain 
                             UniProt stable entry identifier in column 1 and 
                             log2 fold changes between two treatment groups in column 2."),
                         hr(),
                         ## HSA genome UPDATED ++ ##
                         # replace pseudomons example with homo sapien example data
                         # data generated from random shuffle of RNAseq data found at:
                         # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208637
                         
                         #downloadButton('Example', "Pseudomonas PA14 example data", class = "wrap"),
                         downloadButton('Example','HSA example data',class='wrap'),
                         ## HSA genome UPDATED -- ##
                         tags$head(tags$style(HTML(".wrap { white-space: normal; word-break: keep-all;  }"))),
                       
                         hr(),
                         uiOutput('error'),
                         tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }",
                                         '#Text {color: red;}')),
                         
                         
                         # Horizontal line ----
                         tags$hr(),
                         
                         # Input: Select strain ----
                         ## HSA genome UPDATED ++ ##
                         # bacterial strains not relevant for human genome
                         #selectInput('strain', label = "Select strain:", 
                                     #choices = c("Acinetobacter baumanii AYE" = "aby",
                                    #             "Acinetobacter baumanii MDR-ZJ06" = "abz",
                                    #             "Bacteroides fragilis NCTC 9343" = "bfs",
                                    #             "Bacteroides ovatus 3725 D1 iv" = "boa",
                                    #             "Bacteroides thetaiotaomicron VPI-5482" = "bth",
                                     #            "Clostridioides difficile 630" = "cdf",
                                    #             "Clostridioides difficile R20291" = "cdl",
                                    #             "Enterobacter cloacae ATCC 13047" = "enc",
                                    #             "Enterococcus faecium DO" = "efu",
                                    #             "Escherichia coli BL21(DE3)" = "ebd",
                                    #             "Escherichia coli K-12 MG1655" = "eco",
                                    #             "Escherichia coli K-12 W3110" = "ecj",
                                    #             "Escherichia coli O157:H7 EDL933" = "ece",
                                    #             "Escherichia coli O157:H7 Sakai" = "ecs",
                                    #             "Klebsiella pneumoniae MGH 78578" = "kpn",
                                    #             "Prevotella melaninogenica ATCC 25845" = "pmz",
                                    #             "Pseudomonas aeruginosa PA14" = "pau",
                                    #             "Pseudomonas aeruginosa PAO1" = "pae",
                                    #             "Staphylococcus aureus COL" = "sac",
                                    #             "Staphylococcus aureus NCTC8325" = "sao",
                                    #             "Staphylococcus aureus Newman" = "sae",
                                    #             "Staphylococcus aureus USA300 FPR3757" = "saa",
                                    #             "Streptococcus sanguinis SK36" = "ssa")),
                        ## HSA genome UDPATED -- ## 
                         
                         tags$hr(),
                         
                        ## HSA genome UDPATED ++ ##
                        # remove information about strain selection
                         tags$small("Upload a csv file with UniProt protein identifiers and log2 fold changes,
                          then click on 'Run Activation Analysis' that will appear below."),
                         #tags$small("Upload a csv file with UniProt protein identifiers and log2 fold changes,
                                  #  select a strain, then click on 'Run Activation Analysis' that will appear below."),
                        ## HSA genome UDPATED -- ##

                         tags$hr(),
                         
                         uiOutput("action"),
                         
                         tags$hr(),
                         
                         uiOutput("ui"),
                         
                         tags$hr(),
                         
                         downloadButton('Documentation', 'User Manual')
               )),
        
        column(8,  
               br(),
               mainPanel(
                   tabsetPanel(id = 'Results_Output', selected = 'Summary',
                     tabPanel('Summary',
                              tags$hr(),
                              textOutput('KEGG_match'),
                              tags$hr(),
                              textOutput('KEGG_res'),
                              tags$style(type='text/css', '#KEGG_res{color: red;}'),
                              tags$hr(),
                              textOutput('GO_match'),
                              tags$hr(),
                              textOutput('GO_res'),
                              tags$style(type='text/css', '#GO_res{color: red;}')),
                     tabPanel("KEGG Graph",
                                plotOutput("KEGG_plot")),
                     tabPanel("KEGG Table",
                                DT::DTOutput("KEGG_table")),
                     tabPanel("GO Graph",
                                plotOutput("GO_plot")),
                     tabPanel("GO Table",
                                DT::DTOutput("GO_table"))
                       ))
        )#,
        
        ## HSA genome UPDATED ++ ##
        # removed image specific to ESKAPE Act Plus
        #column(1,
        #       img(src="Logo.jpeg", align = "right", height='180px', width='250px'),
        #       img(src="ESKAPE_crop.gif", align = "right", height='600px', width='200px')
        #)
        ## HSA genome UPDATED -- ##
    )
)
)

server <- function(input, output, session) {
  
  First.wd <- getwd()

  # Load lists with KEGG and GO info
  ## HSA genome UPDATED ++ ##
  # RData file created by Create_RData_Object.R
  load(file = "KEGG_GO.Rdata" )
  ## HSA genome UPDATED -- ##
  
  output$Documentation <- downloadHandler(
    ## HSA genome UPDATED ++ ##
    # updated manual for HSAPath
    filename = "HSAPath_User_Manual.pdf",
    content <- function(file) {
      file.copy(from = paste(First.wd, "/", "HSAPath_User_Manual.pdf", sep = ""), to = file)
    }) 
  ## HSA genome UPDATED -- ##
  
  output$Example <- downloadHandler(
    #filename = "pau_Result.csv",
    filename = "HSAPath_demo.csv",
    content <- function(file) {
      #file.copy(from = paste(First.wd, "/", "pau_Result.csv", sep = ""), to = file)
      file.copy(from = paste(First.wd, "/", "HSAPath_demo.csv", sep = ""), to = file)
    }) 
  
  User_data <- eventReactive(input$file1, {
    file <- input$file1
    
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(file$datapath, header = FALSE, stringsAsFactors = FALSE)
  })
  
  observeEvent(input$file1, {
    
    if (dim(User_data())[2] > 2){
      output$Text <- renderText("Error - file has more than 2 columns!")
    } else if (dim(User_data())[2] < 2){
      output$Text <- renderText("Error - file has fewer than 2 columns!")
    } else {
      output$Text <- renderText("File has the correct number of columns.")
      output$action <- renderUI({
        actionButton('run', "Run Activation Analysis",
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
      })
    }
    
    output$error <- renderUI({verbatimTextOutput('Text')})
    
  })
  
  # Run Activation Analysis
  observeEvent(input$run, {
    
      Results <- User_data()
      
      ## HSA genome UPDATED ++ ##
      # conform sub-setting to new kegg.df structure (no need to subset by strain)
      #KEGG_Prots <- KEGG[[input$strain]]$Entry[KEGG[[input$strain]]$Entry %in% Results$V1]
      #KEGG_Paths <- KEGG[[input$strain]]$KEGG_paths[KEGG[[input$strain]]$Entry %in% Results$V1]
      All_Prots<-kegg.df$Entry
      KEGG_Prots <- All_Prots[All_Prots %in% Results$V1]
      KEGG_Paths <- kegg.df$pathway[kegg.df$Entry %in% Results$V1]
      ## HSA genome UPDATED -- ##
      
      output$KEGG_match <- renderText(paste(length(unique(KEGG_Prots)),
                           "of user-provided identifiers matched to one of",
                           ## HSA genome UPDATED ++ ##
                           # updated for new kegg data object structure (no need to subset by strain)
                           #length(unique(KEGG[[input$strain]]$Entry)),
                           length(unique(All_Prots)),
                           ## HSA genome UPDATED -- ##
                           "possible identifiers for KEGG pathway analysis", sep = " "))
      
      if (length(KEGG_Prots) == 0) {
        showModal(modalDialog(title = "ID error - double check selected strain", 
        HTML(paste('The identifiers you provided do not 
        match the UniProt identifiers of the selected strain. <br> 
        These are the first 5 identifiers you entered:',
        paste(Results$V1[1:5], collapse = " "), '<br> Example UniProt identifiers for the 
        selected strain look like this:', 
        paste(unique(KEGG[[input$strain]]$Entry)[1:10], collapse = " ")), sep = " "), 
        size = "l", footer = actionButton('close_ID_error', label = "OK")))
        
        observeEvent(input$close_ID_error, {
          removeModal()
          output$action <- renderUI({
          actionButton('run', "Run Activation Analysis",
          style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
          })
        })
        
    
        # checking to make sure FCs are numbers
      } else if (is.numeric(Results$V2) == FALSE) {
        showModal(modalDialog(title = "Data error", HTML('The second column of the uploaded file 
        needs to contain numeric values only, no text. <br> Please check and correct 
        your input file.'), size = "l", footer = modalButton("OK")))
      
      } else if (sum(is.na(Results$V2)) > 0) {
        showModal(modalDialog(title = "Data error", HTML('The second column of the uploaded file 
        needs to contain numeric values only, no missing values. <br> Please check and correct 
        your input file.'), size = "l", footer = modalButton("OK")))
     
      } else {
      
      setwd(First.wd)
      Session.ID <- gsub(" ", "", Sys.time(), fixed = TRUE)
      Session.ID <- gsub("-", "", Session.ID, fixed = TRUE)
      Session.ID <- gsub(":", "", Session.ID, fixed = TRUE)
      ## HSA genome UPDATED ++ ##
      # update for human organism code
      #Session.ID <- paste(input$strain, Session.ID, sample(1000:9999, 1), sep = "_")
      Session.ID <- paste("hsa", Session.ID, sample(1000:9999, 1), sep = "_")
      ## HSA genome UPDATED -- ##
      
      Session.path <- paste(First.wd, "/results/", Session.ID, sep ="")
      dir.create(Session.path)
      setwd(Session.path)
      
      file.copy(from = paste(First.wd, "/01_README.pdf", sep = ""), to = Session.path)
      
      # turning table back into list:
      KEGG_List <- tapply(KEGG_Prots, KEGG_Paths, c)
      
      # require at least 4 genes per path:
      KEGG_List_Length <- lapply(KEGG_List, length)
      Good_KEGG <- names(KEGG_List)[KEGG_List_Length > 4]
      
      KEGG_Tests <- lapply(Good_KEGG, function (kname){
        x = KEGG_List[[kname]]
        n = sum(unique(KEGG_Prots) %in% x)
        if (n > 0) {
          o = Results[which(unique(Results$V1) %in% x), "V2"]
          s = sum (o > 0)
          binom.test(x = s, p = 0.5, n = n)    
        } else {binom.test(x = 1, p = 0.5, n = 1)}
        
      })
      names(KEGG_Tests) <- Good_KEGG

      P <- unlist(lapply(KEGG_Tests, function(x){x$p.value}))
      Est <- unlist(lapply(KEGG_Tests, function(x){x$estimate}))
      FDR <- p.adjust(P, "fdr")

      KEGG_FC <- lapply(KEGG_List[Good_KEGG], function(x){
        Results[unique(Results$V1) %in% x, "V2"]
      })

      KEGG_FC_median <- lapply(KEGG_FC, median)

      KEGG_Result <- data.frame(
          "Path" = Good_KEGG,
          "Estimate" = round(as.numeric(Est), 2),
          "Median_fold_change" = round(as.numeric(KEGG_FC_median), 2),
          "P_value" = as.numeric(P),
          "FDR" = as.numeric(FDR))
      
      write.csv(KEGG_Result, "KEGG_Results.csv", row.names = FALSE)
      
      Sig <- KEGG_Result
      Sig <- Sig[Sig$FDR < 0.05, ]
      
      
      KEGGSigList <- KEGG_List[Sig$Path]
      
      KEGGSigFC <- lapply(KEGGSigList, function(x){
              Results[which(unique(Results$V1) %in% x), "V2"]
      })
      
      KEGGSigMedians <- lapply(KEGGSigFC, median)
      
      Cols <- c(1:length(KEGGSigMedians))
      Meds <- sort(unlist(KEGGSigMedians)) < 0
      
      Cols[Meds == TRUE] <- "green"
      Cols[Meds == FALSE] <- "magenta"
      
      if(length(KEGGSigList) >= 1) {
        
        if (length(Cols) == 1){
          
          output$KEGG_plot <- renderPlot({
            par(mar = c(5.1, 21, 4, 2))
            boxplot(lapply(KEGGSigList[names(sort(unlist(KEGGSigMedians)))], function(x){
              Results[which(unique(Results$V1) %in% x), "V2"]}),
              xlab = "log2 FC", horizontal = TRUE, las = 1, pch = ".",
              col = Cols,
              main = "Significant KEGG Pathways")
            abline(v = 0, col = "blue", lty = 3, lwd = 2)
            mtext(names(sort(unlist(KEGGSigMedians))), side = 2, las = 1, adj = 1.5)
          }, height = (3+length(KEGGSigList)*0.15)*100, width = 11*75)
          
          pdf("KEGG_Boxplot.pdf", height=3+length(KEGGSigList)*0.15, width=11)
          par(mar = c(5.1, 21, 4, 2))
          boxplot(lapply(KEGGSigList[names(sort(unlist(KEGGSigMedians)))], function(x){
            Results[which(unique(Results$V1) %in% x), "V2"]}),
            xlab = "log2 FC", horizontal = TRUE, las = 1, pch = ".",
            col = Cols,
            main = "Significant KEGG Pathways")
          abline(v = 0, col = "blue", lty = 3, lwd = 2)
          mtext(names(sort(unlist(KEGGSigMedians))), side = 2, las = 1, adj = 1.5)
          dev.off()
          
        } else {
          
          output$KEGG_plot <- renderPlot({
            par(mar = c(5.1, 21, 4, 2))
            boxplot(lapply(KEGGSigList[names(sort(unlist(KEGGSigMedians)))], function(x){
              Results[which(unique(Results$V1) %in% x), "V2"]}),
              xlab = "log2 FC", horizontal = TRUE, las = 1, pch = ".",
              col = Cols,
              main = "Significant KEGG Pathways")
            abline(v = 0, col = "blue", lty = 3, lwd = 2)
          }, height = (3+length(KEGGSigList)*0.15)*100, width = 11*75)
          
          pdf("KEGG_Boxplot.pdf", height=3+length(KEGGSigList)*0.15, width=11)
          par(mar = c(5.1, 21, 4, 2))
          
          boxplot(lapply(KEGGSigList[names(sort(unlist(KEGGSigMedians)))], function(x){
            Results[which(unique(Results$V1) %in% x), "V2"]}),
            xlab = "log2 FC", horizontal = TRUE, las = 1, pch = ".",
            col = Cols,
            main = "Significant KEGG Pathways")
          
          abline(v = 0, col = "blue", lty = 3, lwd = 2)
          dev.off()
        }
        
      }
      
      
      
      # KEGG Mapper Code
      URL_query_sep <- "?"
      URL_base <- "https://www.kegg.jp/kegg-bin/show_pathway"
      URL_query_args_chars <- "&multi_query="
      URL_multi_query_sep <- "%0d%0a" # strangely this is carriage return, line feed
      UpColor <- "magenta"
      DownColor <- "%2354F338"
      DefaultColor <- "white"
      No_Color ="&nocolor=1"
      
      KeggMapperURL <- function(Path_name){
              ## HSA genome UPDATED ++ ##
              # convert PathGeneTable to dataframe and update to conform to kegg.df structure
              #PathGeneTable <- KEGG[[input$strain]] [KEGG[[input$strain]]$KEGG_paths == Path_name,]
              PathGeneTable<-data.frame(Entry=kegg.df[kegg.df$pathway==Path_name,]$Entry,Cross.reference..KEGG.=kegg.df[kegg.df$pathway==Path_name,]$From,KEGG_paths=kegg.df[kegg.df$pathway==Path_name,]$pathway)
              ## HSA genome UPDATED -- ##
              mergedResults <- merge(PathGeneTable, Results, by.x = "Entry", by.y = "V1" )
              KEGGgenesUp <- mergedResults[mergedResults$V2 > 0, "Cross.reference..KEGG."]
              KEGGgenesDown <- mergedResults[mergedResults$V2 <= 0, "Cross.reference..KEGG."]
              
              URL_multi_queryUp <- paste(paste(KEGGgenesUp, UpColor, sep ='+'),
                                         collapse = URL_multi_query_sep)
              URL_multi_queryDown <- paste(paste(KEGGgenesDown, DownColor, sep ='+'),
                                           collapse = URL_multi_query_sep)
              
              ## HSA genome UPDATED ++ ##
              # new human organism code for KEGG mapper URL
              #myId = paste0(input$strain, Path_List[[input$strain]]$ID[which(Path_List[[input$strain]]$Name == Path_name)])
              myId = paste0("hsa", path.list$ID[which(path.list$Name == Path_name)])
              ## HSA genome UPDATED -- ##
              URL_map_arg <- paste("map", myId, sep = '=')
              paste(URL_base,
                    URL_query_sep,
                    URL_map_arg,
                    URL_query_args_chars,
                    paste(URL_multi_queryUp,URL_multi_queryDown, collapse = URL_multi_query_sep ), 
                    No_Color, sep = '')
              
      }
      
      
      myLinks <- lapply(Sig$Path, KeggMapperURL)
      link_text_to_show <- Sig$Path
      
      contents = paste('<a href="', myLinks,'" target="_blank">', link_text_to_show, '</a>', collapse="</br>")
      
      # Combine contents with other HTML and render using rhtmlMetro
      your.html = paste('<!DOCTYPE html><html><head></head><body>', contents, '</body></html>', sep="")
      
      fileConn<-file("Pathways.html")
      writeLines(your.html, fileConn)
      close(fileConn)
    
      
      if (dim(Sig)[1] == 0) {
        output$KEGG_res <- renderText("No significant KEGG pathways found")
      } else {
      
      output$KEGG_res <- renderText(paste(dim(Sig)[1], "significant KEGG pathways 
                         found out of a total of", length(Good_KEGG), sep = " "))
      
      Links_df <- data.frame("Link" = unlist(myLinks),
                             "Name" = Sig$Path)
      
      Sig$Path <- apply(Links_df, 1, function(x){
              paste('<a href="', x[1],'" target="_blank">', 
                    x[2], '</a>', collapse="</br>")
      })
      
      
      Sig <- Sig[order(Sig$Estimate), ]
      

      output$KEGG_table <- DT::renderDT({
        DT::datatable(Sig, options = list(pageLength = 50), 
                      rownames = FALSE, escape = FALSE, selection="none")
      })
  }
      
      
      # GO term activation analysis
      
      # select subset of GO table based on user input:
      ## HSA genome UPDATED ++ ##
      # update to conform to new go.df structure, no need to subset by strain
      #GO_Prots <- GO[[input$strain]]$Entry[GO[[input$strain]]$Entry %in% Results$V1]
      #GO_Paths <- GO[[input$strain]]$Gene.ontology..GO.[GO[[input$strain]]$Entry %in% Results$V1]
      GO_Prots<-go.df$Entry[go.df$Entry %in% Results$V1]
      GO_Paths<-go.df$go.term[go.df$Entry %in% Results$V1]
      ## HSA genome UPDATED -- ##
      output$GO_match <- renderText(paste(length(unique(GO_Prots)),
                                            "of user-provided identifiers matched to one of",
                                            ## HSA genome UPDATED ++ ##
                                            # update for go.df data structure (no need to subset by strain)
                                            #length(unique(GO[[input$strain]]$Entry)),
                                            length(unique(go.df$Entry)),
                                            ## HSA genome UPDATED -- ##
                                            "possible identifiers for GO term analysis", sep = " "))
      
      # turning table back into list:
      GO_List <- tapply(GO_Prots, GO_Paths, c)
      
      # require at least 4 genes per path:
      GO_List_Length <- lapply(GO_List, length)
      Good_GO <- names(GO_List)[GO_List_Length > 4]
      
      GO_Tests <- lapply(Good_GO, function (kname){
        x = GO_List[[kname]]
        n = sum(unique(GO_Prots) %in% x)
        if (n > 0) {
          o = Results[which(unique(Results$V1) %in% x), "V2"]
          s = sum (o > 0)
          binom.test(x = s, p = 0.5, n = n)    
        } else {binom.test(x = 1, p = 0.5, n = 1)}
      })
      names(GO_Tests) <- Good_GO
      
      GO_Tests_P <- unlist(lapply(GO_Tests, function(x){x$p.value}))

      GO_Tests_Est <- unlist(lapply(GO_Tests, function(x){x$estimate}))
      names(GO_Tests_Est) <- names(GO_Tests_P)

      GO_FDR <- p.adjust(GO_Tests_P, "fdr")

      GO_FC <- lapply(GO_List[Good_GO], function(x){
        Results[which(unique(Results$V1) %in% x), "V2"]
      })

      GO_FC_median <- lapply(GO_FC, median)


      GO_Result <- data.frame(
        "GO_term" = names(GO_Tests),
        "Estimate" = round(as.numeric(GO_Tests_Est), 2),
        "Median_fold_change" = round(as.numeric(GO_FC_median), 2),
        "P_value" = as.numeric(GO_Tests_P),
        "FDR" = as.numeric(GO_FDR))

      write.csv(GO_Result, "GO_Results.csv", row.names = FALSE)


      GO_Sig <- GO_Result[GO_Result$FDR < 0.05, ]
      
      if (dim(GO_Sig)[1] == 0) {
        output$GO_res <- renderText("No significant GO terms found")
      } else {
        
        output$GO_res <- renderText(paste(dim(GO_Sig)[1], "significant GO terms
                         found out of a total of", length(Good_GO), sep = " "))
      
      GO_Sig <- GO_Sig[order(GO_Sig$Estimate), ]

      output$GO_table <- DT::renderDT({
        DT::datatable(GO_Sig, options = list(pageLength = 50), rownames = FALSE, selection="none")
      })

      GO_SigList <- GO_List[GO_Sig$GO_term]

      GO_SigFC <- lapply(GO_SigList, function(x){
        Results[which(unique(Results$V1) %in% x), "V2"]
      })

      GO_SigMedians <- lapply(GO_SigFC, median)

      GoCols <- c(1:length(GO_SigMedians))
      GoMeds <- sort(unlist(GO_SigMedians)) < 0

      GoCols[GoMeds == TRUE] <- "green"
      GoCols[GoMeds == FALSE] <- "magenta"
      
      if (length(GoCols) == 1){
        
      output$GO_plot <- renderPlot({
        par(mar = c(5.1, 25, 4, 2))
        boxplot(lapply(GO_SigList[names(sort(unlist(GO_SigMedians)))], function(x){
          Results[which(unique(Results$V1) %in% x), "V2"]}),
          col = GoCols,
          xlab = "log2 FC", horizontal = TRUE, las = 1,
          main = "Significant GO terms")
        abline(v = 0, col = "red", lty = 3, lwd = 2)
        mtext(names(sort(unlist(GO_SigMedians))), side = 2, las = 1, adj = 1.5)
      }, height = (3+length(GO_SigList)*0.15)*100, width = 11*75)

      pdf("GO_Boxplot.pdf", height=3+length(GO_SigList)*0.15, width=11)
      par(mar = c(5.1, 25, 4, 2))
      boxplot(lapply(GO_SigList[names(sort(unlist(GO_SigMedians)))], function(x){
        Results[which(unique(Results$V1) %in% x), "V2"]}),
        col = GoCols,
        xlab = "log2 FC", horizontal = TRUE, las = 1,
        main = "Significant GO terms")
      abline(v = 0, col = "blue", lty = 3, lwd = 2)
      mtext(names(sort(unlist(GO_SigMedians))), side = 2, las = 1, adj = 1.5)
      dev.off()
      
      } else {
      
        output$GO_plot <- renderPlot({
          par(mar = c(5.1, 25, 4, 2))
          boxplot(lapply(GO_SigList[names(sort(unlist(GO_SigMedians)))], function(x){
            Results[which(unique(Results$V1) %in% x), "V2"]}),
            col = GoCols,
            xlab = "log2 FC", horizontal = TRUE, las = 1,
            main = "Significant GO terms")
          abline(v = 0, col = "blue", lty = 3, lwd = 2)
        }, height = (3+length(GO_SigList)*0.15)*100, width = 11*75)
        
        pdf("GOBoxplot.pdf", height=3+length(GO_SigList)*0.15, width=11)
        par(mar = c(5.1, 25, 4, 2))
        boxplot(lapply(GO_SigList[names(sort(unlist(GO_SigMedians)))], function(x){
          Results[which(unique(Results$V1) %in% x), "V2"]}),
          col = GoCols,
          xlab = "log2 FC", horizontal = TRUE, las = 1,
          main = "Significant GO terms")
        abline(v = 0, col = "blue", lty = 3, lwd = 2)
        dev.off()
      }}

      if (length(list.files(pattern = '\\.csv$')) > 0) {
        ## HSA genome UPDATED ++ ##
        # remove strain from output name
        #Composite.Name <- paste('Results_', input$strain, '_', Sys.Date(), '.zip', sep='')
        Composite.Name <- paste('Results_', Sys.Date(), '.zip', sep='')
        ## HSA genome UPDATED -- ##
        SystemCall <- paste("find . \\( -name '*.pdf' -or -name '*.csv' -or -name '*.html' \\) -print | zip ", 
                            Composite.Name, " -@", sep ="")
        system(SystemCall)}
      
      # Log app usage for Date-Time; Organism; # of sig KEGG; # of sig GO:
      ## HSA genome UPDATED ++ ##
      # remove strain from output name
      #write.table(paste(Sys.time(), input$strain, dim(Sig)[1], dim(GO_Sig)[1], sep = ";"),
      write.table(paste(Sys.time(),dim(Sig)[1], dim(GO_Sig)[1], sep = ";"),            
      ## HSA genome UPDATED -- ##
                  file = paste(First.wd, "/ESKAPE_usage_log.txt", sep =""),  
                  append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

      output$ui <- renderUI({
        tagList(
          verbatimTextOutput('done'),
          downloadButton('download', 'Download results'),
          
          hr(),
          actionButton('reset', "Reset", icon("refresh", class = "fa-spin"),
                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
        )
        
      })
      
      
      output$download <- downloadHandler(
        filename = Composite.Name,
        content <- function(file) {
          file.copy(from = paste(Session.path, "/", Composite.Name, sep = ""), to = file)
        },
        contentType = "application/zip")
    
  }
  
      output$action <- renderUI({})  
      
      setwd(First.wd)
      
  })

  
  
  # Reset
  observeEvent(input$reset, {

      shinyjs::reset('file1')
    ## HSA genome UPDATED ++ ##
    # remove strain reset
      #shinyjs::reset('strain')
    ## HSA genome UPDATED -- ##
   
      updateTabsetPanel(session, 'Results_Output', selected = 'Summary')
      
      output$error <- renderUI({})
      output$action <- renderUI({})
      output$ui <- renderUI({})
      
      output$KEGG_match <- renderText({})
      output$KEGG_res <- renderText({})
      output$GO_match <- renderText({})
      output$GO_res <- renderText({})
      
      output$KEGG_plot <- renderPlot({})
      output$KEGG_table <- DT::renderDT({})
      output$GO_plot <- renderPlot({})
      output$GO_table <- DT::renderDT({})
      
      setwd(First.wd)

    })

}

shinyApp(ui = ui, server = server)
