library(shiny)
library(shinythemes)
library(data.table)
library(DT)
# source("library/PPARGC1A_B_CONSTS.R")
.GENE_IDs <- list("PPARGC1A", "PPARGC1B")

GENE_NAMES <- c ("PPARGC1A_1", "PPARGC1A_2", "PPARGC1A_3", "PPARGC1A_4", "PPARGC1A_5", "PPARGC1A_6", "PPARGC1B_1")
GENE_IDS <- append(c((rep(.GENE_IDs[1], 6))), c(.GENE_IDs[2]))
GENE_VERSIONS <- array(1:6, dim = c(7))

PLOT_GENE_ENDPOINT <- "https://epd.epfl.ch/cgi-bin/plot_gene.php?"
GET_SEQUENCE_ENDPOINT <- "https://epd.epfl.ch/cgi-bin/get_sequence.php?"

# MOTIF_LENGTH <- 13
MOTIF_LENGTH <- 8

TF <- "MA0105.4"
TF_NAME <- "NFKB1"
OFFSET <- -5

ROOT_PATH <- "C:/Users/Erik.Butcher/Google Drive (butcer@alumni.wfu.edu)/Engineering/R/Projects/Testing_R_112619/Testing_R"



# source("library/functions.R")
library(tidyverse) ; library(httr) ; library(jsonlite) ; library(rvest) ; library(xlsx) ; library(dplyr)

doSomething <- function(message) {
   print(paste("doing something:", message));
}

getContext <- function(request) {
   return(content(request, as = "text", encoding = "UTF-8"))
}

parseSequenceLocations <- function(html_output) {
   p_content <- c(html_output %>% html_nodes("p:not(div)") %>% html_text(trim=TRUE))
   
   p_content <- str_replace_all(p_content, "::", "|");
   print(paste("content:", p_content ));
   p_pieces_2 <- toString(as.list(strsplit(p_content, ":")[[1]])[2]);
   
   # print(paste("pieces", p_pieces_2));
   
   return(as.list(strsplit(p_pieces_2, ",")[[1]]))
}

getSequence <- function(sequenceLocs, motifLength, transFactor = "MA1117.1", geneId = "PPARGC1B_1", offset = OFFSET) {
   request <- GET(url = GET_SEQUENCE_ENDPOINT,
                  query = list(
                     database = "epdNew_hg",
                     gene_id = geneId,
                     tf = transFactor,
                     from = as.integer(sequenceLocs[1]) + offset,
                     to = (as.integer(sequenceLocs[1]) + motifLength + offset))
   )
   
   html_output <- read_html(request, encoding="UTF-8")
   # html_output %>% xml_structure()
   response_text <- html_output %>% html_node("p") %>% html_text
   response_pieces <- as.list(strsplit(response_text, "\n")[[1]])
   
   return(response_pieces[2])
   
   # p_pieces_2 <- toString(as.list(strsplit(p_content, ":")[[1]])[2])
   # p_consensus_regions <- as.list(strsplit(p_pieces_2, ",")[[1]]) %>% print
}

getSequences <- function(consensus_regions_list, tf, geneId, motif_length) {
   # print(paste("Gene:",GENE_IDS[1],"-","Trans Factor:", TF_NAME))
   for (consensus_regions in consensus_regions_list) {
      for (location in consensus_region) {
         sequence <- getSequence(location, motif_length, tf, geneId)
         print(paste("    location:", location, "-", sequence));
      }
   }
}

getLocationsAndSequences <- function(gene_ids, motifs, from, to, p, offset, motif_length) {
   print("getting genes and locations")
   
   
   result_sequences <- data.frame(NULL)
   for (gene_id in gene_ids ) {
      for (row in 1:nrow(motifs)) {
         transFactor <-  motifs[["motifIds"]][row]
         print(paste("gene:",gene_id, "motif_id", motifs[["motifIds"]][row], "motif_names", motifs[["motifNames"]][row]))
         request <- GET(url = PLOT_GENE_ENDPOINT, 
                        query = list(
                           database = "epdNew_hg",
                           gene_id = gene_id,
                           tf = transFactor,
                           from = toString(from),
                           to = toString(to),
                           motif_lib = "jasparCoreVertebrates",
                           tfname = motifs[["motifNames"]][row],
                           co=p))
         
         print(paste("Request:", request))
         html_output <- read_html(request, encoding="UTF-8")
         print(paste("Response:", html_output))
         p_consensus_regions <- parseSequenceLocations(html_output)
         # print(p_consensus_regions)
         # browser()
         # result_sequences <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("site", "sequence"))
         sequences <- c()
         for (location in p_consensus_regions) {
            sequence <- getSequence(location, motif_length, transFactor, gene_id, offset)
            seq_set <- data.frame(as.list(location, sequence))
            sequences <- append(sequences, sequence)
            # print(paste("    location:", location, "-", sequence));
         }
         
         glimpse(p_consensus_regions)
         glimpse(sequences)
         
         current_sequences <- data.frame(gene= rep(gene_id, times = length(sequences)), motif = rep(motifs[["motifNames"]][row], times = length(sequences)), site = I(p_consensus_regions), sequence = I(sequences))
         result_sequences <- rbind(result_sequences, current_sequences)
         print("complete")
         str(result_sequences)
      }
   }
   return (result_sequences)
}

getGeneMenu <- function() {
   endpoint <- "https://epd.epfl.ch/searchMotifBuildMenu.php?lib=jasparCoreVertebrates"
   
   request <- GET(url = endpoint)
   response <- content(request, as = "text", encoding = "UTF-8")
   geneMenuData <- as.list(strsplit(response, "\n")[[1]])
   
   geneIds <- c()
   geneNames <- c()
   #   geneMenu <- vector("list", length = length(geneMenuData))
   for (i in seq_along(geneMenuData)) {
      geneIds <- append(geneIds, (c(strsplit(toString(geneMenuData[i]), "[[:space:]]")[[1]])[1]))
      geneNames <- append(geneNames, (c(strsplit(toString(geneMenuData[i]), "[[:space:]]")[[1]])[2]))
   }
   geneMenu <- data.frame(geneIds, geneNames)
   return(geneMenu)
}

getMotifMenu <- function() {
   endpoint <- "https://epd.epfl.ch/searchMotifBuildMenu.php?lib=jasparCoreVertebrates"
   
   request <- GET(url = endpoint)
   response <- content(request, as = "text", encoding = "UTF-8")
   motifMenuData <- as.list(strsplit(response, "\n")[[1]])
   
   motifIds <- c()
   motifNames <- c()
   
   for (i in seq_along(motifMenuData)) {
      motifIds <- append(motifIds, (c(strsplit(toString(motifMenuData[i]), "[[:space:]]")[[1]])[1]))
      motifNames <- append(motifNames, (c(strsplit(toString(motifMenuData[i]), "[[:space:]]")[[1]])[2]))
   }
   motifMenu <- data.frame(motifIds, motifNames)
   # remove any `is.NA`
   motifMenu <- motifMenu[complete.cases(motifMenu),]
   return(motifMenu)
}

doSaveToExcel <- function(data, path) {
   res <- write.xlsx(data, path)
}



motifMenu <- getMotifMenu()

motif_ids <- motifMenu[["motifIds"]]
motif_names <- motifMenu[["motifNames"]]

ui <- fluidPage( 
   # theme = shinytheme("darkly"),
   # theme = shinytheme("cybord"),
   # themeSelector(),
   title = 'Sequence One Alpha',
   sidebarLayout(
      sidebarPanel(
         div(style="display: inline-block;vertical-align:top; margin: 0.5em",
             
             div(style="display: inline-block;vertical-align:top; width: 100%;",selectizeInput(
                'gene_names_input', '1. Select Genes', choices = GENE_NAMES, multiple = TRUE
             )),
             div(style="display: inline-block;vertical-align:top; width: 100%", selectizeInput(
                'motif_names_input', '2. Select Motifs', choices = motif_names, multiple = TRUE
             )),
             div(style="display: inline-block;horizontal-align:center; width: 100%",
                 div(style="display: inline-block;vertical-align:top; width: 32.3%; min-width: 80px",numericInput("from_input", "From:", -1000, max = 0, min = -2000)),
                 div(style="display: inline-block;vertical-align:top; width: 32.3%; min-width: 80px",numericInput("to_input", "To:", 1000, max = 2000, min = 0)),
                 div(style="display: inline-block;vertical-align:top; width: 32.3%; min-width: 80px", selectInput("p_input", "P_Value", c(0.05, 0.01, 0.001, 0.0001), selected=0.001))
             ),
             div(style="display: inline-block;vertical-align:center; width: 48%;", sliderInput("offset_input", "Offset",
                                                                                               min = -12, max = 0,
                                                                                               value = OFFSET)),
             div(style="display: inline-block;vertical-align:center; width: 48%;", sliderInput("motif_length_input", "motif length",
                                                                                               min = 5, max = 16,
                                                                                               value = MOTIF_LENGTH)),
             div(style="display: inline-block;vertical-align:center; width: 98%;", textInput("save_path_input", "save path", paste(ROOT_PATH, "/output.xlsx", sep="")),
               tags$head(tags$style("#save_path_input{
                                 font-size: 12px;
                                    }"
                   )
                )),
             div(style="display: inline-block;vertical-align:center; width: auto; margin-right: 10px", actionButton("save_excel_button", "Save to Excel", icon("list-alt"))),
             div(style="display: inline-block;vertical-align:center; width: 48%;", checkboxInput("showStyle", "Show Style Selector", FALSE)),
         )),
      mainPanel(
         # verbatimTextOutput('selection')
         conditionalPanel(
            condition = "input.showStyle",
            themeSelector()
         ),
         fluidRow(
            column(12, DTOutput("sequence_table"))
         ),
         fluidRow(
            helpText('Please continue updating values while the data is loading. The data view will continually update to reflect any changes.')
         ),
      )
   )
)

server <- function(input, output, session) {
   sequenceDT <- data.frame()
   formData <- reactive({
      data <- sapply(c("sequence_table"), function(x) input[[x]])
      data
   })
   
   output$sequence_table <- renderDT({
      sequenceDT <- data.frame()
      sequenceDT <- iris
      # 
      # return(sequenceDT)
      
      motifs <- subset(motifMenu, motifNames %in% input$motif_names_input)
      # print(iris)
      # class(iris)
      # sequenceDT <- data.frame()
      # TODO: Add loading context here
      # names(sequenceDT) <- c("motifIds", "motifNames")
      if (length(motifs[[1]]) > 0 && length(input$gene_names_input) > 0) {
         sequenceDataTable <- getLocationsAndSequences(input$gene_names_input, motifs, input$from_input, input$to_input, input$p_input, input$offset_input, input$motif_length_input)
         sequenceDT <- sequenceDataTable
         return(sequenceDataTable)
      } else {
         return(data.frame(NULL))
      }
   })
   
 
   
   loadData <- function() {
      if (exists("responses")) {
         responses
      }
   }
   
   observeEvent(input$save_excel_button, {
      output <- data.frame()
      # output <- iris
            # sequenceDT <- data.frame()
      # sequenceDT <- iris
      # 
      # return(sequenceDT)
      
      motifs <- subset(motifMenu, motifNames %in% input$motif_names_input)
      # print(iris)
      # class(iris)
      # sequenceDT <- data.frame()
      # TODO: Add loading context here
      # names(sequenceDT) <- c("motifIds", "motifNames")
      if (length(motifs[[1]]) > 0 && length(input$gene_names_input) > 0) {
         sequenceDataTable <- getLocationsAndSequences(input$gene_names_input, motifs, input$from_input, input$to_input, input$p_input, input$offset_input, input$motif_length_input)
         output <- sequenceDataTable
         # return(sequenceDataTable)
      } 
      # else {
      #    return(data.frame(NULL))
      # }
      
      
      # print("saving excel")
      # print(iris)
      save.value <- input$save_excel_button
      
      # todo: make path dynamic
      # save.path <- paste(ROOT_PATH, "/iris_test.xlsx", sep="")
      save.path <- input$save_path_input
      print(paste("Path:", save.path))
      # data(sequenceDT)
      # glimpse(sequenceDT)
      # as.data.frame(t(
      doSaveToExcel(output, save.path)
      # doSaveToExcel(formData, save.path)
      # doSaveToExcel(sequenceDT, "/my_test_sequence_data.xlsx")
      return("Save Successful")
   })
}

shinyApp(ui = ui, server = server)
# app <- shinyApp(ui = ui, server = server)
# runApp(app)
