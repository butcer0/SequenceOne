#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

#CONSTS
.GENE_IDs <- list("PPARGC1A", "PPARGC1B")

GENE_NAMES <- c ("PPARGC1A_1", "PPARGC1A_2", "PPARGC1A_3", "PPARGC1A_4", "PPARGC1A_5", "PPARGC1A_6", "PPARGC1B_1")
GENE_IDS <- append(c((rep(.GENE_IDs[1], 6))), c(.GENE_IDs[2]))
GENE_VERSIONS <- array(1:6, dim = c(7))

PLOT_GENE_ENDPOINT <- "https://epd.epfl.ch/cgi-bin/plot_gene.php?"
GET_SEQUENCE_ENDPOINT <- "https://epd.epfl.ch/cgi-bin/get_sequence.php?"

MOTIF_LENGTH <- 13

TF <- "MA0105.4"
TF_NAME <- "NFKB1"
OFFSET <- -4

#Helper Functions
library(tidyverse) ; library(httr) ; library(jsonlite) ; library(rvest) ; library(xlsx)

doSomething <- function(message) {
   print(paste("doing something:", message));
}

getContext <- function(request) {
   return(content(request, as = "text", encoding = "UTF-8"))
}

parseSequenceLocations <- function(html_output) {
   p_content <- html_output %>% html_nodes("p:not(div)") %>% html_text(trim=TRUE)
   p_pieces_2 <- toString(as.list(strsplit(p_content, ":")[[1]])[2])
   return(as.list(strsplit(p_pieces_2, ",")[[1]]))
}

getSequence <- function(sequenceLocs, motifLength) {
   request <- GET(url = GET_SEQUENCE_ENDPOINT,
                  query = list(
                     database = "epdNew_hg",
                     gene_id = "PPARGC1B_1",
                     tf = "MA0105.4",
                     from = as.integer(sequenceLocs[1]) + OFFSET,
                     to = (as.integer(sequenceLocs[1]) + motifLength + OFFSET))
   )
   
   html_output <- read_html(request, encoding="UTF-8")
   # html_output %>% xml_structure()
   response_text <- html_output %>% html_node("p") %>% html_text
   response_pieces <- as.list(strsplit(response_text, "\n")[[1]])
   
   return(response_pieces[2])
   
   # p_pieces_2 <- toString(as.list(strsplit(p_content, ":")[[1]])[2])
   # p_consensus_regions <- as.list(strsplit(p_pieces_2, ",")[[1]]) %>% print
}

getSequences <- function(consensus_regions_list) {
   # print(paste("Gene:",GENE_IDS[1],"-","Trans Factor:", TF_NAME))
   for (consensus_regions in consensus_regions_list) {
      for (location in consensus_region) {
         sequence <- getSequence(location, MOTIF_LENGTH)
         print(paste("    location:", location, "-", sequence));
      }
   }
}

getLocationsAndSequences <- function(gene_ids, motifs, from, to, p) {
   print("getting genes and locations")
   
   
   result_sequences <- data.frame(NULL)
   for (gene_id in gene_ids ) {
      for (row in 1:nrow(motifs)) {
         print(paste("gene:",gene_id, "motif_id", motifs[["motifIds"]][row], "motif_names", motifs[["motifNames"]][row]))
         request <- GET(url = PLOT_GENE_ENDPOINT, 
                        query = list(
                           database = "epdNew_hg",
                           gene_id = gene_id,
                           tf = motifs[["motifIds"]][row],
                           from = toString(from),
                           to = toString(to),
                           motif_lib = "jasparCoreVertebrates",
                           tfname = motifs[["motifNames"]][row],
                           co=p))
         
         print(paste("Request:", request))
         html_output <- read_html(request, encoding="UTF-8")
         p_consensus_regions <- parseSequenceLocations(html_output)
         # print(p_consensus_regions)
         # browser()
         # result_sequences <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("site", "sequence"))
         sequences <- c()
         for (location in p_consensus_regions) {
            sequence <- getSequence(location, MOTIF_LENGTH)
            seq_set <- data.frame(as.list(location, sequence))
            sequences <- append(sequences, sequence)
            # print(paste("    location:", location, "-", sequence));
         }
         
         print(p_consensus_regions)
         print(sequences)
         
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
   saveToExcel(data, path)
}

# Define UI for application that draws a histogram
motifMenu <- getMotifMenu()

motif_ids <- motifMenu[["motifIds"]]
motif_names <- motifMenu[["motifNames"]]

ui <- fluidPage( 
   theme = shinytheme("darkly"),
   title = 'Gene Sequence Retrieval Tool',
   sidebarLayout(
      sidebarPanel(
         div(style="display: inline-block;vertical-align:top; width: 250px;",selectizeInput(
            'gene_names_input', '1. Select Genes', choices = GENE_NAMES, multiple = TRUE
         )),
         div(style="display: inline-block;vertical-align:top; width: 250px;", selectizeInput(
            'motif_names_input', '2. Select Motifs', choices = motif_names, multiple = TRUE
         )),
         div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("from_input", "From:", -1000, max = 0, min = -2000)),
         div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("to_input", "To:", 1000, max = 2000, min = 0)),
         div(style="display: inline-block;vertical-align:top; width: 150px;", selectInput("p_input", "P_Value",
                                                                                          c(0.05, 0.01, 0.001, 0.0001), selected=0.001)),
      ),
      mainPanel(
         # verbatimTextOutput('selection')
         fluidRow(
            column(12, dataTableOutput("sequence_table"))
         ),
         fluidRow(
            helpText('Please continue updating values while the data is loading. The data view will continually update to reflect any changes.')
         ),
      )
   )
)

server <- function(input, output) {
   output$sequence_table <-   renderDataTable({
      motifs <- subset(motifMenu, motifNames %in% input$motif_names_input)
      # print(iris)
      # class(iris)
      sequenceDT <- data.frame()
      # TODO: Add loading context here
      # names(sequenceDT) <- c("motifIds", "motifNames")
      if (length(motifs[[1]]) > 0 && length(input$gene_names_input) > 0) {
         sequenceDataTable <- getLocationsAndSequences(input$gene_names_input, motifs, input$from_input, input$to_input, input$p_input)
         # print(sequenceDataTable)
         return(sequenceDataTable)
      } else {
         return(invisible(NA))
      }
      
      # print("data table generated")
      # str(sequenceDT)
      # paste('You selected', if (length(input$motif_names_input) == 0) 'nothing' else toString(motifs),
      #       'in the Gene Sequences example.')
      # 
      # if(length (sequenceDT[[1]]) > 0) {
      #   return(sequenceDT)
      # } else {
      #   return(data.table())
      # }
      
      # return(mtcars)
      
      # return(sequenceDT)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

