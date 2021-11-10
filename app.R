library(edgeR) # Used for CPM calculation
library(reshape2) # count reformatting
library(plyr)
library(ggpubr)
library(DT)
library(shinyWidgets)
library(dplyr)

# Script for timeout feature
timeoutSeconds <- 60*15
inactivity <- sprintf("function idleTimer() {
var t = setTimeout(logout, %s);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions
function logout() {
Shiny.setInputValue('timeOut', '%ss')
}
function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
}
}
idleTimer();", timeoutSeconds*1000, timeoutSeconds, timeoutSeconds*1000)

ui <- fluidPage(
  
  shinyUI(fluidPage(tags$script(inactivity))),
  
  tabsetPanel(  
    tabPanel("Processor", fluid = TRUE,
              
             titlePanel("Functional Screen"),
             # Copy the line below to make a text input box
             sidebarPanel(

               fileInput("pcounts", "Upload counts",
                         multiple = TRUE,
                         accept = c('.txt')),
               
               fileInput("samplelookup", "Upload samples",
                         multiple = FALSE,
                         accept=c('.txt')
                         )
               
               # fileInput("comp", "Upload comparison",
               #           multiple = FALSE,
               #           accept=c('.txt')
               # )
             ),
             
             mainPanel(
               h3("Raw Counts"),
               fluidRow(
                 column(12,dataTableOutput('table')
                 )
               )
               
             ) # End tab one main
             
    ), # end tab one
    
    tabPanel('Plot', fluid=TRUE,
            sidebarPanel(
              radioButtons("norm", "Scale",
                           c("Log2 (Raw Counts)" = "Log2 (Raw Counts)",
                             "Log2 (CPM)" = "Log2 (CPM)")),
              
              selectizeInput('gene_selection',
                             choices = NULL,
                             multiple = F,
                             label = HTML('<p style=“color:black;“>Select Gene</p>')),
              pickerInput("group_selection","Select Groups", choices=NULL, options = list(`actions-box` = TRUE),multiple = T),
              downloadButton("downloadPlot", "Download Plot")
              
            ),
             mainPanel(
               plotOutput('plot')
             )
    
    
  )
  
  
  
  

))


server <- function(input, output, session) {

  
  # Counts table
  counts <- reactive({
    
    req(input$pcounts)
    
    path_counts <- input$pcounts$datapath
    # ID should be a combination of Gene and sequence/ID. Indicate which of a "_" split is the gene
    n_order <- 1 # In this example mine is A1BG_CATCTTCTTTCACCTGAACG, so 1 is the gene
    
    # Output directory name
    output_dir <- "Analysis/sgRNA_plots"
    
    
    dir.create(output_dir,
               recursive = T)
    
    temp_counts <- read.table(file = path_counts,
                              sep = "\t", header = T, stringsAsFactors = F)
    counts <- cbind.data.frame(ID = temp_counts$ID,
                               gene = sapply(strsplit(temp_counts$ID, split = "_"),
                                             function(x){x[[n_order]]}),
                               temp_counts[,2:ncol(temp_counts)])
    counts[is.na(counts)] <- 0
    
    
    
    return(counts)
  })
  
  sample_lookup <- reactive({
    req(input$samplelookup)
    path_sample_lookup <- input$samplelookup$datapath
    
    sample_table <- read.table(file =  path_sample_lookup,
                                sep = "\t", header = T, stringsAsFactors = F)
    
    return(sample_table)
  })
  
  selected_groups <- reactive({
    req(sample_lookup())
    
    group_options <- unique(sample_lookup()$Group)
    selection_groups <- group_options # select all groups
    
    return(selection_groups)
  })
  
  sample_lookup_sub <- reactive({
    req(selected_groups())
    req(sample_lookup())
    req(input$norm)
    
    sample_lookup_s <- sample_lookup()[sample_lookup()$Group %in% selected_groups(),]
    
    #lock in order
    sample_lookup_s$Sample <- factor(sample_lookup_s$Sample,
                                       levels = sample_lookup_s$Sample)
    return(
      sample_lookup_s
    )
    
  })
  
  normalization_vals <- reactive({
    req(selected_groups())
    req(sample_lookup_sub())
    req(input$norm)
    req(counts())
  
    counts_ <- data.frame(counts())

    counts_sub <- counts_[,match(
      c("ID", "gene", as.character(sample_lookup_sub()$Sample)),
      colnames(counts_)
    )]
    
    if(input$norm == 'Log2 (CPM)') {
      norm_sub <- cbind(counts_sub[,1:2],
                       log2(cpm(counts_sub[,3:ncol(counts_sub)])+1))
    }
    
    else {
      norm_sub <- cbind(counts_sub[,1:2],
                              log2(counts_sub[,3:ncol(counts_sub)]+1))
    }
    
    
    return(norm_sub)
  })
  
  plotter <- reactive({
    
    req(input$gene_selection)
    req(normalization_vals())
    req(selected_groups())
    req(input$group_selection)
    req(counts())
    
    gene <- input$gene_selection
    
    df <- normalization_vals()
    df_type <- names(normalization_vals())
    
    df_type_saving <- gsub(" ", "_", df_type) # Adjusting for saving file name
    df_type_saving <- gsub("\\(|\\)", "", df_type_saving)
    
    df <- df[df$gene %in% gene, -2]
    
    df <- melt(df)
    colnames(df) <- c("ID",
                      "Sample",
                      "Value")
    
    sample_lookup_sub2 <- sample_lookup_sub()[,c("Sample",
                                              "Group")] # In case there are extra columns
    
    df <- join(df,
               sample_lookup_sub2)
    
    df$Sample <- factor(df$Sample,
                        levels = levels(sample_lookup_sub()$Sample))
    
    #print(df_type)
    
    #df <- dplyr::filter(df, Group == input$group_selection)
    
    df <- df[df$Group %in% input$group_selection, ]
    
    
    #print(df)
    
    
    ggl_facet <- df %>% 
      ggline("Sample", "Value",
             color = "ID",
             legend = "right", legend.title = "sgRNA",
             title = gene,
             font.legend = c(10, 
                             #"bold", 
                             "black"),
             ylab = input$norm) + ylim(0, NA) +
               facet_grid(~Group, scales = "free_x", # Let the x axis vary across facets.
                                   space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                                   switch = "x")     # Move the facet labels to the bottom.)
    
    
    plot(ggl_facet)
    
    
    
  })
  
  output$plot <- renderPlot({
    plotter()
  })
  
  
  # Outputted counts table
  output$table <- renderDT(counts(), filter = 'top', options = list(
    pageLength = 10
  ))
  
  # Gene selection
  gene_options <- reactive({
    genes <- unique(counts()$gene)
  }) 
  
  
  # Update the gene selections for the user
  observeEvent(gene_options(), {
    
    updateSelectizeInput(session, 'gene_selection',
                         choices = as.character(gene_options()),
                         server = TRUE,
                         selected = NULL)
  }, once = FALSE)
  
  observeEvent(selected_groups(), {
    
    updatePickerInput(session, 'group_selection',
                         choices = as.character(selected_groups()))
  }, once = FALSE)
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("FS-", input$gene_selection, ".pdf", sep="")
    },
    content = function(file) {
      ggsave(file, plotter(), device = "pdf", width=11, height=8.5)    
    }
  )

  # Comparison table 
  
  mageck_results <- reactive({
    mgk <- input$comp$datapath
    mgk_table <- read.table(file = mgk,
                              sep = "\t", header = T, stringsAsFactors = F)

    return(mgk_table)
  })
  
  observeEvent(input$timeOut, { 
    # Modified from: https://stackoverflow.com/questions/33839543/shiny-server-session-time-out-doesnt-work
    showModal(modalDialog(
      title = "Inactivity Timeout",
      paste("Session timeout due to", input$timeOut, "inactivity -", Sys.time()),
      footer = NULL
    ))
    stopApp() # I made this change so that the app closes instead of the window only
  })
  

  
}

shinyApp(ui, server)