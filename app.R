library(edgeR) # Used for CPM calculation
library(reshape2) # count reformatting
library(plyr)
library(ggpubr)
library(DT)
library(shinyWidgets)
library(dplyr)
library(ComplexHeatmap)
library(circlize) # Used for colorRamp2 function
library(randomcoloR)

# Script for timeout feature
# timeoutSeconds <- 60*15
# inactivity <- sprintf("function idleTimer() {
# var t = setTimeout(logout, %s);
# window.onmousemove = resetTimer; // catches mouse movements
# window.onmousedown = resetTimer; // catches mouse movements
# window.onclick = resetTimer;     // catches mouse clicks
# window.onscroll = resetTimer;    // catches scrolling
# window.onkeypress = resetTimer;  //catches keyboard actions
# function logout() {
# Shiny.setInputValue('timeOut', '%ss')
# }
# function resetTimer() {
# clearTimeout(t);
# t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
# }
# }
# idleTimer();", timeoutSeconds*1000, timeoutSeconds, timeoutSeconds*1000)

ui <- fluidPage(
  
  #shinyUI(fluidPage(tags$script(inactivity))),
  
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
    
    
  ), # end tab 2
  
  tabPanel('Heatmap', fluid=TRUE,
           sidebarPanel(
             fileInput("HM_genes", "Upload gene list",
                       multiple = FALSE,
                       accept = c("txt",
                                  ".csv")),
             radioButtons("summarize_to_gene", "Summarize to Gene",
                          c("True" = "T",
                            "False" = 'F')),
             tags$style(".btn-file {  
                                         background-color:#878787; 
                                         border-color:#878787; 
                                         }
                                         
                                         .progress-bar {
                                         background-color:#01336d;
                                         }"), # Changes color related to file input, #082959
             downloadButton("Heatmap_example_download", "Download example gene list",
                            class = "butt"),
             br(),
             uiOutput("HM_col_ann_ui"),
             uiOutput("HM_row_ann_ui"),
             uiOutput("HM_row_ann"),
             checkboxGroupInput("HM_options", label = "Select heatmap options to perform:", 
                                choices = list("Z-score transform" = "Z", "Cluster samples" = "col_cluster",
                                               "Cluster genes" = "row_cluster", "Show sample names" = "colnames", 
                                               "Show gene names" = "rownames")
             ),
             actionButton("Heatmap_button", "Generate Plot",
                           styleclass = "success")
           ), 
           
           mainPanel(
             plotOutput("Heatmap_view",  height = "600px")
           )
             
           )
          
           
           
  )
  
  
  
  

)


server <- function(input, output, session) {
  
  # HM_genes_df <- reactive({
  #   req(input$HM_genes)
  #   
    #HM_genes_df <- read.table(file = input$HM_genes$datapath,
  #                            sep = "\t", header = T, stringsAsFactors = F)
  #   
  #   return(HM_genes_df)
  # })
  

  
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
    
    
    #print(counts)
    return(counts)
  })
  
  
  cpm_counts <- reactive({
   # req(input$pcounts)
    req(counts())
    
    cpm_ <- cbind(counts()[,1:2],
                  cpm(counts()[,3:ncol(counts())]))
    
    return(cpm_)
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
    ## all groups
    
    
    #lock in order
    sample_lookup_s$Sample <- factor(sample_lookup_s$Sample,
                                       levels = sample_lookup_s$Sample)
    return(
      sample_lookup_s
    )
    
  })
  
  HM_genes_df <- reactive({
    req(input$HM_genes)
    
    gene_list <- input$HM_genes$datapath
    
    df <- read.table(file = gene_list,
                     sep = "\t", header = T, stringsAsFactors = F)
    
    df <- data.frame(df)
    return(df)
    
  })
  
  
  CPM_table <- reactive({
    req(sample_lookup_sub())
    req(HM_genes_df())
    req(sample_lookup_sub())
    req(cpm_counts())
    req(input$summarize_to_gene)
    
    summarize_to_gene <- input$summarize_to_gene
    print(summarize_to_gene)
    #print(head(sample_lookup_sub()))
    
    cpm_Samplesub <- cpm_counts()[,match(
      c("ID", "gene", as.character(sample_lookup_sub()$Sample)),
      colnames(cpm_counts())
    )]
    

    # Subset cpm_Samplesub() for the inputted genes
    if(summarize_to_gene == 'T')
      {
        cpm_temp <- cpm_Samplesub[cpm_Samplesub$gene %in% HM_genes_df()$Gene,]
        HM_genes_df_temp <- HM_genes_df()[HM_genes_df()$Gene %in% cpm_temp$gene,]
        
        # Want order maintained
        cpm_temp$gene <- factor(cpm_temp$gene,
                                levels = unique(HM_genes_df_temp$Gene))
        cpm_temp <- cpm_temp[order(cpm_temp$gene),] # This will be the same order as input gene-heatmap file
        
        
        # Only care about the gene column (of the first two)
        cpm_temp <- cpm_temp[,-1]
        cpm_temp <- ddply(cpm_temp, "gene", numcolwise(sum))
        rownames(cpm_temp) <- cpm_temp$gene
        cpm_temp <- cpm_temp[,-1]
        cpm_Samplesub_GeneSub <- log2(cpm_temp + 1)
        cpm_Samplesub_GeneSub <- na.omit(cpm_Samplesub_GeneSub) # Remove sgRNA or genes with no variance
    }
    
    else{
      cpm_temp <- cpm_Samplesub[cpm_Samplesub$gene %in% HM_genes_df()$Gene,]
      HM_genes_df_temp <- HM_genes_df()[HM_genes_df()$Gene %in% cpm_temp$gene,]
      # Want order maintained
      cpm_temp$gene <- factor(cpm_temp$gene,
                              levels = unique(HM_genes_df_temp$Gene))
      cpm_temp <- cpm_temp[order(cpm_temp$gene),] # This will be the same order as input gene-heatmap file
      rownames(cpm_temp) <- cpm_temp$ID
      cpm_temp <- cpm_temp[,-c(1,2)]
      cpm_Samplesub_GeneSub <- log2(cpm_temp + 1)
      cpm_Samplesub_GeneSub <- na.omit(cpm_Samplesub_GeneSub) # Remove sgRNA or genes with no variance
    }
    
    #print(head(cpm_Samplesub_GeneSub))
    return(cpm_Samplesub_GeneSub)
    
  })
  
  
  
  HM_genes_final <- reactive({
    req(CPM_table())
    req(HM_genes_df())
    
    if(summarize_to_gene == T){
      HM_genes_df_HM <- HM_genes_df()[HM_genes_df()$Gene %in% rownames(CPM_table()),]
    }else{
      HM_sgRNA_lookup <- cbind.data.frame(sgRNA = rownames(CPM_table()),
                                          Gene = sapply(strsplit(rownames(CPM_table()),
                                                                 split = "_"),
                                                        function(x){x[[1]]}))
      HM_genes_df_HM <- join(HM_sgRNA_lookup,
                             HM_genes_df())
    }
    
    
    return(HM_genes_df_HM)
  })
  
  
  
  ## Where z score transformation occurs
  HM_cpm <- reactive({
    req(CPM_table())
    req(HM_genes_df())
    
    cpm <- CPM_table()
    HM_genes_df_ <- HM_genes_df()
    
    cpm <- cpm[rownames(cpm) %in% HM_genes_df_$Gene,]
    
    # Match the order of the genes (file order applied to the normalized data)
    #genes_in_data <- toupper(HM_genes_df$Gene)[toupper(HM_genes_df$Gene) %in% rownames(cpm)]
    genes_in_data <- HM_genes_df_$Gene[HM_genes_df_$Gene %in% rownames(cpm)]
    cpm <- cpm[match(genes_in_data, rownames(cpm)),]
    
    # Optional Z-score
    if("Z" %in% input$HM_options){
      cpm <-data.frame(t(scale(t(cpm),scale = T,center=T)))
    }
    
    #print(head(cpm))
    return(cpm)
    
  })
  
  output$HM_col_ann_ui <- renderUI({
    req(selected_groups())
    
    df <- sample_lookup_sub()
    
    selectInput("HM_col_ann", "Optionally select sample annotation", choices = as.list(c("None", colnames(df)[2:ncol(df)])),
                selected = "None", multiple = F )
    # Automatically selects the second column
  })
  
  
  output$HM_row_ann_ui <- renderUI({
    df <- HM_genes_df()
    if (ncol(df) > 1){ # If there are additional annotations along with gene names:
      selectInput("HM_row_ann", "Optionally select gene annotation", choices = as.list(c("None", colnames(df)[2:ncol(df)])),
                  selected = "None", multiple = F )
    }else{
      selectInput("HM_row_ann", "Optionally select gene annotation", choices = as.list("None"),
                  selected = "None", multiple = F )
    }
  })
  
  
  ####    Several HM plotting parameters:
  HM_column_names <- reactive({
    if("colnames" %in% input$HM_options){
      T
    }else{
      F
    }
  })
  
  HM_row_names <- reactive({
    if("rownames" %in% input$HM_options){
      T
    }else{
      F
    }
  })
  
  HM_column_cluster <- reactive({
    if("col_cluster" %in% input$HM_options){
      T
    }else{
      F
    }
  })
  
  HM_row_cluster <- reactive({
    if("row_cluster" %in% input$HM_options){
      T
    }else{
      F
    }
  })
   
  # Make Heatmap column annotation file:
  Sample <- reactive({
    req(input$HM_col_ann)
    req(HM_genes_df())
    
    if(input$HM_col_ann == "None"){
      NULL
    }else{
      df <- sample_lookup()
      Samples <- as.vector(df[,input$HM_col_ann])
      samples <- data.frame(Samples) # Helps with asthetics if in two steps
      colnames(samples) <- "Sample"

      # Now, pick divergent colors:
      colorz <- setNames( distinctColorPalette(length(unique(Samples))), 
                          unique(Samples))
      col_HA <- HeatmapAnnotation(df = samples,
                                  col = list(Sample = colorz), # Name of column and the list of value = color
                                  annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
      col_HA
    }
  })
   
   
  Heatmap_object <- eventReactive(input$Heatmap_button,{
    
    cpm_table <- data.matrix(HM_cpm())

    HM <- Heatmap(cpm_table,
                  heatmap_legend_param = list(at = c(floor(min(cpm_table)), 0, ceiling(max(cpm_table)))),
                  cluster_rows = HM_row_cluster(), 
                  cluster_columns = HM_column_cluster(), 
                  #column_title = "METABRIC breast cancer samples",
                  column_title_side = "bottom", # name on bottom
                  name = "Z",
                  show_column_names = HM_column_names(), # INPUT PARAMETER
                  show_row_names = HM_row_names(), # INPUT PARAMETER
                  show_column_dend = HM_column_cluster(), # INPUT PARAMETER, Will match cluster_columns = ?
                  show_row_dend = HM_row_cluster(),
                  top_annotation = Sample()) 
    
    # if/else statement for row annotation here
    if(input$HM_row_ann == "None"){
      HM
    }else{
      df <- HM_genes_df()
      Genes <- as.vector(df[,input$HM_row_ann])
      genes <- data.frame(Genes) # Helps with aesthetics if in two steps
      colnames(genes) <- "Genes"
      # Now, pick divergent colors:
      colorz <- setNames( distinctColorPalette(length(unique(Genes))), 
                          unique(Genes))
      genes <- rowAnnotation(df = genes,
                             col = list(genes = colorz),
                             annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))  
      HM <- HM + genes
      HM
    }
  })
  
  output$Heatmap_view <- renderPlot({
      Heatmap_object()
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
  
  # observeEvent(input$timeOut, { 
  #   # Modified from: https://stackoverflow.com/questions/33839543/shiny-server-session-time-out-doesnt-work
  #   showModal(modalDialog(
  #     title = "Inactivity Timeout",
  #     paste("Session timeout due to", input$timeOut, "inactivity -", Sys.time()),
  #     footer = NULL
  #   ))
  #   stopApp() # I made this change so that the app closes instead of the window only
  # })
  # 

  
}

shinyApp(ui, server)