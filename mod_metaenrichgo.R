################################################################################
#                                                                              #
#                               Master’s Thesis                                #
#                                                                              #
#     Title: TFM: Development of the Meta-Analysis Tool MetaEnrichGO           #
#                                                                              #
#     Author: Maksym Kupchyk Tiurin                                            #
#     Supervisor: Francisco García García                 Year: 2025           #
#     Co-supervisor: Rubén Grillo Risco                                        #
#                                                                              #
#     Academic Institution: CIPF - Centro de Investigación Príncipe Felipe     #
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Loading of the required dependencies 
#-------------------------------------------------------------------------------

library(shiny)                  # enables building interactive web app
library(shinycssloaders)        # add a loading animation to outputs instead
library(shinyWidgets)           # custom input controls and user interface components
library(shinyBS)                # adds additional Twitter Bootstrap components 
library(shinyjs)                # perform common useful JavaScript operations
library(bslib)                  # toolkit for Shiny and R based on Bootstrap
library(dplyr)                  # tool for working with data.frame objects
library(DT)                     # interactive representation of data.frames
library(ggplot2)                # system for declaratively creating graphics
library(processx)               # run system processes in the background
library(plotly)                 # makes interactive, publication-quality graphs
library(zip)                    # allows dwnload more than one file at the same time

#------------------------------------------------------------------------------- 
# The user interface (ui) object controls the layout and appearance of your app
#-------------------------------------------------------------------------------

# R/modules/metaenrichgo/mod_metaEnrichGO.R

source("./modules/metaenrichgo/R/ORA.R")
source("./modules/metaenrichgo/R/meta-analysis.R")

mod_metaEnrichGO_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    tags$head(
      includeCSS("./modules/metaenrichgo/www/style.css")
    ),
    
    tags$head(
      tags$style(HTML("
    /* Wrap header text + icon */
    .dt-header-with-icon {
      display: inline-flex;
      align-items: center;
    }
    .dt-header-label {
      margin-right: 4px;
    }
    /* Little circle with white '?' */
    .dt-info-icon {
      display: inline-block;
      width: 16px;
      height: 16px;
      line-height: 16px;
      text-align: center;
      border-radius: 50%;
      background: #888;
      color: white;
      font-size: 12px;
      cursor: help;
    }
    /* Optional: give a slightly darker bg on hover */
    .dt-info-icon:hover {
      background: #555;
    }
  "))),
    
    page_sidebar(
      id    = ns("metaenrichgo_container"),
      title = "MetaEnrichGO",
      
      sidebar = sidebar(width = "380px",
                        
                        tabsetPanel(
                          id = ns("tabs"),
                          
                          tabPanel("Analysis", useShinyjs(),
                                   
                                   tags$h5("Input", style = "margin-top: 20px;"),
                                   
                                   radioButtons(inputId = ns("inputMethod"),
                                                label   = "Input method",
                                                choices = c("Upload Files" = "files",
                                                            "Paste Genes"  = "paste"),
                                                selected = "files"),
                                   
                                   conditionalPanel(
                                     condition = sprintf("input['%s'] == 'files'", ns("inputMethod")),
                                     fluidRow(
                                       column(10, fileInput(inputId = ns("files"),
                                                            label   = "Choose txt Files",
                                                            multiple= TRUE,
                                                            accept  = c(".txt", ".tsv"))),
                                       column(2, actionButton(inputId = ns("showExample"),
                                                              label   = "",
                                                              icon    = icon("info-circle", lib = "font-awesome"),
                                                              class   = "btn-info"))
                                     ),
                                     
                                     tags$p("Use example data?", style = "margin-bottom: -10px; margin-top: -15px;"),
                                     fluidRow(
                                       column(11,
                                              div(
                                                switchInput(inputId = ns("useExampleData"),
                                                            onLabel  = "Yes",
                                                            offLabel = "No",
                                                            value    = FALSE),
                                                title = HTML("These example human gene lists (4 studies related to lung cancer: GSE10072, GSE19188, GSE63459, GSE75037) contain one gene symbol per line, have no headers, and include duplicates and blank lines to simulate raw data."),
                                                style = "cursor: help;"
                                              )
                                       ),
                                       column(1,
                                              downloadLink(outputId = ns("downloadExample"),
                                                           label    = icon("download"),
                                                           style    = "margin-left:-28px; margin-top:10px;")
                                       )
                                     )
                                   ),
                                   
                                   conditionalPanel(
                                     condition = sprintf("input['%s'] == 'paste'", ns("inputMethod")),
                                     div(
                                       textAreaInput(inputId = ns("pastedGenes"),
                                                     label   = "Paste Gene Lists",
                                                     placeholder = "TP53\nMYC\nEGFR\nBRCA1\n...\n###\nBRCA1\nPTEN\nMYC\n...",
                                                     height  = "160px"),
                                       title = HTML("MetaEnrichGO Paste Format:\n- No header required; only gene symbols.\n- One gene per line, separated by '\\n'.\n- Separate multiple lists with '###' on its own line.\n- Do not include numeric values or additional columns.\n- Raw, unweighted gene rankings only."),
                                       style = "cursor: help;"
                                     )
                                   ),
                                   
                                   tags$hr(style = "margin-top: 20px; margin-bottom: 20px;"),
                                   
                                   tags$h5("Parameters"),
                                   
                                   selectInput(inputId = ns("database"),
                                               label    = "Database",
                                               choices  = c("Gene Ontology (GO)" = "GO",
                                                            "KEGG"              = "KEGG",
                                                            "Reactome"          = "Reactome"),
                                               selected = "GO"),
                                   
                                   conditionalPanel(
                                     condition = sprintf("input['%s'] == 'GO'", ns("database")),
                                     radioButtons(inputId = ns("ontology"),
                                                  label   = "Ontology",
                                                  choices = c("Biological Process (BP)"   = "BP",
                                                              "Molecular Function (MF)"  = "MF",
                                                              "Cellular Component (CC)"  = "CC"),
                                                  selected = "BP")
                                   ),
                                   
                                   selectInput(inputId = ns("organism"),
                                               label    = "Organism",
                                               choices  = c("Homo sapiens"         = "Hsa",
                                                            "Mus musculus"        = "Mmu",
                                                            "Rattus norvegicus"   = "Rno"),
                                               selected = "Hsa"),
                                   
                                   selectInput(inputId = ns("IDtype"),
                                               label    = "Gene ID",
                                               choices  = c("SYMBOL", "ENTREZID", "ENSEMBL"),
                                               selected = "SYMBOL"),
                                   
                                   div(
                                     selectInput(inputId = ns("metaAnalysisMethod"),
                                                 label    = "Meta-Analysis Method",
                                                 choices  = c("Fisher"    = "fisher",
                                                              "Stouffer"  = "stouffer",
                                                              "Tippett"   = "tippett",
                                                              "Wilkinson" = "wilkinson"),
                                                 selected = "fisher"),
                                     title = HTML("Combine p-values from multiple studies:\n- Fisher: Sum of log-transformed p-values, sensitive to small p-values.\n- Stouffer: Combines Z-scores, allows weighting.\n- Tippett: Uses the minimum p-value, emphasizes the most significant study.\n- Wilkinson: Generalized form of Tippett, based on order statistics."),
                                     style = "cursor: help;"
                                   ),
                                   
                                   div(
                                     sliderInput(inputId = ns("termThreshold"),
                                                 label    = "Minimum number of datasets",
                                                 value    = 1,
                                                 min      = 1,
                                                 max      = 4,
                                                 step     = 1),
                                     title = HTML("Gene Appearance Filter:\n - Require a gene to appear in at least this many lists.\n- E.g., if BRCA1 appears in 2/4 lists and you set this to 3, it will be excluded.\n- Setting to 2 will keep it."),
                                     style = "cursor: help;"
                                   )
                          ),
                          
                          tabPanel("Data visualization", br(),
                                   
                                   h5("Table columns"),
                                   
                                   checkboxGroupInput(
                                     inputId = ns("selectedColumns"),
                                     label   = NULL,
                                     choices   = c("ID", "Description", "pvalue", "p.adjust", 
                                                   "ListCount", "ListNames", "GeneCount", "GeneID"),
                                     selected  = c("ID", "Description", "pvalue", "p.adjust", 
                                                   "ListCount", "ListNames", "GeneCount", "GeneID"),
                                     inline    = FALSE),
                                   
                                   tags$hr(style = "margin-top: 20px; margin-bottom: 20px;"),
                                   
                                   h5("Plot"),
                                   
                                   sliderInput(inputId = ns("nPlotTerms"),
                                               label   = "Number of terms to show",
                                               min     = 1,
                                               max     = 20,
                                               value   = 10),
                                   
                                   selectInput(inputId = ns("plotLabel"),
                                               label   = "Y-axis label",
                                               choices = c("ID" = "ID",
                                                           "Description" = "Description",
                                                           "ID + Description" = "Term_Description"),
                                               selected = "Term_Description"),
                                   
                                   selectInput(inputId = ns("plotType"),
                                               label   = "Plot Type",
                                               choices = c("Dot Plot" = "dotplot",
                                                           "Bar Plot" = "barplot"),
                                               selected = "dotplot"),
                                   
                                   colorPickr(inputId = ns("colorLow"),
                                              label   = "Low color",
                                              selected = "#0838a0",
                                              opacity  = TRUE),
                                   
                                   colorPickr(inputId = ns("colorHigh"),
                                              label   = "High color",
                                              selected = "#2ba915",
                                              opacity  = TRUE),
                                   
                                   sliderInput(inputId = ns("textSize"),
                                               label   = "Text Size",
                                               min     = 5,
                                               max     = 20,
                                               value   = 12)
                          )
                        ),
                        
                        actionButton(inputId = ns("runButton"),
                                     label   = "Run Analysis",
                                     style   = "width: 100%; margin-top: 15px;")
      ),
      
      card(
        shinycssloaders::withSpinner(
          dataTableOutput(ns("resultsTable"))
        ), br(),
        
        uiOutput(ns("downloadButtonUI")), br(),
        
        shinyjs::hidden(
          div(id = ns("excludedSection"),
              h5("Excluded terms (below frequency threshold):"),
              dataTableOutput(ns("excludedTable"))
          )
        ), br(),
        
        shinycssloaders::withSpinner(
          plotlyOutput(ns("resultsDotplot"), height = "700px")
        ),
        
        uiOutput(ns("downloadPlotUI"))
      )
    )
  )
}



#-------------------------------------------------------------------------------

# R/modules/metaenrichgo/mod_metaEnrichGO.R

mod_metaEnrichGO_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    exampleFiles <- c("./modules/metaenrichgo/example_data/genes1.txt",
                      "./modules/metaenrichgo/example_data/genes2.txt",
                      "./modules/metaenrichgo/example_data/genes3.txt",
                      "./modules/metaenrichgo/example_data/genes4.txt")
    
    # --------------------------------------------
    # 1) Store uploaded or pasted input
    # --------------------------------------------
    storedInput <- reactiveVal(NULL)
    
    observe({
      req(input$inputMethod)
      
      if (input$inputMethod == "files") {
        if (isTRUE(input$useExampleData)) {
          storedInput(list(type = "files", paths = exampleFiles, names = basename(exampleFiles)))
        } else {
          req(input$files)
          storedInput(list(type = "files", paths = input$files$datapath, names = input$files$name))}
        
      } else if (input$inputMethod == "paste") {
        req(input$pastedGenes)
        storedInput(list(type = "paste", text = input$pastedGenes))}})
    
    
    # --------------------------------------------
    # 2) Read & clean raw gene lists with progress
    # --------------------------------------------
    rawGeneLists <- reactive({
      inputData <- storedInput()
      req(inputData)
      
      withProgress(message = "Loading gene lists...", value = 0, {
        
        incProgress(0.2, detail = "Reading input...")
        raw_lists <- list()
        
        if (inputData$type == "files") {
          for (i in seq_along(inputData$paths)) {
            fp <- inputData$paths[i]
            lines <- readLines(fp, warn = FALSE, encoding = "UTF-8")
            raw_lists[[ inputData$names[i] ]] <- lines}
        } else {
          blocks <- strsplit(inputData$text, "###", fixed = TRUE)[[1]]
          for (i in seq_along(blocks)) {
            raw_lists[[ paste0("List", i) ]] <- unlist(strsplit(blocks[[i]], "\\s+"))}}
        
        
        incProgress(0.7, detail = "Cleaning gene symbols...")
        clean_genes <- function(vec) {
          vec <- trimws(vec)
          vec <- gsub("///.*|\\\\.*", "", vec)
          vec[vec != ""]}
        
        mapping <- lapply(raw_lists, clean_genes)
        
        incProgress(1, detail = "Done.")
        mapping})})
    
    
    # --------------------------------------------
    # 3) Results reactiveVal & single Run handler
    # --------------------------------------------
    resultsVal <- reactiveVal(list(included = NULL, excluded = NULL))
    
    observeEvent(input$runButton, {
      inputData <- storedInput()
      req(inputData)
      
      # 1) File format
      if (inputData$type == "files" && !isTRUE(input$useExampleData)) {
        if (!validate_file_format_metaenrichgo(input, session)) {
          resultsVal(list(included = NULL, excluded = NULL))
          return()}}
      
      # 2) Paste format
      if (inputData$type == "paste") {
        if (!validate_paste_format(input, session)) {
          resultsVal(list(included = NULL, excluded = NULL))
          return()}}
      
      # 3) Clean lists & universe
      lists <- rawGeneLists()
      req(lists)
      genes <- unique(unlist(lists))
      
      # 4) Organism nomenclature
      if (!validate_organism(genes, input$organism, session)) {
        resultsVal(list(included = NULL, excluded = NULL))
        return()}
      
      # 5) Gene ID format
      if (!validate_gene_ids(genes, input$IDtype, session)) {
        resultsVal(list(included = NULL, excluded = NULL))
        return()}
      
      # 6) ORA step
      withProgress(message = "Performing ORA...", value = 0, {
        ORA_data <- switch(input$database,
                           GO       = lapply(lists, ORA_GO, organism = selectedOrganism(),
                                             ontology = input$ontology, pvalue = 1, ID_type = input$IDtype),
                           KEGG     = lapply(lists, ORA_KEGG, organism = selectedOrganism(),
                                             pvalue = 1, ID_type = input$IDtype),
                           Reactome = lapply(lists, ORA_REACTOME, organism = selectedOrganism(),
                                             pvalue = 1, ID_type = input$IDtype))
        incProgress(1, detail = "ORA Completed")})
      
      # 7) Meta‑analysis
      withProgress(message = "Performing Meta‑Analysis...", value = 0, {
        termLists     <- lapply(ORA_data, `[[`, "ID")
        termFrequency <- table(unlist(termLists))
        totalFiles    <- length(lists)
        minFreq       <- input$termThreshold
        commonTerms   <- names(termFrequency)[termFrequency >= minFreq]
        
        METAresults <- data.frame(
          ID          = character(), Description  = character(),
          pvalue      = numeric(),   p.adjust     = numeric(),
          ListCount   = character(),  ListNames   = character(),
          GeneCount   = integer(),   GeneID       = character(),
          stringsAsFactors = FALSE)
        
        totalTerms <- length(commonTerms)
        for (i in seq_along(commonTerms)) {
          term <- commonTerms[i]
          pvs  <- na.omit(sapply(ORA_data, function(res) {
            td <- res[res$ID == term, ]; if (nrow(td)) td$pvalue else NA}))
          
          combP <- switch(input$metaAnalysisMethod,
                          fisher    = fisher_combined_pvalue(pvs),
                          pearson   = pearson_combined_pvalue(pvs),
                          stouffer  = stouffer_combined_pvalue(pvs),
                          tippett   = tippett_combined_pvalue(pvs),
                          wilkinson = wilkinson_combined_pvalue(pvs))
          
          # Description
          desc <- NA
          for (res in ORA_data) {
            td <- res[res$ID == term, ]
            if (nrow(td)) { desc <- td$Description[1]; break }}
          
          # Genes in term
          genesIn <- unique(unlist(lapply(ORA_data, function(res) {
            td <- res[res$ID == term, ]; if (nrow(td)) strsplit(td$geneID, "/")[[1]]})))
          geneStr <- paste(genesIn, collapse = "/ ")
          filesWith <- names(lists)[sapply(lists, function(g) length(intersect(g, genesIn)) > 0)]
          
          METAresults <- rbind(METAresults, data.frame(
            ID          = term,
            Description = desc,
            pvalue      = combP,
            p.adjust    = NA,
            ListCount   = paste0(length(filesWith), "/", totalFiles),
            ListNames   = paste(filesWith, collapse = ", "),
            GeneCount   = length(genesIn),
            GeneID      = geneStr,
            stringsAsFactors = FALSE))
          
          incProgress(1 / totalTerms, detail = paste("Term", i, "of", totalTerms))}
        
        METAresults$p.adjust <- p.adjust(METAresults$pvalue, method = "BH")
        METAresults <- METAresults[order(METAresults$p.adjust), ]
        
        # Excluded
        allT <- names(termFrequency)
        lowT <- setdiff(allT, commonTerms)
        EXCres <- data.frame(
          ID          = character(), Description = character(),
          ListCount   = character(), ListNames   = character(),
          GeneCount   = integer(),   GeneID      = character(),
          stringsAsFactors = FALSE)
        
        for (term in lowT) {
          desc <- NA
          for (res in ORA_data) {
            td <- res[res$ID == term, ]
            if (nrow(td)) { desc <- td$Description[1]; break }}
          
          genesIn <- unique(unlist(lapply(ORA_data, function(res) {
            td <- res[res$ID == term, ]; if (nrow(td)) strsplit(td$geneID, "/")[[1]]})))
          
          filesWith <- names(lists)[sapply(lists, function(g) length(intersect(g, genesIn)) > 0)]
          EXCres <- rbind(EXCres, data.frame(
            ID          = term,
            Description = desc,
            ListCount   = paste0(length(filesWith), "/", totalFiles),
            ListNames   = paste(filesWith, collapse = ", "),
            GeneCount   = length(genesIn),
            GeneID      = paste(genesIn, collapse = "/ "),
            stringsAsFactors = FALSE))}})
      
      # 8) Empty results validation
      if (!check_empty_results(METAresults, session)) {
        resultsVal(list(included = NULL, excluded = NULL))
        return()}
      
      # 9) Store results
      resultsVal(list(included = METAresults, excluded = EXCres))})
    
    
    
    #-----------------------------------------------------------------------
    # Data representation
    #-----------------------------------------------------------------------
    
    output$resultsTable <- DT::renderDataTable({
      res <- resultsVal()$included
      req(res)
      req(input$selectedColumns)
      
      table_data <- res[, input$selectedColumns, drop = FALSE]
      
      if ("GeneID" %in% colnames(table_data)) {
        table_data$GeneID <- vapply(table_data$GeneID, function(txt) {
          esc <- htmltools::htmlEscape(txt)
          if (nchar(txt) > 60) {
            short <- substr(esc, 1, 60)
            sprintf('<span title="%s">%s...</span>', esc, short)
          } else {
            sprintf('<span title="%s">%s</span>', esc, esc)}
        }, character(1))}
      
      
      colTooltips <- list(
        ID = HTML(paste(
          "Stable identifier of the enriched term or pathway.",
          "Examples: GO accessions (e.g. GO:0008150), KEGG codes (e.g. 00010), Reactome IDs (e.g. R-HSA-199420).",
          sep = "\n")),
        
        Description = HTML(paste(
          "Descriptive name conveying the biological role of the term or pathway.",
          "E.g. “cellular response to chemical stimulus” (GO), “Glycolysis / Gluconeogenesis” (KEGG), etc.",
          sep = "\n")),
        
        pvalue = HTML(paste(
          "Combined p-value from enrichment results using Over-Representation Analysis (ORA).",
          "These values are aggregated across multiple input lists",
          "using methods like Fisher's, Stouffer's, Tippett's, or Wilkinson's test.",
          "Each method combines individual p-values from different terms or pathways",
          "to assess the overall statistical significance of a term.",
          sep = "\n")),
        
        p.adjust = HTML(paste(
          "Adjusted p-value controlling the False Discovery Rate (FDR) via Benjamini–Hochberg correction.",
          sep = "\n")),
        
        ListCount = HTML(paste(
          "Number of input lists that contributed genes to this term or pathway.",
          "E.g. a value of 3/4 means that 3 out of 4 input lists",
          "contained genes associated with the term or pathway.",
          sep = "\n")),
        
        ListNames = HTML(paste(
          "Names of the input lists or files that contributed genes to this term.",
          "Shows which sources provided relevant information for the enrichment.",
          sep = "\n")),
        
        GeneCount = HTML(paste(
          "Total number of background genes annotated to this term.",
          sep = "\n")),
        
        GeneID = HTML(paste(
          "List of geneid from the input that map to this term or pathway.",
          "Indicates which specific genes drive the enrichment.",
          sep = "\n")))
      
      
      num_rows <- input$resultsTable_length
      table_h  <- paste0(num_rows * 40, "px")
      
      DT::datatable(
        table_data,
        rownames  = FALSE,
        selection = "multiple",
        filter    = "top",
        escape    = FALSE,
        options   = list(
          pageLength    = 10,
          scrollX       = TRUE,
          scrollY       = table_h,
          autoWidth     = FALSE,
          columnDefs    = list(list(
            targets   = "_all",
            className = "dt-wrap",
            width     = paste0(round(100 / max(1, ncol(table_data))), "%"))),
          headerCallback = htmlwidgets::JS(
            "function(thead, data, start, end, display) {",
            "  var tooltips =", jsonlite::toJSON(colTooltips), ";",
            "  $('th', thead).each(function() {",
            "    var txt = $(this).text().trim();",
            "    if (tooltips[txt]) {",
            "      var label = $('<span>').addClass('dt-header-label').text(txt);",
            "      var icon  = $('<span>?</span>')",
            "                     .addClass('dt-info-icon')",
            "                     .attr('title', tooltips[txt]);",
            "      $(this).empty().append($('<div>').addClass('dt-header-with-icon')",
            "                               .append(label).append(icon));",
            "    }",
            "  });",
            "}"
          ))
      ) %>%
        formatSignif(
          columns = intersect(c("pvalue", "p.adjust"), colnames(table_data)),
          digits  = 4)})
    
    
    # Secondary table that shows the excluded terms based on their appereance  
    output$excludedTable <- DT::renderDataTable({
      res <- resultsVal()$included
      req(res)
      req(input$selectedColumns)
      exc <- resultsVal()$excluded
      
      cols_exc    <- c("ID", "Description", "ListCount", "ListNames", "GeneCount", "GeneID")
      table_data  <- exc[, cols_exc, drop = FALSE]
      num_rows    <- input$excludedTable_length
      table_h     <- paste0(num_rows * 40, "px")
      
      datatable(
        table_data,
        rownames      = FALSE,
        selection     = "multiple",
        filter        = "top",
        options       = list(
          pageLength  = 5,
          scrollX     = TRUE,
          scrollY     = table_h,
          autoWidth   = FALSE,
          columnDefs  = list(list(
            targets   = "_all",
            className = "dt-wrap",
            width     = paste0(round(100 / length(cols_exc)), "%"))))) %>%
        formatSignif(
          columns = intersect(c("pvalue", "p.adjust"), colnames(table_data)),
          digits  = 4)})
    
    
    # Button control
    output$downloadButtonUI <- renderUI({
      req(resultsVal()$included)
      fluidRow(
        column(5, downloadButton(ns("downloadTableCSV"),  "Download CSV", style = "width:100%")),
        column(5, downloadButton(ns("downloadTableTSV"),  "Download TSV", style = "width:100%")),
        column(1, actionButton(ns("toggleExcluded"),
                               label = NULL,
                               icon = icon("eye"),
                               style = "width:100%",
                               title = "Show/hide excluded terms")),
        column(1, downloadButton(ns("downloadExcludedTermsList"),
                                 label = NULL,
                                 icon = icon("download"),
                                 style = "width:100%",
                                 title = "Download excluded terms list (.txt)")))})
    
    
    # Show/hide with animation the secondary table
    observeEvent(input$toggleExcluded, {
      shinyjs::toggle(id = "excludedSection", anim = TRUE, animType = "slide", time = 0.5)})
    
    output$downloadTableCSV <- downloadHandler(
      filename = function() paste0("MetaEnrichGO_results_", Sys.Date(), ".csv"),
      content = function(file) {
        req(resultsVal(), input$selectedColumns, input$resultsTable_rows_all)
        data <- resultsVal()$included
        cols <- intersect(input$selectedColumns, colnames(data))
        filtered_data <- data[input$resultsTable_rows_all, cols, drop = FALSE]
        write.csv(filtered_data, file, row.names = FALSE)})
    
    output$downloadTableTSV <- downloadHandler(
      filename = function() paste0("MetaEnrichGO_results_", Sys.Date(), ".tsv"),
      content = function(file) {
        req(resultsVal(), input$selectedColumns, input$resultsTable_rows_all)
        data <- resultsVal()$included
        cols <- intersect(input$selectedColumns, colnames(data))
        filtered_data <- data[input$resultsTable_rows_all, cols, drop = FALSE]
        write.table(filtered_data, file, row.names = FALSE, sep = "\t")})
    
    
    # EXCLUDED, only gene list
    output$downloadExcludedTermsList <- downloadHandler(
      filename = function() {
        paste0("Excluded_MetaEnrichGO_", Sys.Date(), ".tsv")},
      content = function(file) {
        req(resultsVal())
        exc <- resultsVal()$excluded
        write.table(exc, file,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE)})
    
    
    # Control of showing and hidding some columns
    filtered_data <- reactive({
      req(resultsVal())
      if (!is.null(input$resultsTable_rows_all) &&
          length(input$resultsTable_rows_all) > 0) {
        resultsVal()$included[as.numeric(input$resultsTable_rows_all), ]
      } else {
        resultsVal()$included}})
    
    
    #------------------------
    # Reactive Plot creation
    #-----------------------
    
    plot_obj <- reactive({
      req(filtered_data())
      finalResultsFiltered <- filtered_data() %>%
        arrange(p.adjust) %>%
        head(input$nPlotTerms)
      
      text_size <- input$textSize
      finalResultsFiltered$label_y <- switch(input$plotLabel,
                                             "ID"               = finalResultsFiltered$ID,
                                             "Description"      = finalResultsFiltered$Description,
                                             "Term_Description" = paste(finalResultsFiltered$ID, finalResultsFiltered$Description, sep = " - "))
      
      if (input$plotType == "dotplot") {
        p <- ggplot(finalResultsFiltered, aes(x = GeneCount, y = reorder(label_y, -p.adjust))) +
          geom_point(aes(size = GeneCount, color = p.adjust)) +
          scale_color_gradient(low = input$colorLow, high = input$colorHigh) +
          labs(x = "GeneCount", y = NULL, color = "p.adjust", title = "Dotplot of the Enrichment Analysis") +
          theme_minimal() +
          theme(axis.text.y = element_text(size = text_size), 
                axis.text.x = element_text(size = text_size),
                axis.title.x = element_text(size = text_size + 1),
                axis.title.y = element_text(size = text_size + 1),
                legend.text = element_text(size = text_size - 3),
                legend.title = element_text(size = text_size, margin = margin(b = 20)),
                plot.title = element_text(size = text_size + 2, face = "bold", hjust = 0.5))
      } else if (input$plotType == "barplot") {
        p <- ggplot(finalResultsFiltered, aes(x = GeneCount, y = reorder(label_y, -p.adjust))) +
          geom_bar(stat = "identity", aes(fill = p.adjust)) +
          scale_fill_gradient(low = input$colorLow, high = input$colorHigh) +
          labs(x = "GeneCount", y = NULL, fill = "p.adjust", title = "Barplot of the Enrichment Analysis") +
          theme_minimal() +
          theme(axis.text.x = element_text(size = text_size), 
                axis.text.y = element_text(size = text_size),
                axis.title.x = element_text(size = text_size + 1),
                axis.title.y = element_text(size = text_size + 1),
                legend.text = element_text(size = text_size - 3),
                legend.title = element_text(size = text_size, margin = margin(b = 10)),
                plot.title = element_text(size = text_size + 2, face = "bold", hjust = 0.5))}
      p})
    
    
    # Show the interactive plot
    output$resultsDotplot <- renderPlotly({
      req(plot_obj())
      ggplotly(plot_obj())})
    
    # Multiple download buttons
    output$downloadPlotUI <- renderUI({
      req(resultsVal()$included)
      dropdownButton(
        label = "Download Plot",
        icon = icon("download"),
        circle = FALSE,
        status = "secondary",
        width = "300px",
        tags$h4("Download Format"),
        downloadButton(ns("downloadPlotHTML"), "HTML", class = "btn-block"),
        downloadButton(ns("downloadPlotPNG"), "PNG",   class = "btn-block"),
        downloadButton(ns("downloadPlotJPG"), "JPG",   class = "btn-block"))})
    
    # HTML
    output$downloadPlotHTML <- downloadHandler(
      filename = function() {
        paste("MetaEnrichGO_plot_", Sys.Date(), ".html", sep = "")},
      content = function(file) {
        saveWidget(ggplotly(plot_obj()), file, selfcontained = TRUE)})
    
    # PNG
    output$downloadPlotPNG <- downloadHandler(
      filename = function() {
        paste("MetaEnrichGO_plot_", Sys.Date(), ".png", sep = "")},
      content = function(file) {
        ggsave(file, plot = plot_obj(), device = "png",
               width = 12.5, height = 8.33, units = "in", dpi = 96)})
    
    # JPG
    output$downloadPlotJPG <- downloadHandler(
      filename = function() {
        paste("MetaEnrichGO_plot_", Sys.Date(), ".jpg", sep = "")},
      content = function(file) {
        ggsave(file, plot = plot_obj(), device = "jpeg",
               width = 12.5, height = 8.33, units = "in", dpi = 96)})
    
    
    
    ###############################################################################
    
    #----------------------- NAVIGATION FUNCTIONS --------------------------------#
    
    ###############################################################################
    
    
    # Show the example of file that the app accepts (blue "i" button)
    observeEvent(input$showExample, {
      showModal(
        modalDialog(
          title = tags$div("Example File Contents", style = "text-align: center; font-weight: bold; font-size: 20px;"),
          easy_close = TRUE,
          size = "l",
          fluidRow(
            column(12, 
                   tags$h5("Text Format (.txt)"),
                   HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>COL10A1\nCOL11A1\nGREM1\nMMP1\nMMP12\nSPINK1\nSPP1\nSULF1\nTOX3\nTOP2A</pre>"))),
          br(),
          tags$h4("📌 Important Notes", style = "font-weight: bold; font-size: 15px; color: #007bff;"),
          tags$div(style = "margin-top: 10px;",
                   tags$p(HTML("<b>• Accepted format:</b> Only <code>.txt</code> files are supported for this input type."),
                          style = "margin-bottom: 6px;"),
                   tags$p(HTML("<b>• One gene per line:</b> Each row should contain a single gene name. No additional values allowed."),
                          style = "margin-bottom: 6px;"),
                   tags$p(HTML("<b>• No headers:</b> Do <u>not</u> include any header row. Input should start directly with gene names."),
                          style = "margin-bottom: -6px;")),
          footer = modalButton("Close")))})
    
    
    # Slider information control that shows the number of lists and control que minimum appereance
    observe({
      req(input$inputMethod)
      total_lists <- NULL
      if (input$useExampleData) {
        total_lists <- 4
      } else if (input$inputMethod == "files" && !is.null(input$files)) {
        total_lists <- length(input$files$datapath)
      } else if (input$inputMethod == "paste" && !is.null(input$pastedGenes) && nzchar(input$pastedGenes)) {
        total_lists <- length(strsplit(input$pastedGenes, "###")[[1]])}
      
      if (!is.null(total_lists) && total_lists > 0) {
        isolate({
          updateSliderInput(session, "termThreshold",
                            min   = 1,
                            max   = total_lists,
                            step  = 1,
                            value = total_lists)})}})
    
    
    # Download button that allows to obtain the example data used in the app as ".zip"
    output$downloadExample <- downloadHandler(
      filename = function() {"MetaEnrich_example_data.zip"},
      content  = function(zipfile) {
        tmp <- tempdir()
        files_cp <- file.path(tmp, basename(exampleFiles))
        file.copy(exampleFiles, files_cp, overwrite = TRUE)
        
        owd <- setwd(tmp)
        on.exit(setwd(owd), add = TRUE)
        utils::zip(
          zipfile = zipfile,
          files   = basename(files_cp),
          extras  = "-j")},
      contentType = "application/zip")
    
    
    # Prepare some important data
    selectedOrganism <- reactive({input$organism})
    
    
    # INPUT CONTROL BY RESTARTING THE DATA 
    prev_input_method <- reactiveVal("files")
    observeEvent(input$inputMethod, {
      prev <- prev_input_method()
      
      if (input$inputMethod == "paste" && prev == "files") {
        # CLEAN THE FILES INFORMATION AFTER CHANCHING THE INPUT FORMAT
        shinyjs::reset("files")
        updateSwitchInput(session, "useExampleData", value = FALSE)}
      
      if (input$inputMethod == "files" && prev == "paste") {
        # CLEAN THE PASTE INFORMATION AFTER CHANCHING THE INPUT FORMAT
        updateTextInput(session, "pastedGenes", value = "")}
      
      prev_input_method(input$inputMethod)})
    
    
    #' Validate files for MetaEnrichGO
    #' @param input Shiny input object
    #' @param session Shiny session object
    #' @return TRUE if valid; FALSE if any file had wrong format
    
    validate_file_format_metaenrichgo <- function(input, session) {
      bad_files <- list()
      
      files <- if (input$inputMethod == "files") input$files$datapath else NULL
      file_names <- if (input$inputMethod == "files") input$files$name else NULL
      
      if (is.null(files)) return(TRUE)
      
      for (i in seq_along(files)) {
        lines <- readLines(files[i], warn = FALSE, encoding = "UTF-8")
        lines <- trimws(lines[nzchar(lines)]) 
        
        if (any(grepl("#", lines))) {
          bad_files[[file_names[i]]] <- "Contains invalid characters like '#'."
          next}
        
        if (any(grepl("[,\t]", lines))) {
          bad_files[[file_names[i]]] <- "Files must contain only one gene symbol per line (no commas or tabs)."
          next}}
      
      if (length(bad_files) > 0) {
        showModal(modalDialog(
          title = "⚠️ Format Errors Detected",
          size = "l",
          easyClose = TRUE,
          footer = modalButton("Close"),
          tags$ul(
            lapply(names(bad_files), function(f) {
              tags$li(HTML(paste("<b>", f, "</b>: ", bad_files[[f]])))})),
          tags$hr(),
          fluidRow(
            column(12, tags$h4("Expected Format"),
                   tags$p("Each uploaded file must contain a plain list of gene symbols, one per line, with no headers or separators."),
                   tags$h5("Valid Example"),
                   HTML("<pre>BRCA1\nTP53\nMYC\nEGFR</pre>"),
                   tags$p("Avoid using commas, tabs, or comment characters like '#'.")))))
        return(FALSE)}
      return(TRUE)}
    
    
    #' Validate pasted gene lists (one gene per line, lists separated by ###)
    #' @param input   Shiny input object
    #' @param session Shiny session object
    #' @return TRUE if valid; FALSE (and shows modal) if format is incorrect
    
    validate_paste_format <- function(input, session) {
      if (!(input$inputMethod == "paste" && nzchar(input$pastedGenes))) {
        return(TRUE)}
      
      txt       <- input$pastedGenes
      has_sep   <- grepl("###", txt, fixed = TRUE)
      bad_chars <- grepl("[,;|:\\t\\]", txt)
      
      if (!has_sep || bad_chars) {
        showModal(modalDialog(
          title     = tags$div("⚠️ Paste Format Error",
                               style = "text-align:center;font-weight:bold;font-size:20px;"),
          easyClose = TRUE, size = "l",
          
          tags$h5("Correct format"),
          HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>TP53\nMYC\nEGFR\nBRCA1\n###\nBRCA1\nPTEN\nMYC\n...</pre>"),
          br(),
          
          tags$h5("Possible Issues"),
          tags$div(style = "margin-top:10px;",
                   tags$p("• You must include '###' on its own line to separate multiple lists."),
                   tags$p("• Do NOT use commas (,), tabs, periods (.), hyphens (-), semicolons (;), or pipes (|) in gene lines."),
                   tags$p("• Paste only raw gene symbols, one per line.")),
          br(),
          tags$h5("Recommendations"),
          tags$div(style = "margin-top:10px;",
                   tags$p("• Use '###' to split each ranked list."),
                   tags$p("• Avoid CSV/TSV formats or pasted tables."),
                   tags$p("• Each line must contain exactly one gene symbol.")),
          
          footer = modalButton("Close")))
        
        shinyjs::disable("runButton")
        shinyjs::delay(100, shinyjs::enable("runButton"))
        return(FALSE)}
      return(TRUE)}
    
    
    #' Validate organism nomenclature
    #' @param geneLists Character vector of genes to validate
    #' @param organism  One of "Hsa", "Mmu", "Rno"
    #' @param session   Shiny session (para showModal)
    #' @return TRUE if valid; FALSE (and shows modal) if any invalid genes found
    
    validate_organism <- function(gene_list, organism, session) {
      if (length(gene_list) == 0) return(TRUE)
      
      pattern <- switch(organism,
                        Hsa = "^[A-Z][A-Z0-9._/\\-]*$",
                        Mmu = "^[A-Z][a-z0-9._/\\-]*$",
                        Rno = "^[A-Z][a-z0-9._/\\-]*$",
                        NULL)
      
      if (is.null(pattern)) return(TRUE)
      
      organism_names <- list(
        Hsa = "<i>Homo sapiens</i> (Human)",
        Mmu = "<i>Mus musculus</i> (Mouse)",
        Rno = "<i>Rattus norvegicus</i> (Rat)")
      
      invalid <- unique(gene_list[!grepl(pattern, gene_list)])
      if (length(invalid) > 0) {
        show_genes <- head(invalid, 10)
        extra_note <- if (length(invalid) > 10) paste0("… and ", length(invalid) - 10, " more") else ""
        
        showModal(modalDialog(
          title = tags$div("⚠️ Organism–Gene Nomenclature Issue",
                           style = "text-align:center; font-weight:bold; font-size:20px;"),
          easyClose = TRUE, size = "l",
          
          tags$h5(HTML(paste("Selected organism:", organism_names[[organism]])), style = "margin-bottom:10px;"),
          tags$p("The following gene names do not match the expected format for the selected organism:"),
          tags$pre(paste(c(show_genes, extra_note), collapse = "\n"),
                   style = "background-color:#f8f9fa; padding:10px; border-radius:5px;"),
          
          fluidRow(
            column(6,
                   tags$h5(HTML("Format for <i>Homo sapiens</i> (Hsa - Human)"), style = "margin-top:32px;"),
                   HTML("<pre style='background-color:#f8f9fa; padding:10px; margin-top:15px; border-radius:5px;'>BRCA1\nTPP1\nGREM1\nMMP12\nSPINK1</pre>"),
                   tags$p(HTML("<b>• Expected Format:</b> Fully uppercase (e.g., BRCA1, TPP1)."),
                          style = "margin-bottom: 6px;"),
                   tags$p(HTML("<b>• Official Nomenclature:</b> Gene names for humans are conventionally written in uppercase."),
                          style = "margin-bottom: 6px;")),
            column(6,
                   tags$h5(HTML("Format for <i>Mus musculus</i> (Mmu - Mouse) and <i>Rattus norvegicus</i> (Rno - Rat)"), style = "margin-top:15px;"),
                   HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>Brca1\nTpp1\nGrem1\nMmp12\nSpink1</pre>"),
                   tags$p(HTML("<b>• Expected Format:</b> First letter uppercase, followed by lowercase (e.g., Brca1, Tpp1)."),
                          style = "margin-bottom: 6px;"),
                   tags$p(HTML("<b>• Official Nomenclature:</b> Gene names for mice and rats follow this mixed-case format."),
                          style = "margin-bottom: 6px;"))),
          
          footer = modalButton("Close")))
        return(FALSE)}
      return(TRUE)}
    
    
    #' Validate geneid format
    #' @param gene_list Character vector of genes to validate
    #' @param ID_type  One of "SYMBOL", "ENSEMBL", "ENTREZID"
    #' @param session   Shiny session (para showModal)
    #' @return TRUE if valid; FALSE (and shows modal) if any invalid genes found
    
    validate_gene_ids <- function(gene_vector, ID_type, session) {
      if (length(gene_vector) == 0) return(TRUE)
      
      gene_vector <- trimws(gene_vector)
      gene_vector <- gene_vector[gene_vector != ""]
      
      pattern <- switch(ID_type,
                        SYMBOL   = "^[A-Za-z0-9._/\\-]+$",
                        ENTREZID = "^[0-9]+$",
                        ENSEMBL  = "^ENS[A-Z0-9]+$",
                        stop("Unsupported ID type"))
      
      invalid <- gene_vector[!grepl(pattern, gene_vector)]
      
      if (length(invalid) > 0) {
        show_genes <- head(invalid, 10)
        extra_note <- if (length(invalid) > 10) paste0("… and ", length(invalid) - 10, " more") else ""
        
        showModal(modalDialog(
          title = tags$div("⚠️ Gene ID Format Error",
                           style = "text-align: center; font-weight: bold; font-size: 20px;"),
          easyClose = TRUE,
          size = "l",
          
          tags$h4(paste("Selected ID type:", ID_type), style = "margin-bottom: 10px;"),
          tags$p("The following gene identifiers do not match the expected format:"),
          tags$pre(paste(c(show_genes, extra_note), collapse = "\n"),
                   style = "background-color:#f8f9fa; padding:10px; border-radius:5px;"),
          br(),
          tags$h4("Expected Formats", style = "font-weight:bold;"),
          fluidRow(
            column(4, tags$h5("SYMBOL"),
                   HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>BRCA1\t(or Brca1)\nTPP1\t(or Tpp1)\nGREM1\t(or Grem1)\n...</pre>")),
            column(4, tags$h5("ENTREZID"),
                   HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>12345\n45678\n89012\n...</pre>")),
            column(4, tags$h5("ENSEMBL"),
                   HTML("<pre style='background-color:#f8f9fa; padding:10px; border-radius:5px;'>ENSG00000139618\nENSMUSG00002125\nENSG00000189618\n...</pre>"))),
          br(),
          tags$h4("💡 Recommendation", style = "font-weight:bold;"),
          tags$p("Please convert your gene list to the selected ID type before running ORA."),
          footer = modalButton("Close")))
        return(FALSE) }
      return(TRUE)}
    
    
    #' Check if enrichment results are empty and show modal if so
    #' @param results  Data.frame or object with enrichment results
    #' @param session  Shiny session object
    #' @return TRUE if results are non-empty; FALSE (and shows modal) if empty
    
    check_empty_results <- function(results, session) {
      is_empty <- is.null(results) || (is.data.frame(results) && nrow(results) == 0)
      
      if (is_empty) {
        showModal(modalDialog(
          title = tags$div("⚠️ No Enrichment Results Found", 
                           style = "text-align: center; font-weight: bold; font-size: 20px;"),
          easyClose = TRUE,
          size = "l",
          tags$h4("No biological terms matched your input genes"),
          tags$p("The enrichment analysis did not find any relevant pathways, terms or functions associated with the provided gene list."),
          tags$hr(),
          tags$h4("Suggestions to improve your input"),
          tags$ul(
            tags$li("Check the organism selection (e.g., Human, Mouse or Rat)."),
            tags$li("Verify gene identifiers and formats (e.g., SYMBOL, ENTREZID or ENSEMBL)."),
            tags$li("Ensure your genes are annotated in the selected database."),
            tags$li("Use more genes if possible — small lists often yield no results.")),
          tags$p("If the problem persists, try example data to confirm the tool is working."),
          footer = modalButton("Close")))
        return(FALSE)}
      return(TRUE)}
    
    
  })}