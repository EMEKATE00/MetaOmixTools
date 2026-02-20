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

source("./modules/metaenrichgo/R/ORA.R")
source("./modules/metaenrichgo/R/meta-analysis.R")

mod_metaEnrichGO_ui <- function(id) {
  # Define el namespace para el módulo
  ns <- NS(id)
  
  tagList(
    tags$head(
      includeCSS("modules/metaenrichgo/www/style.css")),
    
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
    
    # Diseño principal de la página con una barra lateral
    page_sidebar(
      id = ns("metaenrichgo_sidebar"),  
      title = "MetaEnrichGO",
      
      sidebar = sidebar(
        width = "380px",
        
        # Paneles con pestañas para organizar los controles
        tabsetPanel(
          
          # Pestaña 1: Análisis
          tabPanel("Analysis", useShinyjs(),
                   
                   tags$h5("Input", style = "margin-top: 20px;"), 
                   
                   # Selección del método de entrada
                   radioButtons(inputId  = ns("inputMethod"),
                                label    = "Input method",
                                choices  = c("Upload Files" = "files",
                                             "Paste Genes" = "paste"),
                                selected = "files"),
                   
                   # Panel condicional para la carga de archivos
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'files'", ns("inputMethod")),
                     fluidRow(
                       column(10, fileInput(inputId   = ns("files"),
                                            label     = "Choose txt Files",
                                            multiple  = TRUE,
                                            accept    = c(".txt", ".tsv"))),
                       column(2, actionButton(inputId = ns("showExample"), 
                                              label   = "", 
                                              icon    = icon("info-circle", lib = "font-awesome"),
                                              class   = "btn-info"))),
                     
                     tags$p("Use example data?", style = "margin-top: -15px; margin-bottom: -10px;"),
                     fluidRow(
                       column(width = 11,
                              div(switchInput(inputId  = ns("useExampleData"),
                                              onLabel  = "Yes",
                                              offLabel = "No",
                                              value    = FALSE),
                                  title = HTML("These example human gene lists (4 studies related to lung cancer: GSE10072, GSE19188, GSE63459, GSE75037) contain one gene symbol per line, have no headers, and include duplicates and blank lines to simulate raw data."),
                                  style = "cursor: help;")),
                       
                       column(width = 1,
                              downloadLink(outputId = ns("downloadExample"),
                                           label    = icon("download"), 
                                           style    = "margin-left:-28px; margin-top:10px;")))),
                   
                   # Panel condicional para pegar genes
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'paste'", ns("inputMethod")),
                     div(textAreaInput(inputId     = ns("pastedGenes"),
                                       label       = "Paste Gene Lists",
                                       placeholder = "TP53\nMYC\nEGFR\nBRCA1\n...\n###\nBRCA1\nPTEN\nMYC\n...",
                                       height      = "160px"),
                         title = HTML("MetaEnrichGO Paste Format:\n- No header required; only gene symbols.\n- One gene per line, separated by '\\n'.\n- Separate multiple lists with '###' on its own line.\n- Do not include numeric values or additional columns.\n- Raw, unweighted gene rankings only."),
                         style = "cursor: help;")),
                   
                   tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
                   
                   
                   tags$h5("Parameters"), 
                   
                   # Selección de la base de datos
                   selectInput(inputId  = ns("database"),
                               label    = "Database",
                               choices  = c("Gene Ontology (GO)" = "GO",
                                            "KEGG"               = "KEGG",
                                            "Reactome"           = "Reactome",
                                            "Personalized"       = "Personalized"),
                               selected = "GO"),
                   
                   # Panel condicional para opciones de GO
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'GO'", ns("database")),
                     radioButtons(inputId    = ns("ontology"),
                                  label      = "Ontology",
                                  choices    = c("Biological Process (BP)" = "BP",
                                                 "Molecular Function (MF)" = "MF",
                                                 "Cellular Component (CC)" = "CC"),
                                  selected   = "BP")),
                   
                   # Panel condicional para anotaciones personalizadas
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'Personalized'", ns("database")),
                     fluidRow(
                       column(width = 10,
                              fileInput(inputId = ns("personalizedAnnotations"), 
                                        label = "Choose annotation file",
                                        multiple = FALSE, 
                                        accept = c(".txt", ".tsv", ".csv"))),
                       column(width = 2,
                              actionButton(inputId = ns("showPersonalizedInfo"), 
                                           label = "",
                                           icon = icon("info-circle", lib = "font-awesome"),
                                           class = "btn-info"),
                              
                              tags$div(
                                downloadLink(outputId = ns("downloadPersonalizedExample"),
                                             label = icon("download"),
                                             style = "margin-left: 4px; font-size: 16px;"),
                                style = "margin-top: 10px;")))),
                   
                   # Panel condicional para bases de datos estándar
                   conditionalPanel(
                     condition = sprintf("input['%s'] != 'Personalized'", ns("database")),
                     selectInput(inputId  = ns("organism"),
                                 label    = "Organism",
                                 choices  = c("Homo sapiens" = "Hsa",
                                              "Mus musculus" = "Mmu",
                                              "Rattus norvegicus" = "Rno"),
                                 selected = "Hsa"),
                     
                     selectInput(inputId  = ns("IDtype"),
                                 label    = "Gene ID",
                                 choices  = c("SYMBOL", "ENTREZID", "ENSEMBL"),
                                 selected = "SYMBOL")),
                   
                   # Selección del método de metaanálisis
                   div(selectInput(inputId  = ns("metaAnalysisMethod"),
                                   label    = "Meta-Analysis Method",
                                   choices  = c("Fisher"    = "fisher",
                                                "Stouffer"  = "stouffer",
                                                "Tippett"   = "tippett",
                                                "Wilkinson" = "wilkinson"),
                                   selected = "fisher"),
                       title = HTML("Combine p-values from multiple studies:\n- Fisher: Sum of log-transformed p-values, sensitive to small p-values.\n- Stouffer: Combines Z-scores, allows weighting.\n- Tippett: Uses the minimum p-value, emphasizes the most significant study.\n- Wilkinson: Generalized form of Tippett, based on order statistics."), 
                       style = "cursor: help;"),
                   
                   # Slider para el umbral de aparición de genes
                   div(sliderInput(inputId = ns("termThreshold"),
                                   label   = "Minimum number of datasets",
                                   value   = 1,
                                   min     = 1,
                                   max     = 4,
                                   step    = 1),
                       title = HTML(
                         "Term Appearance Filter:\n - Require a term to appear in at least this many lists.\n- E.g., if 'nucleus' appears in 2/4 lists and you set this to 3, it will be excluded.\n- Setting to 2 will keep it."),
                       style = "cursor: help;")),
          
          # Pestaña 2: Visualización de datos
          tabPanel("Data visualization", br(),
                   
                   h5("Table columns"),
                   
                   # Checkbox para seleccionar columnas de la tabla
                   checkboxGroupInput(
                     inputId  = ns("selectedColumns"),
                     label    = NULL,
                     choices  = c("ID", "Description", "pvalue", "p.adjust", 
                                  "ListCount", "ListNames", "GeneCount", "GeneID"),
                     selected = c("ID", "Description", "pvalue", "p.adjust", 
                                  "ListCount", "ListNames", "GeneCount", "GeneID"),
                     inline   = FALSE),
                   
                   tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
                   
                   h5("Plot"),
                   
                   # Controles para la personalización del gráfico
                   sliderInput(inputId  = ns("nPlotTerms"),
                               label    = "Number of terms to show",
                               min      = 1,
                               max      = 100,
                               value    = 10),
                   
                   selectInput(inputId  = ns("plotLabel"),
                               label    = "Y-axis label",
                               choices  = c("ID" = "ID",
                                            "Description" = "Description",
                                            "ID + Description" = "Term_Description"),
                               selected = "Term_Description"),
                   
                   selectInput(inputId  = ns("plotType"),
                               label    = "Plot Type",
                               choices  = c("Dot Plot" = "dotplot",
                                            "Bar Plot"  = "barplot"),
                               selected = "dotplot"),
                   
                   colorPickr(inputId   = ns("colorLow"),
                              label     = "Low color",
                              selected  = "#0838a0",
                              opacity   = TRUE),
                   
                   colorPickr(inputId   = ns("colorHigh"),
                              label     = "High color",
                              selected  = "#2ba915",
                              opacity   = TRUE),
                   
                   sliderInput(inputId  = ns("textSize"),
                               label    = "Text Size",
                               min      = 5,
                               max      = 20,
                               value    = 12))),
        
        # Botón para ejecutar el análisis
        actionButton(inputId = ns("runButton"),
                     label   = "Run Analysis",
                     style   = "width: 100%; margin-top: 10px;")),
      
      # Contenedor principal para mostrar los resultados
      card(uiOutput(ns("resultsTabs")))))}



#-------------------------------------------------------------------------------

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
    
    selectedDatabase <- eventReactive(input$runButton, {
      input$database})
    
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
    
    # --------------------------------------------------------------------------
    # 2) Read & clean raw gene lists
    # --------------------------------------------------------------------------
    rawGeneLists <- reactive({
      inputData <- storedInput()
      req(inputData)
      
      withProgress(message = "Processing input data...", value = 0, {
        lists <- list()
        
        # --- CASO 1: ARCHIVOS ---
        if (inputData$type == "files") {
          n_files <- length(inputData$paths)
          
          for (i in seq_along(inputData$paths)) {
            incProgress(1/n_files, detail = paste("Reading file:", inputData$names[i]))
            
            lines <- readLines(inputData$paths[i], warn = FALSE, encoding = "UTF-8")
            genes <- trimws(lines)
            genes <- genes[nzchar(genes)] 
            genes <- unique(genes)        
            
            if (length(genes) > 0) {
              clean_name <- tools::file_path_sans_ext(inputData$names[i])
              lists[[clean_name]] <- genes}}
          
          # --- CASO 2: TEXTO PEGADO ---
        } else {
          incProgress(0.3, detail = "Parsing pasted text...")
          
          raw_text <- inputData$text
          # 1. Dividir por el separador de listas ###
          parts <- strsplit(raw_text, "###", fixed = TRUE)[[1]]
          parts <- trimws(parts)
          parts <- parts[nzchar(parts)] 
          
          n_parts <- length(parts)
          
          if (n_parts > 0) {
            for (i in seq_along(parts)) {
              incProgress(0.7/n_parts, detail = paste("Processing list", i))
              
              # 2. Dividir por saltos de línea
              content <- parts[i]
              lines <- unlist(strsplit(content, "[\r\n]+"))
              
              genes <- trimws(lines)
              genes <- genes[nzchar(genes)]
              genes <- unique(genes)
              
              # Eliminar posibles encabezados tipo "Gene", "ID" si se colaron
              genes <- genes[!tolower(genes) %in% c("gene", "geneid", "symbol", "id")]
              
              if (length(genes) > 0) {
                lists[[paste0("List", i)]] <- genes}}}}
        
        incProgress(1, detail = "Finalizing...")
        
        if (length(lists) == 0) {
          return(NULL) }
        
        return(lists)})})
    
    
    # --------------------------------------------------------------------------
    # Helper: compute GeneCount column
    # --------------------------------------------------------------------------
    compute_gene_count <- function(geneid_vec) {
      vapply(geneid_vec, function(txt) {
        if (is.na(txt) || txt == "" ) return(0L)
        genes <- unlist(strsplit(as.character(txt), "[/\\s]+"))
        genes <- genes[genes != ""]
        length(unique(genes))
      }, integer(1))}
    
    
    # --------------------------------------------------------------------------
    # totalLists reactive (CRITICAL FIX)
    # --------------------------------------------------------------------------
    totalLists <- reactive({
      req(rawGeneLists())
      length(rawGeneLists())})
    
    
    # --------------------------------------------
    # 3) Results reactiveVal & runContext
    # --------------------------------------------
    resultsVal <- reactiveVal(list(included = NULL, excluded = NULL))
    runContext  <- reactiveVal(NULL)   
    
    observeEvent(input$runButton, {
      inputData <- storedInput()
      
      ###' @INPUT_NO_DATA_Validation
      if (!validate_input_exists(inputData)) {
        resultsVal(list(included = NULL, excluded = NULL))
        return()}
      
      req(inputData)
      n_lists <- totalLists()
      lists   <- rawGeneLists()
      genes   <- unique(unlist(lists))
      
      withProgress(message = "Initializing analysis...", value = 0, {
        
        # ----------------------------------------------------------
        # 1) Step: Validations
        # ----------------------------------------------------------
        incProgress(0.10, detail = "Validating input structure and formatting...")
        
        is_personalized <- (input$database == "Personalized")
        
        ###' @DATABASE_FILE_CHECK
        if (!check_database_exists(session, input$database, input$personalizedAnnotations)) {
          resultsVal(list(included = NULL, excluded = NULL))
          return()}
        
        ###' @ANNOTATION_FILE_FORMAT_CHECK
        if (is_personalized) {
          incProgress(0.05, detail = "Checking personalized annotation file format...")
          if (!check_annotation_file_format(input$personalizedAnnotations$datapath, session)) {
            resultsVal(list(included = NULL, excluded = NULL))
            return()}}
        
        # Column selections for future display
        incProgress(0.05, detail = "Preparing result table structure...")
        if (is_personalized) {
          all_cols <- c("Description", "pvalue", "p.adjust", "ListCount", "ListNames",
                        "FeatureCount", "FeatureID")
          default_selected <- if (n_lists > 1) all_cols else setdiff(all_cols, c("ListCount", "ListNames"))
          updateSelectInput(session, "plotLabel",
                            choices = c("Description" = "Description"),
                            selected = "Description")
        } else {
          all_cols <- c("ID", "Description", "pvalue", "p.adjust", "ListCount", "ListNames",
                        "GeneCount", "GeneID")
          default_selected <- if (n_lists > 1) all_cols else setdiff(all_cols, c("ListCount", "ListNames"))
          updateSelectInput(session, "plotLabel",
                            choices = c("ID" = "ID",
                                        "Description" = "Description",
                                        "ID + Description" = "Term_Description"),
                            selected = "Term_Description")
        }
        
        updateCheckboxGroupInput(session, "selectedColumns",
                                 choices  = all_cols,
                                 selected = default_selected)
        
        ###' @FILE_FORMAT_ERROR 
        if (inputData$type == "files" && !isTRUE(input$useExampleData)) {
          incProgress(0.05, detail = "Validating file format...")
          if (!validate_file_format_metaenrichgo(input, session)) {
            resultsVal(list(included = NULL, excluded = NULL))
            return()}}
        
        ###' @PASTE_FORMAT_ERROR
        if (inputData$type == "paste") {
          incProgress(0.05, detail = "Checking pasted gene list format...")
          if (!validate_paste_format(input, session)) {
            resultsVal(list(included = NULL, excluded = NULL))
            return()}}
        
        ### OTHER VALIDATION (IDs & Organism)
        if (!is_personalized) {
          incProgress(0.10, detail = "Validating gene identifiers...")
          
          ###' @GENE_ID_FORMAT_ERROR 
          checked_genes_ids <- validate_gene_ids(genes, input$IDtype, session)
          if (is.null(checked_genes_ids)) {
            resultsVal(list(included = NULL, excluded = NULL))
            return()}
          genes <- checked_genes_ids 
          
          ###' @ORGANISM_FORMAT_ERROR 
          incProgress(0.05, detail = "Checking organism compatibility...")
          checked_genes_org <- validate_organism(genes, input$organism, input$IDtype, session)
          if (is.null(checked_genes_org)) {
            resultsVal(list(included = NULL, excluded = NULL))
            return()}
          genes <- checked_genes_org}
        
        
        # ----------------------------------------------------------
        # 2) ORA computation
        # ----------------------------------------------------------
        incProgress(0.05, detail = "Starting ORA computation for all lists...")
        
        ORA_raw <- switch(input$database,
                          GO       = lapply(lists, ORA_GO, 
                                            organism = input$organism, 
                                            ontology = input$ontology, 
                                            pvalue = 1, 
                                            ID_type = input$IDtype),
                          KEGG     = lapply(lists, ORA_KEGG, 
                                            organism = input$organism, 
                                            pvalue = 1, 
                                            ID_type = input$IDtype),
                          Reactome = lapply(lists, ORA_REACTOME, 
                                            organism = input$organism, 
                                            pvalue = 1, 
                                            ID_type = input$IDtype),
                          Personalized = lapply(lists, ORA_PERSONALIZED, 
                                                annotation_file = input$personalizedAnnotations$datapath))
        
        ###' @MAPPING_FAILURE_CHECK 
        if (!check_mapping_failure(ORA_raw, session)) {
          resultsVal(list(included = NULL, excluded = NULL))
          return()}
        
        ORA_data <- withProgress(message = "Processing results...", value = 0, {
          total_l <- length(lists)
          processed <- list()
          
          for (i in seq_along(ORA_raw)) {
            incProgress(1/total_l, detail = paste("Formatting results", i, "of", total_l, "..."))
            
            df <- ORA_raw[[i]]
            
            if (is.null(df) || nrow(df) == 0) {
              df <- data.frame(ID=character(), Description=character(), pvalue=numeric(), p.adjust=numeric(), 
                               GeneID=character(), GeneCount=character(), stringsAsFactors=FALSE)
            } else {
              if (is_personalized) {
                df$ID <- df$Description
                df$GeneID <- df$FeatureID
                df$GeneCount <- df$FeatureCount
                df$FeatureID <- NULL; df$FeatureCount <- NULL
              } else {
                if ("geneID" %in% colnames(df) && !("GeneID" %in% colnames(df))) colnames(df)[colnames(df) == "geneID"] <- "GeneID"
                if (!("GeneID" %in% colnames(df))) df$GeneID <- ""
                if (!("GeneCount" %in% colnames(df))) df$GeneCount <- compute_gene_count(df$GeneID)}
              
              if (!("pvalue" %in% colnames(df))) df$pvalue <- NA_real_
              if (!("Description" %in% colnames(df))) df$Description <- NA_character_
              if (!("ID" %in% colnames(df))) df$ID <- NA_character_
              if (!("p.adjust" %in% colnames(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")}
            processed[[i]] <- df}
          processed})
        
        if (is.null(names(lists)) || any(names(lists) == ""))
          names(ORA_data) <- paste0("List", seq_along(ORA_data))
        else
          names(ORA_data) <- names(lists)
        
        runContext(list(ORA_data = ORA_data, lists = lists))
        
        # ----------------------------------------------------------
        # 3) Meta-Analysis (if >1 list)
        # ----------------------------------------------------------
        if (n_lists > 1) {
          
          META <- withProgress(message = "Performing Meta-Analysis...", value = 0, {
            
            termLists <- lapply(ORA_data, `[[`, "ID")
            termFrequency <- table(unlist(termLists))
            totalFiles    <- n_lists
            minFreq       <- as.numeric(input$termThreshold)
            commonTerms   <- names(termFrequency)[termFrequency >= minFreq]
            totalTerms    <- length(commonTerms)
            
            incProgress(0.05, detail = paste("Found", totalTerms, "shared terms..."))
            
            METAresults <- data.frame(
              ID = character(0), Description = character(0),
              pvalue = numeric(0), p.adjust = numeric(0),
              ListCount = character(0), ListNames = character(0),
              GeneCount = integer(0), GeneID = character(0),
              stringsAsFactors = FALSE)
            
            for (i in seq_along(commonTerms)) {
              incProgress(0.90 / max(totalTerms, 1),
                          detail = paste("Combining p-values for term", i, "of", totalTerms, "..."))
              
              term <- commonTerms[i]
              pvs <- na.omit(sapply(ORA_data, function(res) {
                td <- res[res$ID == term, ]
                if (nrow(td)) td$pvalue else NA}))
              
              combP <- switch(input$metaAnalysisMethod,
                              fisher    = fisher_combined_pvalue(pvs),
                              stouffer  = stouffer_combined_pvalue(pvs),
                              tippett   = tippett_combined_pvalue(pvs),
                              wilkinson = wilkinson_combined_pvalue(pvs))
              
              desc <- NA_character_
              for (res in ORA_data) {
                td <- res[res$ID == term, ]
                if (nrow(td)) { desc <- td$Description[1]; break }}
              
              genesIn <- unique(unlist(lapply(ORA_data, function(res) {
                td <- res[res$ID == term, ]
                if (nrow(td)) strsplit(as.character(td$GeneID), "[/\\s]+")[[1]] else character(0)})))
              genesIn <- genesIn[genesIn != ""]
              
              filesWith <- names(ORA_data)[sapply(ORA_data, function(res) term %in% res$ID)]
              
              METAresults <- rbind(METAresults, data.frame(
                ID = term,
                Description = desc,
                pvalue = combP,
                p.adjust = NA_real_,
                ListCount = paste0(length(filesWith), "/", totalFiles),
                ListNames = paste(filesWith, collapse = ", "),
                GeneCount = length(genesIn),
                GeneID = paste(genesIn, collapse = "/ "),
                stringsAsFactors = FALSE))}
            
            METAresults$p.adjust <- p.adjust(METAresults$pvalue, method = "BH")
            METAresults <- METAresults[order(METAresults$p.adjust), ]
            
            ###' @STRICT_THRESHOLD_CHECK
            if (!check_strict_threshold(METAresults, termFrequency, minFreq, session)) {
              return(NULL) }
            
            allT <- names(termFrequency)
            lowT <- setdiff(allT, commonTerms)
            
            EXCres <- data.frame(
              ID = character(0), Description = character(0),
              ListCount = character(0), ListNames = character(0),
              GeneCount = integer(0), GeneID = character(0),
              stringsAsFactors = FALSE)
            
            for (term in lowT) {
              tdDesc <- NA
              for (res in ORA_data) {
                td <- res[res$ID == term, ]
                if (nrow(td)) { tdDesc <- td$Description[1]; break }}
              
              genesIn <- unique(unlist(lapply(ORA_data, function(res) {
                td <- res[res$ID == term, ]
                if (nrow(td)) strsplit(as.character(td$GeneID), "[/\\s]+")[[1]] else character(0)})))
              genesIn <- genesIn[genesIn != ""]
              
              filesWith <- names(ORA_data)[sapply(ORA_data, function(res) term %in% res$ID)]
              
              EXCres <- rbind(EXCres, data.frame(
                ID = term,
                Description = tdDesc,
                ListCount = paste0(length(filesWith), "/", totalFiles),
                ListNames = paste(filesWith, collapse = ", "),
                GeneCount = length(genesIn),
                GeneID = paste(genesIn, collapse = "/ "),
                stringsAsFactors = FALSE))}
            
            list(included = METAresults, excluded = EXCres)})
          
          if (is.null(META)) {
            resultsVal(list(included = NULL, excluded = NULL)) 
            return() }
          
          METAresults <- META$included
          EXCres      <- META$excluded
          
          if (is_personalized) {
            colnames(METAresults)[colnames(METAresults) == "GeneCount"] <- "FeatureCount"
            colnames(METAresults)[colnames(METAresults) == "GeneID"] <- "FeatureID"
            METAresults$ID <- NULL
            
            colnames(EXCres)[colnames(EXCres) == "GeneCount"] <- "FeatureCount"
            colnames(EXCres)[colnames(EXCres) == "GeneID"] <- "FeatureID"
            EXCres$ID <- NULL}
          

          resultsVal(list(included = METAresults, excluded = EXCres))
          return()}
        
        # ----------------------------------------------------------
        # 4) Single-list case
        # ----------------------------------------------------------
        incProgress(0.20, detail = "Processing final results for single list...")
        df <- ORA_data[[1]]
        
        if ("geneID" %in% colnames(df) && !("GeneID" %in% colnames(df)))
          colnames(df)[colnames(df) == "geneID"] <- "GeneID"
        if (!("p.adjust" %in% colnames(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
        if (!("GeneCount" %in% colnames(df))) df$GeneCount <- compute_gene_count(df$GeneID)
        
        list_name <- if (!is.null(names(lists)) && length(names(lists)) >= 1) names(lists)[1] else "List1"
        df$ListCount <- "1/1"
        df$ListNames <- list_name
        
        if (is_personalized) {
          colnames(df)[colnames(df) == "GeneCount"] <- "FeatureCount"
          colnames(df)[colnames(df) == "GeneID"]    <- "FeatureID"
          df$ID <- NULL}
        
        expected <- all_cols
        missing  <- setdiff(expected, colnames(df))
        if (length(missing) > 0)
          for (cname in missing) df[[cname]] <- NA
        
        df <- df[, intersect(expected, colnames(df)), drop = FALSE]
        
        exc_cols <- if (is_personalized) {
          c("Description", "ListCount", "ListNames", "FeatureCount", "FeatureID")
        } else {
          c("ID", "Description", "ListCount", "ListNames", "GeneCount", "GeneID")}
        
        EXC_empty <- df[0, exc_cols, drop = FALSE]
        
        resultsVal(list(included = df, excluded = EXC_empty))
      }) # end global withProgress
    })
    
    
    
    #-----------------------------------------------------------------------
    # Generación dinámica de pestañas (Combined + 1 pestaña por archivo)
    #-----------------------------------------------------------------------
    output$resultsTabs <- renderUI({
      rc <- runContext()
      req(rc)
      ORA_list <- rc$ORA_data
      list_names <- names(ORA_list)
      if (is.null(list_names) || length(list_names) == 0) {
        list_names <- paste0("List", seq_along(ORA_list))}
      
      tabs <- list()
      
      # --- PESTAÑA COMBINED ---
      if (length(list_names) > 1) {
        tabs[[length(tabs) + 1]] <- tabPanel(
          "Combined", br(),
          shinycssloaders::withSpinner(DTOutput(ns("resultsTable"))), br(),
          uiOutput(ns("downloadButtonUI")), br(), 
          shinyjs::hidden(
            div(id = ns("excludedSection"),
                h5("Excluded terms (below frequency threshold):"),
                DTOutput(ns("excludedTable")))
          ), br(),
          uiOutput(ns("resultsDotplotContainer")), 
          br(),
          uiOutput(ns("downloadPlotUI")))}
      
      # --- PESTAÑAS INDIVIDUALES ---
      for (nm in list_names) {
        safe <- gsub("[^[:alnum:]_]", "_", nm)
        
        tabs[[length(tabs) + 1]] <- tabPanel(
          title = nm, br(),
          shinycssloaders::withSpinner(DTOutput(ns(paste0("resultsTable_ind_", safe)))), br(),
          fluidRow(
            column(6, downloadButton(ns(paste0("downloadTableCSV_ind_", safe)), "Download CSV", style = "width:100%")),
            column(6, downloadButton(ns(paste0("downloadTableTSV_ind_", safe)), "Download TSV", style = "width:100%"))
          ), br(),
          uiOutput(ns(paste0("resultsDotplotContainer_ind_", safe))),
          br(),
          uiOutput(ns(paste0("downloadPlotUI_ind_", safe))))}
      
      do.call(tabsetPanel, c(id = ns("dynamicResults"), tabs))})
    
    
    output$resultsDotplotContainer <- renderUI({
      req(plot_obj())
      plot_height <- max(600, input$nPlotTerms * 25)
      shinycssloaders::withSpinner(
        plotlyOutput(ns("resultsDotplot"), height = paste0(plot_height, "px")))})
    
    #-----------------------------------------------------------------------
    # Crear outputs dinámicos por cada archivo usando runContext()
    #-----------------------------------------------------------------------
    observeEvent(runContext(), {
      rc <- runContext()
      req(rc)
      ORA_list <- rc$ORA_data
      list_names <- names(ORA_list)
      if (is.null(list_names) || length(list_names) == 0) {
        list_names <- paste0("List", seq_along(ORA_list))
        names(ORA_list) <- list_names}
      
      db_type <- selectedDatabase()
      is_personalized <- !is.null(db_type) && db_type == "Personalized"
      
      if (is_personalized) {
        colTooltips <- list(
          Description = HTML("Descriptive name for the user-defined category."),
          pvalue = HTML("p-value from the Over-Representation Analysis."),
          p.adjust = HTML("Adjusted p-value using Benjamini–Hochberg."),
          FeatureCount = HTML("Number of features from your list mapped to this term."),
          FeatureID = HTML("List of feature IDs from your list mapped to this term."))
      } else {
        colTooltips <- list(
          ID = "Stable identifier of the enriched term or pathway.",
          Description = "Descriptive name conveying the biological role of the term.",
          pvalue = "p-value from the Over-Representation Analysis.",
          p.adjust = "Adjusted p-value using Benjamini–Hochberg.",
          GeneCount = "Number of genes from your list mapped to this term.",
          GeneID = "List of gene IDs from your list that map to this term.")}
      
      compute_gene_count <- function(vec) {
        vapply(vec, function(txt) {
          if (is.na(txt) || txt == "") return(0L)
          genes <- unlist(strsplit(as.character(txt), "[/\\s]+"))
          genes <- genes[genes != ""]
          length(unique(genes))
        }, integer(1))}
      
      format_geneid_html <- function(vec) {
        vapply(vec, function(txt) {
          if (is.na(txt) || txt == "") return("")
          genes <- unlist(strsplit(as.character(txt), "[/\\s]+"))
          genes <- genes[genes != ""]
          full_text <- paste(htmltools::htmlEscape(genes), collapse = ", ")
          display_text <- if (nchar(full_text) > 100) paste0(substr(full_text, 1, 100), "...") else full_text
          sprintf('<span style="display:none;">%s</span><span title="%s">%s</span>', full_text, full_text, display_text)
        }, character(1))}
      
      for (nm in list_names) {
        safe <- gsub("[^[:alnum:]_]", "_", nm)
        df <- ORA_list[[nm]]
        
        if (!("GeneID" %in% colnames(df))) df$GeneID <- ""
        if (!("GeneCount" %in% colnames(df))) df$GeneCount <- compute_gene_count(df$GeneID)
        if (!("p.adjust" %in% colnames(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
        
        if (is_personalized) {
          colnames(df)[colnames(df) == "GeneCount"] <- "FeatureCount"
          colnames(df)[colnames(df) == "GeneID"] <- "FeatureID"
          df$ID <- NULL }
        
        container_id <- paste0("resultsDotplotContainer_ind_", safe)        
        table_id <- paste0("resultsTable_ind_", safe)
        plot_id  <- paste0("resultsDotplot_ind_", safe)
        csv_id   <- paste0("downloadTableCSV_ind_", safe)
        tsv_id   <- paste0("downloadTableTSV_ind_", safe)
        html_id  <- paste0("downloadPlotHTML_ind_", safe)
        png_id   <- paste0("downloadPlotPNG_ind_", safe)
        jpg_id   <- paste0("downloadPlotJPG_ind_", safe)
        dropdown_ui_id <- paste0("downloadPlotUI_ind_", safe)
        
        local({
          df_local  <- df
          nm_local  <- nm
          tid <- table_id
          pid <- plot_id
          cid <- csv_id
          cid_ui <- container_id
          sid <- tsv_id
          hid <- html_id
          pid2 <- png_id
          jid <- jpg_id
          drop_ui <- dropdown_ui_id
          
          output[[tid]] <- DT::renderDataTable({
            req(df_local)
            allowed_cols <- intersect(input$selectedColumns, colnames(df_local))
            table_data <- df_local[, allowed_cols, drop = FALSE]
            
            if (is_personalized && "FeatureID" %in% colnames(table_data)) {
              table_data$FeatureID <- format_geneid_html(table_data$FeatureID)
            } else if ("GeneID" %in% colnames(table_data)) {
              table_data$GeneID <- format_geneid_html(table_data$GeneID)}
            
            DT::datatable(
              table_data,
              rownames = FALSE,
              selection = "multiple",
              filter = "top",
              escape = FALSE,
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                scrollY = paste0(input$resultsTable_length * 40, "px"),
                autoWidth = FALSE,
                columnDefs = list(list(
                  targets = "_all",
                  className = "dt-wrap",
                  width = paste0(round(100 / max(1, ncol(table_data))), "%")
                )),
                headerCallback = htmlwidgets::JS(
                  "function(thead, data, start, end, display) {",
                  "var tooltips =", jsonlite::toJSON(colTooltips), ";",
                  "$('th', thead).each(function() {",
                  "var txt = $(this).text().trim();",
                  "if (tooltips[txt]) {",
                  "var label = $('<span>').addClass('dt-header-label').text(txt);",
                  "var icon = $('<span>?</span>').addClass('dt-info-icon').attr('title', tooltips[txt]);",
                  "$(this).empty().append($('<div>').addClass('dt-header-with-icon').append(label).append(icon));",
                  "}",
                  "});",
                  "}"
                ))
            ) %>% formatSignif(columns = intersect(c("pvalue","p.adjust"),
                                                   colnames(table_data)), digits = 4)})
          
          filtered_data_ind <- reactive({
            rows <- input[[paste0(tid, "_rows_all")]]
            if (is.null(rows)) rows <- seq_len(nrow(df_local))
            df_local[rows, , drop = FALSE]})
          
          plot_obj_ind <- reactive({
            req(filtered_data_ind())
            finalResultsFiltered <- filtered_data_ind() %>%
              arrange(p.adjust) %>%
              head(input$nPlotTerms)
            
            text_size <- input$textSize
            
            if (is_personalized) {
              finalResultsFiltered$label_y <- finalResultsFiltered$Description
              p <- ggplot(finalResultsFiltered,
                          aes(x = FeatureCount,
                              y = reorder(label_y, -p.adjust),
                              text = paste0(
                                "Description: ", Description,
                                "<br>FeatureCount: ", FeatureCount,
                                "<br>p.adjust: ", signif(p.adjust, 3))))
            } else {
              finalResultsFiltered$label_y <- switch(input$plotLabel,
                                                     "ID" = finalResultsFiltered$ID,
                                                     "Description" = finalResultsFiltered$Description,
                                                     "Term_Description" =
                                                       paste(finalResultsFiltered$ID,
                                                             finalResultsFiltered$Description, sep = " - "))
              
              p <- ggplot(finalResultsFiltered,
                          aes(x = GeneCount,
                              y = reorder(label_y, -p.adjust),
                              text = paste0(
                                "ID: ", ID,
                                "<br>Description: ", Description,
                                "<br>GeneCount: ", GeneCount,
                                "<br>p.adjust: ", signif(p.adjust, 3))))}
            
            if (input$plotType == "dotplot") {
              if (is_personalized) {
                p <- p + aes(size = FeatureCount, color = p.adjust) +
                  labs(x = "FeatureCount", y = NULL, color = "p.adjust", title = paste("Dotplot -", nm_local))
              } else {
                p <- p + aes(size = GeneCount, color = p.adjust) +
                  labs(x = "GeneCount", y = NULL, color = "p.adjust", title = paste("Dotplot -", nm_local))}
              p <- p + geom_point() + scale_color_gradient(low = input$colorLow, high = input$colorHigh)
            } else { # Barplot
              if (is_personalized) {
                p <- p + aes(fill = p.adjust) +
                  labs(x = "FeatureCount", y = NULL, fill = "p.adjust", title = paste("Barplot -", nm_local))
              } else {
                p <- p + aes(fill = p.adjust) +
                  labs(x = "GeneCount", y = NULL, fill = "p.adjust", title = paste("Barplot -", nm_local))}
              p <- p + geom_bar(stat = "identity") + scale_fill_gradient(low = input$colorLow, high = input$colorHigh)}
            
            p + theme_minimal() +
              theme(axis.text.y  = element_text(size = text_size),
                    axis.text.x  = element_text(size = text_size),
                    axis.title.x = element_text(size = text_size + 1),
                    axis.title.y = element_text(size = text_size + 1),
                    legend.text  = element_text(size = text_size - 3),
                    legend.title = element_text(size = text_size, margin = margin(b = 20)),
                    plot.title   = element_text(size = text_size + 2,
                                                face = "bold", hjust = 0.5))})
          
          output[[cid_ui]] <- renderUI({
            req(plot_obj_ind())
            
            plot_height <- max(600, input$nPlotTerms * 25)
            shinycssloaders::withSpinner(
              plotlyOutput(ns(pid), height = paste0(plot_height, "px")))})
          
          output[[pid]] <- renderPlotly({
            req(plot_obj_ind())
            ggplotly(plot_obj_ind(), tooltip = "text")})
          
          output[[cid]] <- downloadHandler(
            filename = function()
              paste0(gsub("[^A-Za-z0-9\\-]", "_", nm_local), "_", Sys.Date(), ".csv"),
            content = function(file) {
              rows_all <- input[[paste0(tid, "_rows_all")]]
              if (is.null(rows_all)) rows_all <- seq_len(nrow(df_local))
              cols <- intersect(input$selectedColumns, colnames(df_local))
              if (length(cols) == 0) cols <- colnames(df_local)
              write.csv(df_local[as.numeric(rows_all), cols, drop = FALSE],
                        file, row.names = FALSE)})
          
          output[[sid]] <- downloadHandler(
            filename = function()
              paste0(gsub("[^A-Za-z0-9\\-]", "_", nm_local), "_", Sys.Date(), ".tsv"),
            content = function(file) {
              rows_all <- input[[paste0(tid, "_rows_all")]]
              if (is.null(rows_all)) rows_all <- seq_len(nrow(df_local))
              cols <- intersect(input$selectedColumns, colnames(df_local))
              if (length(cols) == 0) cols <- colnames(df_local)
              write.table(df_local[as.numeric(rows_all), cols, drop = FALSE],
                          file, row.names = FALSE, sep = "\t")})
          
          output[[drop_ui]] <- renderUI({
            shinyWidgets::dropdownButton(
              label = "Download Plot",
              icon = icon("download"),
              circle = FALSE,
              status = "secondary",
              width = "200px",
              tags$h5("Select format"),
              downloadButton(ns(hid), "HTML", class = "btn-block"),
              downloadButton(ns(pid2), "PNG", class = "btn-block"),
              downloadButton(ns(jid),  "JPG", class = "btn-block"))})
          
          output[[hid]] <- downloadHandler(
            filename = function()
              paste0("MetaEnrichGO_plot_", gsub("[^A-Za-z0-9\\-]", "_", nm_local),
                     "_", Sys.Date(), ".html"),
            content = function(file) {
              plot_h_px <- max(600, input$nPlotTerms * 25)
              plt <- ggplotly(plot_obj_ind(), 
                              height = plot_h_px, 
                              width = 1200, 
                              tooltip = "text")
              plt <- htmlwidgets::prependContent(plt, tags$style(HTML(paste0(
                ".plotly { height: ", plot_h_px, "px !important; width: 1200px !important; }"))))
              
              htmlwidgets::saveWidget(plt, file, selfcontained = TRUE)})
          output[[pid2]] <- downloadHandler(
            filename = function()
              paste0("MetaEnrichGO_plot_", gsub("[^A-Za-z0-9\\-]", "_", nm_local),
                     "_", Sys.Date(), ".png"),
            content = function(file) {
              plot_h_in <- max(8, input$nPlotTerms * 0.35)
              ggsave(file, plot = plot_obj_ind(),
                     device = "png", width = 12.5, height = plot_h_in,
                     units = "in", dpi = 96)})
          
          output[[jid]] <- downloadHandler(
            filename = function()
              paste0("MetaEnrichGO_plot_", gsub("[^A-Za-z0-9\\-]", "_", nm_local),
                     "_", Sys.Date(), ".jpg"),
            content = function(file) {
              plot_h_in <- max(8, input$nPlotTerms * 0.35)
              ggsave(file, plot = plot_obj_ind(),
                     device = "jpeg", width = 12.5, height = plot_h_in,
                     units = "in", dpi = 96)})})}})
    
    
    #-----------------------------------------------------------------------
    # Data representation (COMBINED) - mantenemos tus outputs previos
    #-----------------------------------------------------------------------
    output$resultsTable <- DT::renderDataTable({
      res <- resultsVal()$included
      req(res)
      req(input$selectedColumns)
      db_type <- selectedDatabase()
      is_personalized <- !is.null(db_type) && db_type == "Personalized"
      
      allowed_cols <- intersect(input$selectedColumns, colnames(res))
      table_data <- res[, allowed_cols, drop = FALSE]
      
      target_column <- NULL
      
      # Check for "GeneID" first, then "FeatureID"
      if ("GeneID" %in% colnames(table_data)) {
        target_column <- "GeneID"
      } else if ("FeatureID" %in% colnames(table_data)) {
        target_column <- "FeatureID"}
      
      if (!is.null(target_column)) {
        table_data[[target_column]] <- vapply(table_data[[target_column]], function(txt) {
          if (is.na(txt) || txt == "") return("")
          genes <- unlist(strsplit(as.character(txt), "[/\\s]+"))
          genes <- genes[genes != ""]
          esc_genes <- htmltools::htmlEscape(genes)
          full_text <- paste(esc_genes, collapse = ", ")
          display_text <- if (nchar(full_text) > 100) paste0(substr(full_text, 1, 100), "...") else full_text
          sprintf('<span style="display:none;">%s</span><span title="%s">%s</span>',
                  full_text, full_text, display_text)
        }, character(1))}
      
      if (is_personalized) {
        colTooltips <- list(
          Description = HTML("Descriptive name for the user-defined category."),
          pvalue = HTML("Combined p-value from ORA."),
          p.adjust = HTML("Adjusted p-value via Benjamini–Hochberg."),
          ListCount = HTML("Number of input lists that contributed features to this term."),
          ListNames = HTML("Names of the input lists containing features that contributed to this term"),
          FeatureCount = HTML("Total number of unique features from your input mapped to this term."),
          FeatureID = HTML("List of feature IDs from your input mapped to this term."))
      } else {
        colTooltips <- list(
          ID = HTML("Stable identifier of the enriched term or pathway."),
          Description = HTML("Descriptive name conveying the biological role of the term."),
          pvalue = HTML("Combined p-value from ORA."),
          p.adjust = HTML("Adjusted p-value via Benjamini–Hochberg."),
          ListCount = HTML("Number of input lists that contributed genes to this term."),
          ListNames = HTML("Names of the input lists containing genes that contributed to this term."),
          GeneCount = HTML("Total number of unique genes from your input mapped to this term."),
          GeneID = HTML("List of gene IDs from your input mapped to this term."))}
      
      num_rows <- input$resultsTable_length
      table_h  <- paste0(num_rows * 40, "px")
      
      DT::datatable(
        table_data,
        rownames  = FALSE,
        selection = "multiple",
        filter    = "top",
        escape    = FALSE,
        options   = list(
          pageLength = 10,
          scrollX    = TRUE,
          scrollY    = table_h,
          autoWidth  = FALSE,
          columnDefs = list(list(
            targets = "_all",
            className = "dt-wrap",
            width = paste0(round(100 / max(1, ncol(table_data))), "%")
          )),
          headerCallback = htmlwidgets::JS(
            "function(thead, data, start, end, display) {",
            "  var tooltips =", jsonlite::toJSON(colTooltips), ";",
            "  $('th', thead).each(function() {",
            "    var txt = $(this).text().trim();",
            "    if (tooltips[txt]) {",
            "      var label = $('<span>').addClass('dt-header-label').text(txt);",
            "      var icon  = $('<span>?</span>').addClass('dt-info-icon').attr('title', tooltips[txt]);",
            "      $(this).empty().append($('<div>').addClass('dt-header-with-icon').append(label).append(icon));",
            "    }",
            "  });",
            "}"))) %>% formatSignif(columns = intersect(c("pvalue", "p.adjust"), colnames(table_data)), digits = 4)})
    
    # Secondary (excluded) table
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
    
    output$downloadButtonUI <- renderUI({
      req(resultsVal()$included)
      fluidRow(
        column(5, downloadButton(ns("downloadTableCSV"), "Download CSV", style = "width:100%")),
        column(5, downloadButton(ns("downloadTableTSV"), "Download TSV", style = "width:100%")),
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
    
    output$downloadExcludedTermsList <- downloadHandler(
      filename = function() paste0("Excluded_MetaEnrichGO_", Sys.Date(), ".tsv"),
      content = function(file) {
        req(resultsVal())
        exc <- resultsVal()$excluded
        write.table(exc, file, sep = "\t", row.names = FALSE, quote = FALSE)})
    
    filtered_data <- reactive({
      req(resultsVal())
      if (!is.null(input$resultsTable_rows_all) &&
          length(input$resultsTable_rows_all) > 0) {
        resultsVal()$included[as.numeric(input$resultsTable_rows_all), ]
      } else {
        resultsVal()$included}})
    
    plot_obj <- reactive({
      req(filtered_data())
      finalResultsFiltered <- filtered_data() %>%
        arrange(p.adjust) %>%
        head(input$nPlotTerms)
      
      db_type <- selectedDatabase()
      is_personalized <- !is.null(db_type) && db_type == "Personalized"
      text_size <- input$textSize
      
      if (is_personalized) {
        finalResultsFiltered$label_y <- finalResultsFiltered$Description
        
        p <- ggplot(finalResultsFiltered,
                    aes(x = FeatureCount,
                        y = reorder(label_y, -p.adjust),
                        text = paste0(
                          "Description: ", Description,
                          "<br>FeatureCount: ", FeatureCount,
                          "<br>p.adjust: ", signif(p.adjust, 3))))
      } else {
        finalResultsFiltered$label_y <- switch(
          input$plotLabel,
          "ID"               = finalResultsFiltered$ID,
          "Description"      = finalResultsFiltered$Description,
          "Term_Description" = paste(finalResultsFiltered$ID, finalResultsFiltered$Description, sep = " - "))
        
        p <- ggplot(finalResultsFiltered,
                    aes(x = GeneCount,
                        y = reorder(label_y, -p.adjust),
                        text = paste0(
                          "ID: ", ID,
                          "<br>Description: ", Description,
                          "<br>GeneCount: ", GeneCount,
                          "<br>p.adjust: ", signif(p.adjust, 3))))}
      
      if (input$plotType == "dotplot") {
        if (is_personalized) {
          p <- p + aes(size = FeatureCount, color = p.adjust) +
            labs(x = "FeatureCount", y = NULL, color = "p.adjust", title = "Dotplot of the Enrichment Analysis")
        } else {
          p <- p + aes(size = GeneCount, color = p.adjust) +
            labs(x = "GeneCount", y = NULL, color = "p.adjust", title = "Dotplot of the Enrichment Analysis")}
        p <- p + geom_point() + scale_color_gradient(low = input$colorLow, high = input$colorHigh)
      } else { # Barplot
        if (is_personalized) {
          p <- p + aes(fill = p.adjust) +
            labs(x = "FeatureCount", y = NULL, fill = "p.adjust", title = "Barplot of the Enrichment Analysis")
        } else {
          p <- p + aes(fill = p.adjust) +
            labs(x = "GeneCount", y = NULL, fill = "p.adjust", title = "Barplot of the Enrichment Analysis")}
        p <- p + geom_bar(stat = "identity") + scale_fill_gradient(low = input$colorLow, high = input$colorHigh)}
      
      p <- p + theme_minimal() +
        theme(axis.text.y  = element_text(size = text_size),
              axis.text.x  = element_text(size = text_size),
              axis.title.x = element_text(size = text_size + 1),
              axis.title.y = element_text(size = text_size + 1),
              legend.text  = element_text(size = text_size - 3),
              legend.title = element_text(size = text_size, margin = margin(b = 20)),
              plot.title   = element_text(size = text_size + 2, face = "bold", hjust = 0.5))
      p})
    
    output$resultsDotplot <- renderPlotly({
      req(plot_obj())
      ggplotly(plot_obj(), tooltip = "text")})
    
    
    # Combined plot downloads UI
    output$downloadPlotUI <- renderUI({
      req(resultsVal()$included)
      shinyWidgets::dropdownButton(
        label = "Download Plot",
        icon = icon("download"),
        circle = FALSE,
        status = "secondary",
        width = "300px",
        tags$h4("Download Format"),
        downloadButton(ns("downloadPlotHTML"), "HTML", class = "btn-block"),
        downloadButton(ns("downloadPlotPNG"), "PNG",   class = "btn-block"),
        downloadButton(ns("downloadPlotJPG"), "JPG",   class = "btn-block"))})
    
    output$downloadPlotHTML <- downloadHandler(
      filename = function() paste0("MetaEnrichGO_plot_", Sys.Date(), ".html"),
      content = function(file) {
        plot_h_px <- max(600, input$nPlotTerms * 25)
        plt <- ggplotly(plot_obj(), 
                        height = plot_h_px, 
                        width = 1200, 
                        tooltip = "text")
        plt <- htmlwidgets::prependContent(plt, tags$style(HTML(paste0(
          ".plotly { height: ", plot_h_px, "px !important; width: 1200px !important; }"))))
        htmlwidgets::saveWidget(plt, file, selfcontained = TRUE)})
    
    output$downloadPlotPNG <- downloadHandler(
      filename = function() paste0("MetaEnrichGO_plot_", Sys.Date(), ".png"),
      content = function(file) {
        plot_h_in <- max(8, input$nPlotTerms * 0.35)
        ggsave(file, plot = plot_obj(), device = "png", 
               width = 12.5, height = plot_h_in, units = "in", dpi = 96)})
    
    output$downloadPlotJPG <- downloadHandler(
      filename = function() paste0("MetaEnrichGO_plot_", Sys.Date(), ".jpg"),
      content = function(file) {
        plot_h_in <- max(8, input$nPlotTerms * 0.35)
        ggsave(file, plot = plot_obj(), device = "jpeg", 
               width = 12.5, height = plot_h_in, units = "in", dpi = 96)})
    
    
    
    
    ###############################################################################
    
    #----------------------- NAVIGATION FUNCTIONS --------------------------------#
    
    ###############################################################################
    
    
    # Show the example of file that the app accepts (blue "i" button)
    observeEvent(input$showExample, {
      showModal(
        modalDialog(
          title = tags$div(icon("file-alt"), " Example File Contents", 
                           style = "text-align: center; font-weight: bold; font-size: 20px; color: #2980b9;"),
          easy_close = TRUE,
          size = "m",
          
          tags$div(
            style = "margin-bottom: 15px;",
            tags$h5(icon("code"), " Text Format (.txt)", style="color: #34495e; font-weight:bold; margin-top:0;"),
            HTML("<pre style='background-color:#f8f9fa; border: 1px solid #e9ecef; padding:12px; border-radius:5px; font-family: monospace; color: #2c3e50;'>COL10A1\nCOL11A1\nGREM1\nMMP1\nMMP12\nSPINK1\nSPP1\nSULF1\nTOX3\nTOP2A</pre>")),
          
          tags$div(
            style = "background-color: #eef6fb; padding: 15px; border-radius: 5px; border-left: 5px solid #2980b9;",
            tags$h5(icon("info-circle"), " Important Notes", style="color: #2980b9; font-weight:bold; margin-top:0;"),
            tags$ul(
              style = "margin-bottom: 0; padding-left: 20px;",
              tags$li(HTML("<b>Accepted format:</b> Only <code>.txt</code> files are supported for this input type."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>One gene per line:</b> Each row should contain a single gene name. No additional values allowed."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>No headers:</b> Do <u>not</u> include any header row. Input should start directly with gene names.")))),
          
          footer = modalButton("Close")))})
    
    
    # Show the example of file that the app accepts for annotations (blue "i" button)
    observeEvent(input$showPersonalizedInfo, {
      showModal(
        modalDialog(
          title = tags$div(icon("database"), " Personalized Annotations", 
                           style = "text-align: center; font-weight: bold; font-size: 20px; color: #2980b9;"),
          easy_close = TRUE,
          size = "l", 
          tags$div(
            style = "margin-bottom: 20px;",
            tags$h5(icon("table"), " Two-Column Format (Tab-Separated)", style="color: #34495e; font-weight:bold; margin-top:0;"),
            HTML('<pre style="background-color:#f8f9fa; border: 1px solid #e9ecef; padding:15px; border-radius:5px; overflow-x: auto; font-family: monospace; color: #2c3e50;">
FeatureID\tAnnotationValue
CYP2D7\tMembrane
CYP2D7\tMulti-pass membrane protein
CYP2D7\tCytoplasm
PIGBOS1\tMitochondrion outer membrane
PIGBOS1\tSingle-pass membrane protein
MT-RNR1\tNucleus
MT-RNR1\tMitochondrion
CD300H\tMembrane
CD300H\tSingle-pass type I membrane protein
NBDY\tCytoplasm</pre>')
          ),
          
          tags$div(
            style = "background-color: #eef6fb; padding: 15px; border-radius: 5px; border-left: 5px solid #2980b9; margin-bottom: 20px;",
            tags$h5(icon("clipboard-check"), " File Requirements", style="color: #2980b9; font-weight:bold; margin-top:0;"),
            tags$ul(
              style = "margin-bottom: 0; padding-left: 20px;",
              tags$li(HTML("<b>Single File, Two Columns:</b> Upload exactly one annotation file containing two columns."), style="margin-bottom: 6px;"),
              tags$li(HTML("<b>Required Header:</b> The file <b>must</b> include a header row (e.g., <code>FeatureID\\tAnnotationValue</code>)."), style="margin-bottom: 6px;"),
              tags$li(HTML("<b>Column Content:</b> 1st Column = Features (Genes/IDs), 2nd Column = Annotation (Location/Function)."), style="margin-bottom: 6px;"),
              tags$li(HTML("<b>Accepted Separator:</b> Tab-separated (<code>.tsv</code>, <code>.txt</code>) or Comma-separated (<code>.csv</code>).")))),
          
          tags$div(
            style = "background-color: #fff5f5; padding: 15px; border-radius: 5px; border-left: 5px solid #c0392b;",
            tags$h5(icon("exclamation-triangle"), " CRITICAL WARNING: Nomenclature", style="color: #c0392b; font-weight:bold; margin-top:0;"),
            tags$p(HTML("Ensure the feature names (IDs) in this annotation file <b>exactly match</b> the nomenclature (case-sensitivity, version) used in your input lists to avoid recognition issues."), 
                   style = "color: #7f2a21; margin-bottom: 0;")),
          
          footer = modalButton("Close")))})
    
    
    
    # 1) Reactive único que calcula el número total de listas de entrada
    totalLists <- reactive({
      req(input$inputMethod)
      if (isTRUE(input$useExampleData)) {
        4
      } else if (input$inputMethod == "files" && !is.null(input$files)) {
        length(input$files$datapath)
      } else if (input$inputMethod == "paste" &&
                 !is.null(input$pastedGenes) && nzchar(input$pastedGenes)) {
        txt <- input$pastedGenes
        if (grepl("###", txt, fixed = TRUE)) {
          length(strsplit(txt, "###", fixed = TRUE)[[1]])
        } else {
          1}} else {
            0}})
    
    # 2) Actualizar el slider usando totalLists()
    observe({n <- totalLists()
    if (n > 0) {
      isolate({
        updateSliderInput(session, "termThreshold",
                          min   = 1,
                          max   = n,
                          step  = 1,
                          value = n)})}})
    
    # 3) Mostrar/ocultar controles según el nº de listas
    observe({n <- totalLists()
    shinyjs::toggle(id = "metaAnalysisMethod", condition = n > 1)
    shinyjs::toggle(id = "termThreshold",     condition = n > 1)})
    
    
    # Descarga de archivo de anotación de ejemplo
    output$downloadPersonalizedExample <- downloadHandler(
      filename = function() {
        "protein_location_example.txt"},
      content = function(file) {
        src_path <- "./example_data/protein_location.txt"
        if (!file.exists(src_path)) {
          stop("Error: El archivo de ejemplo no se encontró en la ruta especificada: ", src_path)}
        file.copy(src_path, file)})
    
    
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
        shinyjs::reset("files")
        updateSwitchInput(session, "useExampleData", value = FALSE)}
      
      if (input$inputMethod == "files" && prev == "paste") {
        updateTextInput(session, "pastedGenes", value = "")}
      
      prev_input_method(input$inputMethod)})
    
    
    
    # ============================================================================== #
    #                                 VALIDATIONS                                    #
    # ============================================================================== #
    
    # ------------------------------------------------------------------------------
    #' @1.INPUT_EXISTENCE_CHECK
    # ------------------------------------------------------------------------------
    validate_input_exists <- function(input_data) {
      if (is.null(input_data) || 
          (input_data$type == "files" && length(input_data$paths) == 0) ||
          (input_data$type == "paste" && (!nzchar(input_data$text) || input_data$text == ""))) {
        
        showModal(modalDialog(
          title = tags$div(icon("exclamation-circle"), " No Input Data Detected", 
                           style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
          easyClose = TRUE, size = "m",
          
          tags$div(
            class = "alert alert-danger", role = "alert",
            tags$h4(icon("exclamation-triangle"), " No Data Detected"),
            tags$p("The analysis cannot be initiated because no gene lists have been provided.")),
          
          tags$div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #c0392b; margin-top: 15px;",
            tags$h5("Action Required:", style="font-weight:bold; margin-top:0;"),
            tags$ul(
              tags$li("Upload gene lists using ", tags$b(".txt files"), " or use the 'Paste Genes' input method."),
              tags$li(style="margin-top: 5px;", 
                      "For specifications regarding file structure, consult the ", 
                      tags$span(style="color: #17a2b8; font-weight:bold;", icon("info-circle"), " Information button"), 
                      " located next to the file upload control."))),
          
          footer = modalButton("Close")))
        return(FALSE)}
      return(TRUE)}
    
    
    
    # ------------------------------------------------------------------------------
    #' @2.DATABASE_FILE_CHECK
    # ------------------------------------------------------------------------------
    check_database_exists <- function(session, db_type, annot_file) {
      if (!is.null(db_type) && db_type == "Personalized" && is.null(annot_file)) {
        
        showModal(modalDialog(
          title = tags$div(icon("exclamation-circle"), " Annotation File Required", 
                           style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
          easyClose = TRUE, size = "m",
          
          tags$div(
            class = "alert alert-danger", role = "alert",
            tags$h4(icon("file-upload"), " No Annotation File"),
            tags$p("The analysis cannot be initiated because the 'Personalized' database requires a custom annotation file.")),
          
          tags$div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #c0392b; margin-top: 15px;",
            tags$h5("Action Required:", style="font-weight:bold; margin-top:0;"),
            tags$ul(
              tags$li("Upload your custom annotation file using the provided file input."),
              tags$li(style="margin-top: 5px;", 
                      "Ensure the file follows the required 2-column format (FeatureID and Description)."),
              tags$li(style="margin-top: 5px;",
                      "Consult the ", tags$span(style="color: #17a2b8; font-weight:bold;", icon("info-circle"), " Information button"), 
                      " for format examples."))),
          
          footer = modalButton("Close")))
        return(FALSE)}
      return(TRUE)}
    
    
    
    # ==============================================================================
    #' @3.MAPPING_FAILURE_CHECK
    # ==============================================================================
    check_mapping_failure <- function(ora_results, session) {
      
      found_valid_result <- FALSE
      
      for (i in seq_along(ora_results)) {
        res <- ora_results[[i]]
        
        if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
          first_col_values <- res[[1]]
          has_content <- any(!is.na(first_col_values) & as.character(first_col_values) != "")
          
          if (has_content) {
            found_valid_result <- TRUE}}}
      
      if (!found_valid_result) {
        showModal(modalDialog(
          title = tags$div(icon("exclamation-circle"), " Gene Mapping Failure", 
                           style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
          easyClose = TRUE, size = "l",
          tags$div(
            class = "alert alert-danger", role = "alert",
            tags$h4(icon("ban"), " Analysis Terminated"),
            tags$p("The enrichment analysis could not be performed because none of the provided identifiers mapped successfully to the selected database.")),
          tags$div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #c0392b; margin-top: 15px;",
            tags$h5("Diagnostic Report:", style="font-weight:bold; margin-top:0;"),
            tags$ul(
              tags$li("The analysis returned empty results for all input lists."),
              tags$li("If using 'Personalized', check that your Feature IDs match exactly."),
              tags$li("Check formatting: headers, hidden spaces, or incorrect delimiters."))),
          footer = modalButton("Close")))
        return(FALSE)}
      return(TRUE)}
    
    
    
    # ==============================================================================
    #' @4.FILE_FORMAT_ERROR
    # ==============================================================================
    validate_file_format_metaenrichgo <- function(input, session) {
      bad_files <- list()
      files <- if (input$inputMethod == "files") input$files$datapath else NULL
      file_names <- if (input$inputMethod == "files") input$files$name else NULL
      if (is.null(files)) return(TRUE)
      
      forbidden_checks <- list(
        "#"  = list(name = "Header/Comment (#)", msg = "Remove lines starting with '#'."),
        ","  = list(name = "Comma (,)",          msg = "File looks like CSV."),
        "|"  = list(name = "Pipes (|)",          msg = "Remove lines with '|'."), 
        "\t" = list(name = "Tab (\\t)",          msg = "Tabs detected (multiple columns?)."),
        ";"  = list(name = "Semicolon (;)",      msg = "Semicolons detected."))
      
      for (i in seq_along(files)) {
        lines <- readLines(files[i], warn = FALSE, encoding = "UTF-8")
        file_errors_html <- character() 
        
        for (char in names(forbidden_checks)) {
          bad_indices <- grep(char, lines, fixed = TRUE)
          
          if (length(bad_indices) > 0) {
            check_info <- forbidden_checks[[char]]
            
            preview_text <- head(lines[bad_indices], 3)
            preview_text <- ifelse(nchar(preview_text) > 50, paste0(substr(preview_text, 1, 50), "..."), preview_text)
            
            error_block <- paste0(
              "<div style='margin-bottom: 10px; border-bottom: 1px solid #eee; padding-bottom: 5px;'>",
              "<strong>", check_info$name, ":</strong> <span style='color:#c0392b;'>", check_info$msg, "</span>",
              "<div style='margin-top:4px; background-color:#fff3cd; padding:5px; border-radius:3px; font-family:monospace; font-size:0.85em; color:#856404;'>",
              paste0("Line ", head(bad_indices, 3), ": ", htmltools::htmlEscape(preview_text), collapse = "<br>"),
              if(length(bad_indices) > 3) paste0("<br><em>... (+", length(bad_indices)-3, " more lines)</em>") else "", 
              "</div></div>")
            
            file_errors_html <- c(file_errors_html, error_block)}}
        
        if (length(file_errors_html) > 0) {
          bad_files[[file_names[i]]] <- paste(file_errors_html, collapse = "")}}
      
      if (length(bad_files) > 0) {
        showModal(modalDialog(
          title = tags$div(icon("exclamation-circle"), " Format Error: File Input", 
                           style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
          size = "l", easyClose = TRUE,
          
          tags$div(class = "alert alert-danger", role = "alert",
                   tags$h4(icon("file-excel"), " Invalid File Format"),
                   tags$p("The analysis stopped because forbidden characters were detected.")),
          
          tags$div(style="max-height: 350px; overflow-y: auto; border: 1px solid #e2e2e2; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
                   lapply(names(bad_files), function(f) {
                     tags$div(style = "margin-bottom: 15px;",
                              tags$h5(icon("file-alt"), tags$b(f), style="color: #c0392b; margin-top:0; border-bottom: 2px solid #c0392b; padding-bottom:5px;"),
                              HTML(bad_files[[f]]))})),
          
          tags$div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #34495e;",
            tags$h5("Expected Format:", style="font-weight:bold; margin-top:0;"),
            tags$ul(tags$li("Plain text file (.txt) with one gene per line."),
                    tags$li("No headers, comments (#), commas, or tabs."))),
          
          footer = modalButton("Close and Fix Files")))
        return(FALSE)}
      return(TRUE)}
    
    # ==============================================================================
    #' @5.PASTE_FORMAT_ERROR 
    # ==============================================================================
    validate_paste_format <- function(input, session) {
      if (!(input$inputMethod == "paste" && nzchar(input$pastedGenes))) return(TRUE)
      
      txt <- input$pastedGenes
      lines <- unlist(strsplit(txt, "\n"))
      
      forbidden <- list(
        ","  = "Commas (,)", 
        "\t" = "Tabs (hidden)", 
        ";"  = "Semicolons (;)", 
        "|"  = "Pipes (|)", 
        ":"  = "Colons (:)",
        "#"  = "Hashtags/Headers (#)")
      
      all_errors_html <- character() 
      
      for (char in names(forbidden)) {
        found_issue <- FALSE
        bad_indices <- integer(0)
        
        if (char == "#") {
          candidates <- grep("#", lines, fixed = TRUE)
          real_errors <- candidates[!grepl("###", lines[candidates], fixed = TRUE)]
          
          if (length(real_errors) > 0) {
            found_issue <- TRUE
            bad_indices <- real_errors}
          
        } else {
          if (grepl(char, txt, fixed = TRUE)) { 
            found_issue <- TRUE
            bad_indices <- grep(char, lines, fixed = TRUE)}}
        
        if (found_issue) {
          error_name <- forbidden[[char]]
          
          preview_lines <- head(lines[bad_indices], 3)
          preview_html <- paste(htmltools::htmlEscape(preview_lines), collapse = "<br>")
          
          snippet <- paste0(
            "<div style='margin-bottom: 10px; border-bottom: 1px solid #eee; padding-bottom: 5px;'>",
            "<strong>Found:</strong> ", error_name,
            "<div style='margin-top:4px; background-color:#fff3cd; padding:5px; border-radius:3px; font-family:monospace; font-size:0.85em; color:#856404;'>",
            preview_html,
            if(length(bad_indices) > 3) paste0("<br><em>... (+", length(bad_indices)-3, " more lines)</em>") else "",
            "</div></div>")
          
          all_errors_html <- c(all_errors_html, snippet)}}
      
      if (length(all_errors_html) > 0) {
        
        showModal(modalDialog(
          title = tags$div(icon("exclamation-circle"), " Format Error: Paste Input", 
                           style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
          easyClose = TRUE, size = "l",
          
          tags$div(class = "alert alert-danger", role = "alert",
                   tags$h4(icon("times-circle"), " Analysis Stopped"),
                   tags$p("Multiple formatting issues detected. Please clean your input.")),
          
          tags$div(style="max-height: 350px; overflow-y: auto; border: 1px solid #e2e2e2; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
                   tags$h5(icon("keyboard"), tags$b("Issues Detected in Text"), style="color: #c0392b; margin-top:0; border-bottom: 2px solid #c0392b; padding-bottom:5px;"),
                   HTML(paste(all_errors_html, collapse=""))), 
          
          tags$div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #34495e;",
            tags$h5("Expected Format", style="margin-top:0; font-weight:bold;"),
            tags$ul(tags$li("Single column of gene symbols."),
                    tags$li("Separate multiple lists using ", tags$code("###"), "."),
                    tags$li("No headers, hidden tabs, or extra columns."),
                    tags$li("Tip: Paste into Notepad first to clean hidden formatting."))),
          
          footer = modalButton("Close and Fix")))
        
        shinyjs::disable("runButton"); shinyjs::delay(1000, shinyjs::enable("runButton"))
        return(FALSE)}
      return(TRUE)}
    
    # ------------------------------------------------------------------------------
    #' @6.ORGANISM_FORMAT_CHECK
    # ------------------------------------------------------------------------------
    validate_organism <- function(gene_list, organism, ID_type, session) {
      if (length(gene_list) == 0) return(gene_list)
      
      pattern <- switch(ID_type,
                        SYMBOL = switch(organism, Hsa="^[A-Z][A-Z0-9._/\\-]*$", Mmu="^[A-Z][a-z0-9._/\\-]*$", Rno="^[A-Z][a-z0-9._/\\-]*$", NULL),
                        ENTREZID = "^[0-9]+$",
                        ENSEMBL = switch(organism, Hsa="^ENSG[0-9]{11}$", Mmu="^ENSMUSG[0-9]{11}$", Rno="^ENSRNOG[0-9]{11}$", NULL), NULL)
      
      if (is.null(pattern)) return(gene_list)
      
      is_valid <- grepl(pattern, gene_list)
      valid_genes <- gene_list[is_valid]; invalid_genes <- gene_list[!is_valid]
      percent_invalid <- (length(invalid_genes) / length(gene_list)) * 100
      org_name <- switch(organism, Hsa="Homo sapiens", Mmu="Mus musculus", Rno="Rattus norvegicus")
      
      if (length(invalid_genes) == 0) return(gene_list)
      
      example_box <- tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #34495e; margin-top: 15px;",
        tags$h5("Correct Nomenclature:", style = "margin-top:0; font-weight:bold;"),
        tags$ul(tags$li(tags$b("Human:"), " All Uppercase (TP53)"), tags$li(tags$b("Mouse/Rat:"), " Title Case (Tp53)")))
      
      if (percent_invalid >= 50) {
        showModal(modalDialog(
          title = tags$div(icon("ban"), " Analysis Stopped: Data Mismatch", style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
          easyClose = TRUE, size = "l",
          tags$div(class = "alert alert-danger", role = "alert",
                   tags$h4(icon("exclamation-triangle"), paste0(round(percent_invalid, 1), "% invalid genes")),
                   tags$p(paste0("Most genes do not match the nomenclature for ", org_name, "."))),
          tags$div(style="margin:15px 0; padding:10px; border:1px solid #f5c6cb; background:#f8d7da; border-radius:5px;",
                   tags$b("Invalid Examples: "), paste(head(unique(invalid_genes), 5), collapse=", ")),
          example_box, footer = modalButton("Close and Fix Data")))
        return(NULL)} 
      
      showModal(modalDialog(
        title = tags$div(icon("exclamation-circle"), " Nomenclature Warning", style = "color: #e67e22; font-weight: bold; font-size: 18px;"),
        easyClose = TRUE, size = "l",
        tags$div(class = "alert alert-warning", role = "alert",
                 tags$h4(icon("info-circle"), paste(length(invalid_genes), "genes will be excluded")),
                 tags$p(paste0("The analysis will continue, but some genes do not match ", org_name, "."))),
        tags$div(style="max-height: 100px; overflow-y: auto; background: #fff3cd; padding: 10px; border: 1px solid #ffeeba; border-radius: 4px;",
                 tags$code(paste(head(unique(invalid_genes), 20), collapse=", "))),
        example_box, footer = modalButton("I understand, continue")))
      
      return(valid_genes)}
    
    # ------------------------------------------------------------------------------
    #' @7.GENE_ID_FORMAT_CHECK
    # ------------------------------------------------------------------------------
    validate_gene_ids <- function(gene_vector, ID_type, session) {
      if (length(gene_vector) == 0) return(gene_vector)
      gene_list <- as.character(gene_vector)
      ensembl_pattern <- "^ENS[A-Z]{0,4}[0-9]{6,}$"
      pattern <- switch(ID_type, SYMBOL="^[A-Za-z][A-Za-z0-9._/\\-]*$", ENTREZID="^[0-9]+$", ENSEMBL=ensembl_pattern, stop("Error"))
      
      is_valid <- grepl(pattern, gene_list)
      if (ID_type == "SYMBOL") is_valid <- is_valid & !grepl(ensembl_pattern, gene_list)
      
      valid_genes <- gene_list[is_valid]; invalid_genes <- gene_list[!is_valid]
      percent_invalid <- (length(invalid_genes) / length(gene_list)) * 100
      
      if (length(invalid_genes) == 0) return(gene_list)
      
      example_box <- tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #34495e; margin-top: 15px;",
        tags$h5("Expected ID Formats:", style = "margin-top:0; font-weight:bold;"),
        tags$ul(tags$li(tags$b("SYMBOL:"), " e.g., TP53, Brca1"), 
                tags$li(tags$b("ENTREZID:"), " Digits only (e.g., 7157)"),
                tags$li(tags$b("ENSEMBL:"), " Starts with 'ENS'")))
      
      if (percent_invalid >= 50) {
        showModal(modalDialog(
          title = tags$div(icon("ban"), " Analysis Stopped: ID Mismatch", style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
          easyClose = TRUE, size = "l",
          tags$div(class = "alert alert-danger", role = "alert",
                   tags$h4(icon("times-circle"), paste0(round(percent_invalid, 1), "% invalid IDs")),
                   tags$p(paste0("Data does not match the selected ID type: ", ID_type, "."))),
          tags$div(style="margin:15px 0; padding:10px; border:1px solid #f5c6cb; background:#f8d7da; border-radius:5px;",
                   tags$b("Invalid Examples: "), paste(head(unique(invalid_genes), 5), collapse=", ")),
          example_box, footer = modalButton("Close and Fix Data")))
        return(NULL)}
      
      showModal(modalDialog(
        title = tags$div(icon("exclamation-circle"), " ID Format Warning", style = "color: #e67e22; font-weight: bold; font-size: 18px;"),
        easyClose = TRUE, size = "l",
        tags$div(class = "alert alert-warning", role = "alert",
                 tags$h4(icon("info-circle"), paste(length(invalid_genes), "IDs will be excluded")),
                 tags$p("The analysis will continue with valid IDs only.")),
        tags$div(style="max-height: 100px; overflow-y: auto; background: #fff3cd; padding: 10px; border: 1px solid #ffeeba; border-radius: 4px;",
                 tags$code(paste(head(unique(invalid_genes), 20), collapse=", "))),
        example_box, footer = modalButton("I understand, continue")))
      
      return(valid_genes)}
    
    # ==============================================================================
    #' @8.ANNOTATION_FILE_FORMAT_CHECK
    # ==============================================================================
    check_annotation_file_format <- function(file_path, session) {
      
      show_error_modal <- function(title_text, error_msg, actual_preview = NULL) {
        preview_html <- ""
        if (!is.null(actual_preview)) {
          preview_text <- paste(as.character(actual_preview), collapse="\n")
          preview_html <- paste0(
            "<div style='margin-top:5px; background-color:#fff3cd; padding:5px; border-radius:3px; font-family:monospace; font-size:0.9em; color:#856404;'>",
            "<strong>File Preview (First 3 lines):</strong><br>",
            htmltools::htmlEscape(preview_text),
            "</div>")}
        
        showModal(modalDialog(
          title = tags$div(icon("table"), " Annotation File Error", 
                           style = "text-align: center; font-weight: bold; font-size: 20px; color: #c0392b;"),
          size = "l", easyClose = TRUE,
          tags$div(class = "alert alert-danger", role = "alert",
                   tags$h4(icon("exclamation-triangle"), " Format Mismatch"),
                   tags$p(HTML(error_msg))),
          tags$div(style="max-height: 300px; overflow-y: auto; border: 1px solid #e2e2e2; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
                   tags$div(style = "margin-bottom: 10px;",
                            tags$h5(icon("file-alt"), tags$b("Annotation File"), style="color: #c0392b; margin-top:0;"),
                            HTML(preview_html))),
          tags$div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #34495e;",
            tags$h5("Expected Format", style="margin-top:0; font-weight:bold;"),
            tags$ul(tags$li("Exactly ", tags$b("two columns"), "."),
                    tags$li("Headers are required (e.g., FeatureID, Description)."),
                    tags$li("If headers are missing, the app creates them automatically, but format must be consistent."))),
          footer = modalButton("Close and Fix File")))}
      
      tryCatch({
        raw_lines <- readLines(file_path, warn = FALSE)
        if (length(raw_lines) == 0) {
          show_error_modal("Read Error", "The file appears to be empty.", NULL)
          return(FALSE)}
        writeLines(raw_lines, file_path)
        sep <- if(any(grepl("\t", raw_lines))) "\t" else ","
        df <- utils::read.table(file_path, header = TRUE, sep = sep, nrows = 5, 
                                fill = TRUE, stringsAsFactors = FALSE, 
                                quote = "\"", comment.char = "")
        
        if (ncol(df) != 2) {
          msg <- paste0("The file appears to have <b>", ncol(df), " columns</b>. The tool strictly requires <b>2 columns</b>.")
          show_error_modal("Format Error", msg, head(raw_lines, 3))
          return(FALSE)}
        
        full_df <- utils::read.table(file_path, header = TRUE, sep = sep, 
                                     fill = TRUE, stringsAsFactors = FALSE, 
                                     quote = "\"", comment.char = "")
        
        colnames(full_df) <- c("FeatureID", "Description")
        utils::write.table(full_df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
        return(TRUE)
        
      }, error = function(e) {
        safe_preview <- tryCatch(readLines(file_path, n=3, warn=FALSE), error=function(x) "Unreadable content")
        show_error_modal("Read Error", paste0("Critical error reading file: ", e$message), safe_preview)
        return(FALSE)})}
    
    
    # ==============================================================================
    #' @9.STRICT_THRESHOLD_CHECK
    # ==============================================================================
    check_strict_threshold <- function(meta_results, term_frequency, threshold, session) {
      
      if (nrow(meta_results) > 0) return(TRUE)
      if (length(term_frequency) > 0) {
        
        max_freq_found <- max(term_frequency)
        total_potential <- length(term_frequency)
        
        showModal(modalDialog(
          title = tags$div(icon("filter"), " No Significant Results Found", 
                           style = "color: #e67e22; font-weight: bold; font-size: 18px;"),
          easyClose = TRUE, size = "m",
          
          tags$div(
            class = "alert alert-warning", role = "alert",
            tags$h4(icon("search-minus"), " Threshold Too Strict"),
            tags$p(HTML(paste0("The analysis generated <b>", total_potential, "</b> enriched terms, but none appeared in <b>", threshold, "</b> lists.")))),
          
          tags$div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #e67e22; margin-top: 15px;",
            tags$h5("Analysis Stopped:", style="font-weight:bold; margin-top:0;"),
            tags$ul(
              tags$li(paste("Current Threshold:", threshold)),
              tags$li(paste("Max Overlap Found:", max_freq_found)),
              tags$li("No table will be displayed."))),
          
          tags$p(style="margin-top:15px; font-weight: bold; text-align: center; color: #333;", 
                 paste("Suggestion: Set 'Minimum number of datasets' to", max_freq_found, "or lower.")),
          
          footer = modalButton("Close")))
        
        return(FALSE)}
      return(TRUE)}
    
    
    
    
    
  })}