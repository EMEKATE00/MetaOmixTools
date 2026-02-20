################################################################################
#                                                                              #
#                               Master’s Thesis                                #
#                                                                              #
#     Title: TFM: Development of the Meta-Analysis Tool MetaRank               #
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
library(zip)                    # makes zipped downloadable files 

#------------------------------------------------------------------------------- 
# The user interface (ui) object controls the layout and appearance of your app
#-------------------------------------------------------------------------------

source("./modules/metarank/R/MetaRank_Functions.R")
source("./modules/metarank/R/ORA.R")

mod_metaRank_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$head(
      includeCSS("./modules/metarank/www/style.css")),
    
    page_sidebar(
      id    = ns("metarank_container"),
      title = "MetaRank",
      
      sidebar = sidebar(width = "380px", useShinyjs(),
                        
                        # Select the type of analysis: 
                        div(selectInput(
                          inputId = ns("selectedPackage"),
                          label   = "Select Meta-Analysis type",
                          choices = c("Weighted (ranked by scores)"       = "Rankprodpack",
                                      "Unweighted (simple ranked lists)"   = "RRApack"),
                          selected = "Rankprodpack"),
                          title = HTML(
                            "Choose between:\n- Weighted: Gene lists are ranked based on numerical scores that represent their importance or association with a given condition. These scores could be p-values, fold-changes, or other statistical measures that quantify the relevance of each gene.\n- Unweighted: Gene lists are simply ranked by their order in different datasets, without assigning any numerical score or weight. This method assumes that the genes’ relative positions across datasets are the key factor, rather than a statistical score.\n- The Weighted method is more suitable when you have numerical data that reflects the significance of each gene. The Unweighted method is ideal for cases where you are only interested in the relative rankings of genes across different datasets."),
                          style = "cursor: help; display: inline-block; width: 100%;"),
                        
                        tabsetPanel(
                          id = ns("tabs"),
                          
                          # PESTAÑA 1
                          tabPanel(title = "Analysis",
                                   
                                   #' @DATA_INPUT_METHODS
                                   tags$h5("Input", style = "margin-top: 20px;"),
                                   
                                   radioButtons(
                                     inputId = ns("inputMethod"),
                                     label   = "Input Method",
                                     choices = c("Upload Files" = "files", "Paste Genes" = "paste"),
                                     selected = "files"),
                                   
                                   # Conditional: Files + Rankprodpack
                                   #--------------------------------------------------------------
                                   conditionalPanel(
                                     condition = sprintf(
                                       "input['%s'] == 'files' && input['%s'] == 'Rankprodpack'",
                                       ns("inputMethod"), ns("selectedPackage")),
                                     fluidRow(
                                       column(10, fileInput(inputId   = ns("files"),
                                                            label     = "Upload gene rankings (RP)",
                                                            multiple  = TRUE,
                                                            accept    = c(".txt", ".tsv", ".csv"))),
                                       column(2, actionButton(inputId = ns("showExample"), label = "", 
                                                              icon     = icon("info-circle", lib = "font-awesome"),
                                                              class    = "btn-info"))),
                                     
                                     tags$p("Use example data?" , style = "margin-bottom: -10px; margin-top: -15px;"),
                                     fluidRow(
                                       column(10, div(switchInput(
                                         inputId  = ns("useExampleData"),
                                         onLabel  = "Yes",
                                         offLabel = "No",
                                         value    = FALSE),
                                         title = HTML(
                                           "These example human gene lists (4 studies related to lung cancer: GSE10072, GSE19188, GSE63459, GSE75037) contain gene symbols and p-values, have no headers, and include duplicates and blank lines to simulate raw data."),
                                         style = "cursor: help;")),
                                       column(2, downloadLink(
                                         outputId = ns("downloadExampleRank"),
                                         label    = icon("download"),
                                         style    = "margin-top: 10px;")))),
                                   
                                   
                                   # Condicional: Files + RRApack
                                   #--------------------------------------------------------------
                                   conditionalPanel(
                                     condition = sprintf(
                                       "input['%s'] == 'files' && input['%s'] == 'RRApack'",
                                       ns("inputMethod"), ns("selectedPackage")),
                                     fluidRow(
                                       column(10, fileInput(inputId  = ns("files_rra"),
                                                            label    = "Upload gene rankings (RRA)",
                                                            multiple = TRUE,
                                                            accept   = c(".txt"))),
                                       column(2, actionButton(inputId = ns("showExampleRRA"), label = "", 
                                                              icon    = icon("info-circle", lib = "font-awesome"),
                                                              class   = "btn-info"))),
                                     
                                     tags$p("Use example data?", style = "margin-top: -15px; margin-bottom: -10px;"),
                                     fluidRow(
                                       column(10, div(switchInput(
                                         inputId  = ns("useExampleData_rra"),
                                         onLabel  = "Yes",
                                         offLabel = "No",
                                         value    = FALSE),
                                         title = HTML(
                                           "These example human gene lists (4 studies related to lung cancer: GSE10072, GSE19188, GSE63459, GSE75037) contain gene symbols only (no p-values), have no headers, include duplicates and blank lines to simulate raw data."),
                                         style = "cursor: help;")),
                                       column(2, downloadLink(
                                         outputId = ns("downloadExampleRRA"),
                                         label    = icon("download"),
                                         style    = "margin-top: 10px;")))),
                                   
                                   
                                   # Condicional: Paste + Rankprodpack
                                   #--------------------------------------------------------------
                                   conditionalPanel(
                                     condition = sprintf(
                                       "input['%s'] == 'paste' && input['%s'] == 'Rankprodpack'",
                                       ns("inputMethod"), ns("selectedPackage")),
                                     radioButtons(
                                       inputId = ns("pasteFormat"),
                                       label   = "Select Paste Format",
                                       choices = c("CSV" = "csv", "TSV" = "tsv"),
                                       selected = "csv"),
                                     
                                     conditionalPanel(
                                       condition = sprintf("input['%s'] == 'csv'", ns("pasteFormat")),
                                       div(
                                         textAreaInput(
                                           inputId     = ns("pasteText"),
                                           label       = "Paste your data",
                                           width       = "300px",
                                           height      = "200px",
                                           placeholder = "Myc,1.2e-06\nGapdh,3.5e-04\nIl2ra,0.0021\n...\n###\nGapdh,5.4e-07\nIl6,9.8e-04\nNanog,0.0062\n..."),
                                         title = HTML("CSV Paste Format: \n- Each block must start with Gene,Stat.data header. \n- Use comma ',' as field separator and '.' as decimal point.\n- Separate multiple lists with '###' on its own line.\n- Do not use commas as decimal separators.\n- For full details, see the format instructions."),
                                         style = "cursor: help;")),
                                     
                                     conditionalPanel(
                                       condition = sprintf("input['%s'] == 'tsv'", ns("pasteFormat")),
                                       div(
                                         textAreaInput(
                                           inputId     = ns("pasteText1"),
                                           label       = "Paste your data",
                                           width       = "300px",
                                           height      = "200px",
                                           placeholder = "Myc\t1.2e-06\nGapdh\t3.5e-04\nIl2ra\t0.0021\n...\n###\nGapdh\t5.4e-07\nIl6\t9.8e-04\nNanog\t0.0062\n..."),
                                         title = HTML(
                                           "TSV Paste Format:\n- Each block must start with 'Gene\\tStat.data' header.\n- Use tabulation '\\t' as field separator and '.' as decimal point.\n- Separate multiple lists with '###' on its own line.\n- Do not use commas or semicolons as separators.\n- For full details, see the format instructions."),
                                         style = "cursor: help;"))),
                                   
                                   
                                   # Condicional: Paste + RRApack
                                   #--------------------------------------------------------------
                                   conditionalPanel(
                                     condition = sprintf(
                                       "input['%s'] == 'paste' && input['%s'] == 'RRApack'",
                                       ns("inputMethod"), ns("selectedPackage")),
                                     div(
                                       textAreaInput(
                                         inputId     = ns("pasteText_rra"),
                                         label       = "Paste your ranked gene lists",
                                         width       = "300px",
                                         height      = "200px",
                                         placeholder = "TP53\nMYC\nEGFR\nBRCA1\n...\n###\nBRCA1\nPTEN\nMYC\n..."),
                                       title = HTML(
                                         "RRA Paste Format:\n- No header required; only gene symbols.\n- One gene per line, separated by '\\n'.\n- Separate multiple lists with '###' on its own line.\n- Do not include numeric values or additional columns.\n- Raw, unweighted gene rankings only."),
                                       style = "cursor: help;")),
                                   
                                   tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
                                   
                                   
                                   
                                   #' @PARAMETERS
                                   tags$h5("Parameters"),
                                   
                                   # -----------------------------
                                   # Conditional: Rankprodpack
                                   # -----------------------------
                                   conditionalPanel(
                                     condition = sprintf("input['%s'] == 'Rankprodpack'", ns("selectedPackage")),
                                     
                                     # Meta-Analysis method
                                     div(style = "margin-top: -5px; margin-bottom: -5px;",
                                         selectInput(
                                           inputId = ns("RankBasedMethod"),
                                           label   = "Rank-Based Meta-Analysis Methods",
                                           choices = c("RankProd basic" = "rp", "RankProd advance" = "rp.adv"),
                                           selected = "rp"),
                                         title = HTML(
                                           "Choose your RankProd method:\n- RankProd basic: Basic analysis without additional grouping.\n- RankProd advance: Advanced RankProd allowing an 'origin' vector to group your lists by metadata (e.g., data source, date of acquisition, technology)."),
                                         style = "cursor: help; display: inline-block;"),
                                     
                                     # Origin
                                     conditionalPanel(
                                       condition = sprintf("input['%s'] == 'rp.adv'", ns("RankBasedMethod")),
                                       fluidRow(
                                         column(10, textInput(
                                           inputId     = ns("origin"),
                                           label       = "Origin",
                                           placeholder = "1,1,2,2,1")),
                                         column(2,  actionButton(
                                           inputId     = ns("showOriginInfo"),
                                           label       = "",
                                           icon        = icon("info-circle", lib = "font-awesome"),
                                           class       = "btn-info")))),
                                     
                                     # Minimum treshold
                                     div(style = "margin-top: -5px; margin-bottom: -5px;",
                                         sliderInput(
                                           inputId = ns("gene_app_filter_rp"),
                                           label   = "Minimum number of datasets",
                                           min     = 1,
                                           max     = 4,
                                           value   = 4,
                                           step    = 1,
                                           ticks   = TRUE),
                                         title = HTML(
                                           "Gene Appearance Filter:\n - Require a gene to appear in at least this many lists.\n- E.g., if BRCA1 appears in 2/4 lists and you set this to 3, it will be excluded.\n- Setting to 2 will keep it."),
                                         style = "cursor: help;"),
                                     
                                     # Ranking direction
                                     div(
                                       style = "margin-top: -5px; margin-bottom: -5px;",
                                       radioButtons(
                                         inputId  = ns("ranking_direction"),
                                         label    = "Ranking Direction",
                                         choices  = c(
                                           "Ascending (lower values = better)"  = "asc",
                                           "Descending (higher values = better)" = "desc"),
                                         selected = "asc",
                                         inline   = TRUE),
                                       title = HTML("..."),
                                       style = "cursor: help;"),
                                     
                                     # NA management
                                     div(style = "margin-top: -5px; margin-bottom: -5px;",
                                         selectInput(
                                           inputId  = ns("NA_management"),
                                           label    = "NA Management",
                                           choices  = c(
                                             "Impute NA"   = "impute",
                                             "Ignore NA"   = "ignore",
                                             "Penalize NA" = "penalize"),
                                           selected = "impute"),
                                         title = HTML(
                                           "Missing Value Handling:\n- Impute NA: Replace missing with the gene’s median across lists (extra penalty optional).\n- Ignore NA: Analyze only present values (extra penalty optional).\n- Penalize NA: Replace missing with the worst value (no extra penalty)."),
                                         style = "cursor: help;"),
                                     
                                     # Extra Penalization
                                     conditionalPanel(
                                       condition = sprintf(
                                         "input['%s'] == 'impute' || input['%s'] == 'ignore'",
                                         ns("NA_management"), ns("NA_management")),
                                       tags$p("Apply an extra penalization?", style = "margin-bottom: -10px; margin-top: -5px;"),
                                       div(switchInput(
                                         inputId  = ns("extra_penalty_switch"),
                                         onLabel  = "Yes",
                                         offLabel = "No",
                                         value    = FALSE),
                                         title = HTML(
                                           "Extra Penalization:\n- After initial ranking, adds a penalty based on missing appearances.\n- AdjustedRank = Rank + ((total_lists − Count) × (max_rank / total_lists)).\n- Genes in fewer lists receive a larger penalty."),
                                         style = "cursor: help;"))),
                                   
                                   # -----------------------------
                                   # Conditional: RRApack
                                   # -----------------------------
                                   conditionalPanel(
                                     condition = sprintf("input['%s'] == 'RRApack'", ns("selectedPackage")),
                                     
                                     # Meta-Analysis Method
                                     div(selectInput(
                                       inputId = ns("rra_method"),
                                       label   = "Aggregation Method",
                                       choices = c(
                                         "Robust Rank Aggregation" = "RRA",
                                         "Min"                      = "min",
                                         "Geometric Mean"           = "geom.mean",
                                         "Arithmetic Mean"          = "mean",
                                         "Median"                   = "median",
                                         "Stuart"                   = "stuart"),
                                       selected = "RRA"),
                                       title = HTML(
                                         "Aggregation Methods:\n- RRA: Probabilistic rank aggregation robust to noise.\n- Min: Takes the minimum rank across lists.\n- Geometric Mean: Geometric mean of ranks.\n- Arithmetic Mean: Arithmetic mean of ranks.\n- Median: Median of ranks.\n- Stuart: Stuart’s method for order-statistics aggregation."),
                                       style = "cursor: help;"),
                                     
                                     # Minimum Treshold
                                     div(sliderInput(
                                       inputId = ns("gene_app_filter_rra"),
                                       label   = "Minimum number of datasets",
                                       min     = 1,
                                       max     = 4,
                                       value   = 4,
                                       step    = 1,
                                       ticks   = TRUE),
                                       title = HTML(
                                         "Gene Appearance Filter:\n- Require a gene to appear in at least this many lists.\n- E.g., setting to 2 keeps genes in ≥2/4 lists."),
                                       style = "cursor: help;"),
                                     
                                     # NA managemnet
                                     div(selectInput(
                                       inputId = ns("rra_na_handling"),
                                       label   = "Missing Value Handling",
                                       choices = c(
                                         "Ignore missing values"   = "ignore",
                                         "Penalize missing values" = "penalize"),
                                       selected = "ignore"),
                                       title = HTML(
                                         "Missing Value Handling:\n- Ignore: Proceed with only present values.\n- Penalize: Substitute missing with the worst observed rank."),
                                       style = "cursor: help;")),
                                   
                                   tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
                                   
                                   #--------------------------------------------------------------
                                   # Analysis options
                                   #--------------------------------------------------------------
                                   
                                   #' @ENRICHMENT_ANALYSIS
                                   tags$h5("Enrichment analysis"),
                                   
                                   # Selector y/n
                                   div(checkboxInput(inputId = ns("enableEnrichment"),
                                                     label = "Perform Enrichment Analysis?",
                                                     value = FALSE),
                                       title = HTML("Check this box to run an Over-Representation Analysis (ORA) on your meta-analysis results. It will take the top N genes you’ve selected (up to 100) from the consensus ranking and test for pathway or GO term enrichment."),
                                       style = "cursor: help;"),
                                   
                                   conditionalPanel(
                                     condition = sprintf("input['%s'] == true", ns("enableEnrichment")),
                                     
                                     # Gene number used in the analysis
                                     div(fluidRow(column(width = 8,
                                                         sliderInput(inputId = ns("topGenes"),
                                                                     label   = "Number of Top Genes",
                                                                     min     = 1,
                                                                     max     = 100,
                                                                     value   = 30,
                                                                     step    = 1)),
                                                  column(width = 4,  
                                                         numericInput(inputId = ns("topGenes_numeric"),
                                                                      label   = NULL,
                                                                      value   = 30,
                                                                      min     = 1,
                                                                      max     = 100,
                                                                      step    = 1)),
                                                  
                                                  title = HTML("Select how many of the top-ranked genes from the meta-analysis consensus list you want to use for the over-representation analysis (ORA)."),
                                                  style = "cursor: help;")),
                                     
                                     # Database
                                     selectInput(inputId = ns("oraDatabase"),
                                                 label = "ORA Database",
                                                 choices = c("Gene Ontology (GO)" = "go", "KEGG" = "kegg", "Reactome" = "reactome", "Personalized" = "Personalized"),
                                                 selected = "go"),
                                     
                                     conditionalPanel(
                                       condition = sprintf("input['%s'] == 'go'", ns("oraDatabase")),
                                       
                                       # GO Options
                                       radioButtons(inputId = ns("ontology"),
                                                    label = "Ontology",
                                                    choices = c("Biological Process (BP)" = "bp",
                                                                "Molecular Function (MF)" = "mf",
                                                                "Cellular Component (CC)" = "cc"),
                                                    selected = "bp")),
                                     
                                     # Personalized DB
                                     conditionalPanel(
                                       condition = sprintf("input['%s'] == 'Personalized'", ns("oraDatabase")),
                                       fluidRow(
                                         column(width = 10,
                                                fileInput(
                                                  inputId = ns("personalizedAnnotations"),  
                                                  label = "Choose annotation file",
                                                  multiple = FALSE,              
                                                  accept = c(".txt", ".tsv", ".csv"))),
                                         column(width = 2,
                                                actionButton(
                                                  inputId = ns("showPersonalizedInfo"),      
                                                  label = "",
                                                  icon = icon("info-circle", lib = "font-awesome"),
                                                  class = "btn-info"),
                                                br(), 
                                                tags$div(
                                                  downloadLink(
                                                    outputId = ns("downloadPersonalizedExample"),
                                                    label = icon("download"),
                                                    style = "margin-left: 4px; font-size: 16px;"),
                                                  style = "margin-top: 10px;")))),
                                     
                                     # Organism
                                     conditionalPanel(
                                       condition = sprintf("input['%s'] != 'Personalized'", ns("oraDatabase")),
                                       selectInput(inputId = ns("species"),
                                                   label = "Species",
                                                   choices = c("Homo sapiens" = "Hsa",
                                                               "Mus musculus" = "Mmu",
                                                               "Rattus norvegicus" = "Rno"),
                                                   selected = "Hsa"),
                                       # ID Type
                                       selectInput(inputId = ns("IDtype"),
                                                   label = "Gene ID",
                                                   choices = c("SYMBOL", "ENTREZID", "ENSEMBL"),
                                                   selected = "SYMBOL")))),
                          
                          
                          tabPanel(
                            
                            # PESTAÑA 2
                            title = "Data Visualization", br(),
                            
                            #--------------------------------------------------------------------
                            # Meta-analysis table setting
                            #--------------------------------------------------------------------
                            tags$h5("Meta-analysis table setting"),
                            conditionalPanel(
                              condition = sprintf("input['%s'] == 'Rankprodpack'", ns("selectedPackage")),
                              checkboxGroupInput(
                                inputId = ns("selectedColumnsResultsRP"),
                                label = "Select Columns to Display",
                                choices = c("GeneID", "Rank", "FileCount", "FileNames", "GenePositions", 
                                            "RP_stat", "PFP", "pvalue", "p.adjust"),
                                selected = c("GeneID", "Rank", "FileCount", "FileNames", "GenePositions", 
                                             "RP_stat", "PFP", "pvalue", "p.adjust"),
                                inline = FALSE)),
                            
                            conditionalPanel(
                              condition = sprintf("input['%s'] == 'RRApack'", ns("selectedPackage")),
                              checkboxGroupInput(
                                inputId = ns("selectedColumnsResultsRRA"),
                                label = "Select Columns to Display",
                                choices = c("GeneID", "Rank", "Score", "p.adjust", "FileCount", "FileNames", "GenePositions"),
                                selected = c("GeneID", "Rank", "Score", "p.adjust", "FileCount", "FileNames", "GenePositions"),
                                inline = FALSE)),
                            
                            tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
                            
                            # -------------------------------------------------------------------
                            # Upset Plot Settings
                            # -------------------------------------------------------------------
                            
                            tags$h5("Upset Plot Settings"),
                            sliderInput(
                              inputId = ns("upsetTextSize"),
                              label   = "Upset Text Size",
                              min     = 15,
                              max     = 40,
                              value   = 20),
                            
                            colorPickr(
                              inputId  = ns("upsetSetsColor"),
                              label    = "Sets Bar Color (horizontal)",
                              selected = "#2ba915",
                              opacity  = TRUE),
                            
                            colorPickr(
                              inputId  = ns("upsetIntersectionColor"),
                              label    = "Intersection Bar Color (vertical)",
                              selected = "#0838a0",
                              opacity  = TRUE),
                            
                            tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
                            
                            # -------------------------------------------------------------------
                            # Heatmap Settings
                            # -------------------------------------------------------------------
                            tags$h5("Heatmap Settings"),
                            sliderInput(
                              inputId = ns("heatmapTickSize"),
                              label   = "Heatmap Tick Size",
                              min     = 5,
                              max     = 30,
                              value   = 12),
                            
                            sliderInput(
                              inputId = ns("heatmapTitleSize"),
                              label   = "Heatmap Title Size",
                              min     = 10,
                              max     = 30,
                              value   = 17),
                            
                            selectInput(
                              inputId = ns("heatmapColorscale"),
                              label   = "Heatmap Colorscale",
                              choices = c("Viridis", "Cividis", "Portland"),
                              selected = "Viridis"),
                            
                            tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
                            
                            # -------------------------------------------------------------------
                            # Enrichment Plot Settings (Dot/Bar Plot)
                            # -------------------------------------------------------------------
                            tags$h5("Enrichment table Settings"),
                            checkboxGroupInput(
                              inputId = ns("selectedColumnsRank"),
                              label   = "Select Columns to Display in the Table",
                              choices = c("ID", "Description", "GeneRatio", "BgRatio", 
                                          "pvalue", "p.adjust", "Count", "geneID"),
                              selected = c("ID", "Description", "GeneRatio", "BgRatio", 
                                           "pvalue", "p.adjust", "Count", "geneID"),
                              inline   = FALSE),
                            
                            tags$hr(style = "margin-top: 10px; margin-bottom: 10px;"),
                            
                            tags$h5("Enrichment Plot Settings"),
                            sliderInput(
                              inputId = ns("nPlotTerms"),
                              label   = "Number of Terms to Show",
                              min     = 1,
                              max     = 100,
                              value   = 10),
                            
                            selectInput(
                              inputId = ns("plotLabel"),
                              label   = "Y-Axis Label",
                              choices = c(
                                "Term"         = "ID",
                                "Description"  = "Description",
                                "Term + Description" = "ID_Description"),
                              selected = "ID_Description"),
                            
                            selectInput(
                              inputId = ns("plotType"),
                              label   = "Plot Type",
                              choices = c("Dot Plot" = "dotplot", "Bar Plot" = "barplot"),
                              selected = "dotplot"),
                            
                            colorPickr(
                              inputId  = ns("colorLow"),
                              label    = "Low Color",
                              selected = "#0838a0",
                              opacity  = TRUE),
                            
                            colorPickr(
                              inputId  = ns("colorHigh"),
                              label    = "High Color",
                              selected = "#2ba915",
                              opacity  = TRUE),
                            
                            sliderInput(
                              inputId = ns("textSize"),
                              label   = "Text Size",
                              min     = 5,
                              max     = 20,
                              value   = 12))),
                        
                        #' @RUN_BUTTON
                        actionButton(inputId = ns("Run"), label = "Run Analysis", style = "width: 100%; margin-top: 15px;")),
      
      # -------------------------------------------------------------------------
      # Display results in the form of graphs and tables
      # -------------------------------------------------------------------------
      
      card(
        tabsetPanel(id = ns("mainTabs"),
                    
                    # Pestaña 1
                    tabPanel("MetaRanking Results", br(),
                             # Main Table
                             shinycssloaders::withSpinner(
                               DT::dataTableOutput(ns("Analysis_results_table"))), br(),
                             # Download buttons
                             uiOutput(ns("downloadButtonUI")), br(),
                             # Excluded table
                             shinyjs::hidden(
                               div(id = ns("excludedGenesSection"), br(),
                                   h5("Excluded terms (below frequency threshold):"),
                                   DT::dataTableOutput(ns("excludedGenesTable")))), br(),
                             # UpsetPlot
                             shinycssloaders::withSpinner(
                               plotOutput(ns("upsetPlot"), height = "600px")), br(),
                             # UpsetPlot download
                             uiOutput(ns("downloadButtonUpset")), br(),
                             # Heatmap
                             shinycssloaders::withSpinner(
                               plotlyOutput(ns("heatmapPlot"), height = "800px")), br(),
                             # Heatmap download
                             uiOutput(ns("downloadHeatmapUI"))),
                    
                    # Pestaña 2
                    tabPanel("Enrichment Results", br(),
                             # Main Enrichment Table
                             shinycssloaders::withSpinner(
                               DT::dataTableOutput(ns("enrichmentTable"))), br(),
                             # Download Buttons
                             uiOutput(ns("downloadButtonUI2")), br(),
                             # Dot/Bar plot
                             shinycssloaders::withSpinner(
                               uiOutput(ns("enrichmentDotPlotUI"))),
                             # Dot/Bar plot download
                             uiOutput(ns("downloadPlotUI3")))))))}



mod_metaRank_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    example_files <- c(
      "./modules/metarank/example_data/data1_1.tsv",
      "./modules/metarank/example_data/data2_2.tsv",
      "./modules/metarank/example_data/data3_3.tsv",
      "./modules/metarank/example_data/data4_4.tsv")
    
    example_files_rra <- c(
      "./modules/metarank/example_data/gene1.txt",
      "./modules/metarank/example_data/gene2.txt",
      "./modules/metarank/example_data/gene3.txt",
      "./modules/metarank/example_data/gene4.txt")
    
    annotation_dir <- file.path(getwd(), "database_annotations")
    
    #----------------------------------------------------------------------------
    # 1) CAPTURA DE INPUT (storedInput)
    #----------------------------------------------------------------------------
    storedInput <- reactiveVal(NULL)
    
    observeEvent(input$Run, {
      req(input$inputMethod, input$selectedPackage)
      
      pkg <- input$selectedPackage
      
      # ------------------------------------------------------------------------
      # 1. PREPARACIÓN DE DATOS PRELIMINAR PARA VALIDACIÓN
      # ------------------------------------------------------------------------
      check_data <- list()
      is_example <- FALSE
      
      if (input$inputMethod == "files") {
        is_example <- if (pkg == "Rankprodpack") isTRUE(input$useExampleData) else isTRUE(input$useExampleData_rra)
        
        if (is_example) {
          check_data <- list(type = "files", paths = "example_data_active")
        } else {
          files_obj <- if (pkg == "Rankprodpack") input$files else input$files_rra
          check_data <- list(type = "files", paths = files_obj$datapath)
        }
        
      } else {
        txt_to_check <- ""
        if (pkg == "Rankprodpack") {
          txt_to_check <- if (input$pasteFormat == "csv") input$pasteText else input$pasteText1
        } else {
          txt_to_check <- input$pasteText_rra}
        check_data <- list(type = "paste", text = txt_to_check)}
      
      # ------------------------------------------------------------------------
      # 2. VALIDACIÓN 1: EXISTENCIA
      # ------------------------------------------------------------------------
      if (!validate_input_exists(check_data)) {
        storedInput(NULL)
        return() }
      
      # ------------------------------------------------------------------------
      # 3A. VALIDACIÓN RRA
      # ------------------------------------------------------------------------
      if (pkg == "RRApack" && !is_example) {
        if (input$inputMethod == "files") {
          if (!validate_file_format_rra(input, session)) {
            storedInput(NULL); return()
          }
        } else {
          if (!validate_paste_format_rra(input, session)) {
            storedInput(NULL); return()}}}
      
      # ------------------------------------------------------------------------
      # 3B. VALIDACIÓN RANKPROD
      # ------------------------------------------------------------------------
      if (pkg == "Rankprodpack" && !is_example) {
        if (input$inputMethod == "files") {
          if (!validate_file_format_rp(input, session)) {
            storedInput(NULL); return()
          }
        } else {
          if (!validate_paste_format_rp(input, session)) {
            storedInput(NULL); return()}}}
      
      # ------------------------------------------------------------------------
      # 4. ALMACENAMIENTO DE DATOS
      # ------------------------------------------------------------------------
      if (input$inputMethod == "files") {
        if (pkg == "Rankprodpack") {
          if (isTRUE(input$useExampleData)) {
            storedInput(list(type = "files", package = pkg, paths = example_files, names = basename(example_files), is_example = TRUE))
          } else {
            storedInput(list(type = "files", package = pkg, paths = input$files$datapath, names = input$files$name, is_example = FALSE))
          }
        } else {
          if (isTRUE(input$useExampleData_rra)) {
            storedInput(list(type = "files", package = pkg, paths = example_files_rra, names = basename(example_files_rra), is_example = TRUE))
          } else {
            storedInput(list(type = "files", package = pkg, paths = input$files_rra$datapath, names = input$files_rra$name, is_example = FALSE))
          }}
        
      } else if (input$inputMethod == "paste") {
        if (pkg == "Rankprodpack") {
          req(input$pasteFormat)
          txt <- if(input$pasteFormat == "csv") input$pasteText else input$pasteText1
          storedInput(list(type = "paste", package = pkg, text = txt, format = input$pasteFormat))
        } else {
          storedInput(list(type = "paste", package = pkg, text = input$pasteText_rra))
        }}})
    
    
    #----------------------------------------------------------------------------
    # 2) PROCESAMIENTO DE DATOS (gene_data) 
    #----------------------------------------------------------------------------
    gene_data <- reactive({
      input$Run 
      
      inputData <- storedInput()
      req(inputData)
      
      withProgress(message = "Processing input data...", value = 0, {
        lists <- list()
        
        if (inputData$package == "Rankprodpack") {
          
          if (inputData$type == "files") {
            n_files <- length(inputData$paths)
            for (i in seq_along(inputData$paths)) {
              incProgress(1/n_files, detail = paste("Reading file:", inputData$names[i]))
              
              content <- tryCatch({
                ext <- tolower(tools::file_ext(inputData$names[i]))
                sep_char <- if(ext == "csv") "," else "\t"
                
                df <- read.table(inputData$paths[i], header = FALSE, sep = sep_char, 
                                 stringsAsFactors = FALSE, fill = TRUE, 
                                 quote = "", comment.char = "")
                
                if(ncol(df) >= 2) {
                  df <- df[, 1:2]
                  colnames(df) <- c("Gene", "Stat.data")
                  
                  df <- df[!is.na(df$Gene) & trimws(df$Gene) != "", ]
                  df$Gene <- sapply(strsplit(as.character(df$Gene), "//"), `[`, 1)
                  
                  if(is.character(df$Stat.data)) df$Stat.data <- gsub(",", ".", df$Stat.data)
                  df$Stat.data <- suppressWarnings(as.numeric(df$Stat.data))
                  df <- df[!is.na(df$Stat.data), ]
                  df <- df[!duplicated(df$Gene), ]
                  df
                } else { NULL }
              }, error = function(e) NULL)
              
              if (!is.null(content) && nrow(content) > 0) {
                clean_name <- inputData$names[i]
                lists[[clean_name]] <- content}}
            
          } else {
            incProgress(0.5, detail = "Parsing pasted text...")
            raw_text <- inputData$text
            fmt <- inputData$format
            blocks <- strsplit(raw_text, "###", fixed = TRUE)[[1]]
            n_blocks <- length(blocks)
            sep_char <- if(fmt == "csv") "," else "\t"
            
            for(i in seq_along(blocks)) {
              incProgress(0.5/n_blocks, detail = paste("Processing list", i))
              b <- trimws(blocks[i])
              if(nzchar(b)) {
                df <- tryCatch({
                  read.table(text = b, header = FALSE, sep = sep_char, 
                             stringsAsFactors = FALSE, fill = TRUE, quote = "")
                }, error = function(e) NULL)
                
                if(!is.null(df) && ncol(df) >= 2) {
                  colnames(df)[1:2] <- c("Gene", "Stat.data")
                  
                  df <- df[!is.na(df$Gene) & trimws(df$Gene) != "", ]
                  df$Gene <- sapply(strsplit(as.character(df$Gene), "//"), `[`, 1)
                  if(is.character(df$Stat.data)) df$Stat.data <- gsub(",", ".", df$Stat.data)
                  df$Stat.data <- suppressWarnings(as.numeric(df$Stat.data))
                  df <- df[!is.na(df$Stat.data), ]
                  df <- df[!duplicated(df$Gene), ]
                  
                  if(nrow(df) > 0) {
                    lists[[paste0("List_", i)]] <- df
                  }}}}}
          
          # --- LOGIC PARA RRA (Vectores) ---
        } else {
          
          if (inputData$type == "files") {
            # Archivos RRA
            n_files <- length(inputData$paths)
            for (i in seq_along(inputData$paths)) {
              incProgress(1/n_files, detail = paste("Reading file:", inputData$names[i]))
              lines <- readLines(inputData$paths[i], warn = FALSE, encoding = "UTF-8")
              genes <- trimws(lines)
              genes <- genes[nzchar(genes)]
              genes <- unique(genes)
              if (length(genes) > 0) {
                lists[[inputData$names[i]]] <- genes
              }}
            
          } else {
            # Paste RRA
            incProgress(0.5, detail = "Parsing pasted text...")
            raw_text <- inputData$text
            blocks <- strsplit(raw_text, "###", fixed = TRUE)[[1]]
            n_blocks <- length(blocks)
            
            for(i in seq_along(blocks)) {
              incProgress(0.5/n_blocks, detail = paste("Processing list", i))
              lines <- unlist(strsplit(blocks[i], "\n"))
              genes <- trimws(lines)
              genes <- genes[nzchar(genes)]
              genes <- unique(genes)
              if (length(genes) > 0) {
                lists[[paste0("List_", i)]] <- genes
              }}}}
        
        incProgress(1, detail = "Finalizing data import...")
        
        if (!validate_min_lists_metarank(lists)) {
          return(NULL)
        }
        
        return(lists)
      })})
    
    
    #----------------------------------------------------------------------------
    # DATA PREPARING 
    #----------------------------------------------------------------------------
    prepared_data <- reactive({
      
      withProgress(message = "Preparing data...", value = 0, {
        incProgress(0.10, detail = "Filtering valid input datasets...")
        raw <- gene_data()
        req(raw)
        is_rra <- input$selectedPackage == "RRApack"
        
        is_valid <- if (is_rra) function(x) length(x) > 0 else function(x) is.data.frame(x) && nrow(x) > 0
        
        data_list <- Filter(is_valid, raw)
        names_list <- names(data_list)
        req(length(data_list) > 0)
        
        incProgress(0.20, detail = "Extracting genes from input lists...")
        all_genes_list <- if (is_rra) data_list else
          lapply(data_list, function(df) if ("Gene" %in% colnames(df)) df$Gene else character(0))
        
        all_genes_list <- Filter(function(x) length(x) > 0, all_genes_list)
        req(length(all_genes_list) > 0)
        
        incProgress(0.25, detail = "Validating gene overlap threshold...")
        
        # ----------------------------------------------------------------------
        # VALIDACIÓN DEL THRESHOLD
        # ----------------------------------------------------------------------
        validated_threshold <- validate_gene_threshold(all_genes_list, names_list, input, session)
        
        if (is.null(validated_threshold)) {
          return(NULL) 
        }
        
        min_appearance <- validated_threshold
        app_df <- count_gene_appearance_rra(all_genes_list, names_list)
        genes_keep <- app_df$Gene[app_df$Count >= min_appearance]
        
        incProgress(0.20, detail = "Filtering genes by appearance threshold...")
        
        valid <- if (is_rra) {
          lapply(data_list, function(lst) intersect(lst, genes_keep))
        } else {
          lapply(data_list, function(df)
            if ("Gene" %in% colnames(df)) filter(df, Gene %in% genes_keep) else df)}
        
        excluded <- if (is_rra) {
          lapply(data_list, function(lst) data.frame(Gene = setdiff(lst, genes_keep)))
        } else {
          lapply(data_list, function(df)
            if ("Gene" %in% colnames(df)) filter(df, !Gene %in% genes_keep) else df)}
        
        incProgress(0.10, detail = "Data preparation completed.")
        list(valid = valid, excluded = excluded, files = names_list)})})
    
    
    #----------------------------------------------------------------------------
    # PERFORM THE META-ANALYSIS 
    #----------------------------------------------------------------------------
    combined_results <- eventReactive(input$Run, {
      
      raw <- gene_data()
      data <- prepared_data()
      req(data$valid)
      
      if (!validate_origin_rp(input, session, data)) {
        return(NULL)
      }
      
      withProgress(message = "Running meta-analysis...", value = 0, {
        incProgress(0.10, detail = "Preparing input matrices...")
        
        if (input$selectedPackage == "Rankprodpack") {
          direction <- if (input$ranking_direction == "asc") "ascending" else "descending"
          current_direction <- input$ranking_direction
          na_rm <- input$NA_management != "ignore"
          statdata <- generate_statdata_matrix(
            gene_data = data$valid,
            na_management = input$NA_management,
            ranking_direction = direction)
          
          cl <- rep(1, ncol(statdata))
          incProgress(0.20, detail = "Computing RankProd statistics...")
          if (input$RankBasedMethod == "rp") {
            res <- perform_RP(
              statdata, cl,
              logged = TRUE,
              na.rm = na_rm,
              ranking_direction = direction,
              progress = incProgress)
            
          } else {
            # Asumimos origen válido
            origin <- as.numeric(strsplit(gsub("\\s+", "", input$origin), ",")[[1]])
            res <- perform_RP_advance(
              statdata, cl, origin,
              logged = TRUE,
              na.rm = na_rm,
              ranking_direction = direction,
              progress = incProgress)}
          
          incProgress(0.20, detail = "Computing aggregate gene statistics...")
          app <- count_gene_appearance_with_files(data$valid, data$files)
          apply_penalty <- input$extra_penalty_switch &&
            input$NA_management %in% c("impute", "ignore")
          tbl <- if (apply_penalty) {
            create_combined_rank_table2(res, app)
          } else {
            create_combined_rank_table(res, app)}
          
        } else {
          incProgress(0.35, detail = "Running RRA analysis pipeline...")
          full_flag <- input$rra_na_handling == "ignore"
          tbl <- run_rra_analysis(data$valid, method = input$rra_method, full = full_flag)}
        
        incProgress(0.15, detail = "Merging gene coordinates...")
        pos_df <- get_gene_positions(raw, ranking_direction = current_direction)      
        tbl <- dplyr::left_join(tbl, pos_df, by = "GeneID")
        
        incProgress(0.10, detail = "Meta-analysis completed.")
        tbl})})
    
    
    #----------------------------------------------------------------------------
    # DATA SELECT FOR THE ENRICHMENT ANALYSIS
    #----------------------------------------------------------------------------
    preview_genes <- reactive({
      input$inputMethod
      input$selectedPackage
      input$useExampleData; input$useExampleData_rra
      input$files; input$files_rra
      input$pasteText; input$pasteText1; input$pasteText_rra
      input$pasteFormat
      
      pkg <- input$selectedPackage
      lists <- list()
      
      tryCatch({
        # --- CASO 1: RANKPROD (RP) ---
        if (pkg == "Rankprodpack") {
          
          if (input$inputMethod == "files") {
            # Archivos RP
            paths <- if (isTRUE(input$useExampleData)) example_files else input$files$datapath
            names <- if (isTRUE(input$useExampleData)) basename(example_files) else input$files$name
            req(paths)
            
            for (i in seq_along(paths)) {
              ext <- tolower(tools::file_ext(names[i]))
              sep_char <- if(ext == "csv") "," else "\t"
              df <- tryCatch(read.table(paths[i], header = FALSE, sep = sep_char, fill=TRUE, quote=""), error=function(e) NULL)
              if (!is.null(df) && ncol(df) >= 1) {
                g <- as.character(df[[1]])
                g <- g[!is.na(g) & trimws(g) != ""]
                g <- sapply(strsplit(g, "//"), `[`, 1)
                lists[[paste0("File_RP_", i)]] <- unique(g)
              }
            }
            
          } else {
            # Paste RP
            req(input$pasteFormat)
            txt <- if(input$pasteFormat == "csv") input$pasteText else input$pasteText1
            sep_char <- if(input$pasteFormat == "csv") "," else "\t"
            
            if (!is.null(txt) && nzchar(txt)) {
              blocks <- strsplit(txt, "###", fixed = TRUE)[[1]]
              for (i in seq_along(blocks)) {
                b <- trimws(blocks[i])
                if (nzchar(b) && !grepl("^#", b)) {
                  df <- tryCatch(read.table(text = b, header = FALSE, sep = sep_char, fill=TRUE, quote=""), error=function(e) NULL)
                  if (!is.null(df) && ncol(df) >= 1) {
                    g <- as.character(df[[1]])
                    g <- g[!is.na(g) & trimws(g) != ""]
                    g <- sapply(strsplit(g, "//"), `[`, 1)
                    lists[[paste0("Paste_RP_", i)]] <- unique(g)
                  }
                }
              }
            }
          }
          
          # --- CASO 2: RRA ---
        } else {
          
          if (input$inputMethod == "files") {
            paths <- if (isTRUE(input$useExampleData_rra)) example_files_rra else input$files_rra$datapath
            req(paths)
            for (i in seq_along(paths)) {
              l <- readLines(paths[i], warn = FALSE)
              l <- trimws(l); l <- l[nzchar(l)]
              lists[[paste0("File_RRA_", i)]] <- unique(l)
            }
            
          } else {
            txt <- input$pasteText_rra
            if (!is.null(txt) && nzchar(txt)) {
              blocks <- strsplit(txt, "###", fixed = TRUE)[[1]]
              for (i in seq_along(blocks)) {
                b <- trimws(blocks[i])
                if (nzchar(b)) {
                  con <- textConnection(b)
                  lines <- readLines(con, warn = FALSE)
                  close(con)
                  
                  l <- trimws(lines)
                  l <- l[nzchar(l)]   
                  l <- unique(l)      
                  
                  if (length(l) > 0) {
                    lists[[paste0("Paste_RRA_", i)]] <- l
                  }}}}}}
        return(lists)
      }, error = function(e) return(NULL))
    })
    
    
    observe({
      # --- DEPENDENCIAS REACTIVAS ---
      data_lists <- preview_genes()
      pkg <- input$selectedPackage
      
      slider_id <- if (pkg == "Rankprodpack") "gene_app_filter_rp" else "gene_app_filter_rra"
      min_app <- input[[slider_id]]
      
      # --- LÓGICA ---
      if (is.null(data_lists) || length(data_lists) == 0) {
        updateSliderInput(session, "topGenes", max = 100, value = 30)
        updateNumericInput(session, "topGenes_numeric", max = 100, value = 30)
        return()
      }
      
      if (is.null(min_app)) min_app <- 1
      
      all_genes_flat <- unlist(data_lists)
      gene_counts <- table(all_genes_flat)
      valid_genes_count <- sum(gene_counts >= min_app)
      max_slider <- max(1, valid_genes_count)
      current_val <- isolate(input$topGenes)
      safe_current <- if (is.null(current_val)) 30 else current_val
      new_val <- min(safe_current, max_slider)
      updateSliderInput(session, "topGenes", 
                        min = 1, 
                        max = max_slider, 
                        value = new_val)
      
      updateNumericInput(session, "topGenes_numeric", 
                         min = 1, 
                         max = max_slider, 
                         value = new_val)
    })
    observeEvent(input$topGenes, {
      val <- input$topGenes
      if (!is.null(val) && !is.na(val)) {
        updateNumericInput(session, "topGenes_numeric", value = val)}
    }, ignoreInit = TRUE)
    
    observeEvent(input$topGenes_numeric, {
      req(input$topGenes)
      val <- input$topGenes_numeric
      if (!is.null(val)) {
        updateSliderInput(session, "topGenes", value = val)
      }
    }, ignoreInit = TRUE)
    
    
    selected_genes <- reactive({
      tbl <- combined_results()
      req(tbl)
      gene_col <- if("GeneID" %in% colnames(tbl)) "GeneID" else if("Gene" %in% colnames(tbl)) "Gene" else colnames(tbl)[1]
      nsel <- min(input$topGenes, nrow(tbl))
      head(tbl[[gene_col]], nsel)
    })
    
    
    #----------------------------------------------------------------------------
    # ENRICHMENT ANALYSIS OF THE RANKED GENES
    #----------------------------------------------------------------------------
    enrichment_results <- reactiveVal(NULL)
    analysis_params <- reactiveVal(list())
    selected_cols_proxy <- reactiveVal(NULL)
    
    observeEvent(input$Run, {
      req(input$oraDatabase, selected_genes)
      
      definitive_cols <- if (input$oraDatabase == "Personalized") {
        c("Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust",
          "FeatureCount", "FeatureID")
      } else {
        c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust",
          "GeneCount", "GeneID")}
      
      analysis_params(list(
        database = input$oraDatabase,
        columns = definitive_cols))
      
      if (input$oraDatabase == "Personalized") {
        # 1A. Comprobar si existe el archivo (NUEVO)
        if (!check_database_exists_metarank(input$oraDatabase, input$personalizedAnnotations)) {
          enrichment_results(NULL)
          return()
        }
        
        # 1B. Comprobar el formato del archivo
        if (!check_annotation_file_format(input$personalizedAnnotations$datapath, session)) {
          enrichment_results(NULL)
          return()
        }
      }
      
      if (input$enableEnrichment) {
        withProgress(message = "Performing enrichment...", value = 0, {
          
          incProgress(0.10, detail = "Validating gene identifiers...")
          genes_to_test <- selected_genes()
          
          # 2. Validation: Gene IDs & Organism (only if NOT Personalized)
          if (input$oraDatabase != "Personalized") {
            full_results <- combined_results()
            gene_file_map <- NULL
            if (!is.null(full_results) && "GeneID" %in% colnames(full_results) && "FileNames" %in% colnames(full_results)) {
              gene_file_map <- full_results[, c("GeneID", "FileNames")]
            }
            
            # --- VALIDACIÓN 1: IDs ---
            genes_to_test <- validate_gene_ids(genes_to_test, input$IDtype, session, gene_file_map)
            if (is.null(genes_to_test)) { enrichment_results(NULL); return() } # Si retorna NULL (>50% error), paramos
            
            # --- VALIDACIÓN 2: Organismo ---
            genes_to_test <- validate_organism(genes_to_test, input$species, input$IDtype, session, gene_file_map)
            if (is.null(genes_to_test)) { enrichment_results(NULL); return() } # Si retorna NULL (>50% error), paramos
            
            if (length(genes_to_test) == 0) {
              showModal(modalDialog(
                title=tags$div(icon("ban"), " No Genes Left", style="color: #c0392b; font-weight:bold; font-size:18px;"),
                tags$div(class="alert alert-danger", "No valid genes remained after filtering by ID Type and Organism."), 
                easyClose=TRUE))
              enrichment_results(NULL); return()
            }
          }
          
          incProgress(0.40, detail = "Running enrichment test...")
          
          # Execution
          result <- switch(
            input$oraDatabase,
            "go" = ORA_GO(gene_list = genes_to_test,
                          organism = input$species,
                          ontology = toupper(input$ontology),
                          pvalue = 1, 
                          ID_type = input$IDtype),
            "kegg" = ORA_KEGG(gene_list = genes_to_test,
                              organism = input$species,
                              pvalue = 1,
                              ID_type = input$IDtype),
            "reactome" = ORA_REACTOME(gene_list = genes_to_test,
                                      organism = input$species,
                                      pvalue = 1,
                                      ID_type = input$IDtype),
            "Personalized" = ORA_PERSONALIZED(gene_list = genes_to_test,
                                              annotation_file = input$personalizedAnnotations$datapath))
          
          incProgress(0.30, detail = "Finalizing enrichment output...")
          
          # 3. Validation: Check Empty Results
          if (!check_empty_enrichment_results(result, session)) {
            enrichment_results(NULL) 
            return()
          }
          
          enrichment_results(result)
          selected_cols_proxy(definitive_cols)
          
          incProgress(0.20, detail = "Enrichment completed.")})
        
      } else {
        enrichment_results(NULL)
        selected_cols_proxy(NULL)}})
    
    
    
    #----------------------------------------------------------------------------
    # DATA REPRESENTATION OF META-ANALYSIS (CARD)
    #----------------------------------------------------------------------------
    
    # RANKPROD ANALYSIS (TABLE + 2 DOWNLOAD BUTTONS + 2 BUTTONS OF EXCLUDED GENES)
    #----------------------------------------------------------------------------
    output$Analysis_results_table <- DT::renderDataTable({
      req(combined_results()) 
      req(input$selectedPackage) 
      
      selected_columns <- if (input$selectedPackage == "Rankprodpack") {
        input$selectedColumnsResultsRP
      } else if (input$selectedPackage == "RRApack") {
        input$selectedColumnsResultsRRA
      } else {
        return(data.frame(Message = "No package selected"))}
      
      data <- combined_results()
      if (nrow(data) == 0) return(data.frame(Message = "No data available"))
      
      selected_data <- data[, selected_columns, drop = FALSE]
      if (nrow(selected_data) == 0) return(data.frame(Message = "No data available"))
      
      num_rows <- input$Analysis_results_table_length
      table_height <- paste0(num_rows * 40, "px")
      
      colTooltips4 <- list(
        GeneID = "Unique gene identifier, which can be a HUGO symbol (e.g., TP53), an Entrez ID (e.g., 7157), or an Ensembl ID (e.g., ENSG00000141510), depending on the selected input format.",
        Rank = "Consensus ranking of the gene across all input lists; lower values indicate higher overall relevance or consistency among the lists.",
        FileCount = "Number of input gene lists in which this gene appears; a higher count suggests greater consistency across datasets.",
        FileNames = "Names of the input files where the gene was found, separated by spaces; useful for identifying the sources supporting the gene's relevance.",
        RP_stat = "Rank Product statistic calculated to assess the significance of the gene's ranking across multiple lists; lower values suggest higher significance.",
        PFP = "Estimated Proportion of False Positives, analogous to the False Discovery Rate (FDR); lower values indicate more reliable findings.",
        pvalue = "Raw p-value from the meta-analysis, indicating the probability of observing the gene's ranking by chance; lower values suggest higher significance.",
        p.adjust = "Adjusted p-value accounting for multiple hypothesis testing using the Benjamini-Hochberg method; helps control the FDR.",
        GenePositions = "Rank positions of the gene in each individual input list; provides insight into the gene's performance across different datasets.",
        Score = "Aggregate score computed to reflect the gene's overall ranking consistency; lower scores indicate higher agreement across lists.")
      
      dt <- datatable(
        selected_data,
        rownames  = FALSE,
        selection = "multiple",
        filter    = "top",
        escape    = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          scrollY = table_height,
          autoWidth = FALSE,
          columnDefs = list(list(
            targets = "_all",
            className = "dt-wrap",
            width = paste0(round(100 / max(1, ncol(selected_data))), "%")
          )),
          headerCallback = htmlwidgets::JS(
            "function(thead, data, start, end, display){",
            "  var tooltips = ", jsonlite::toJSON(colTooltips4), ";",
            "  $('th', thead).each(function(){",
            "    var txt = $(this).text().trim();",
            "    if(tooltips[txt]){",
            "      var label = $('<span>').addClass('dt-header-label').text(txt);",
            "      var icon  = $('<span>?</span>')",
            "                     .addClass('dt-info-icon')",
            "                     .attr('title', tooltips[txt]);",
            "      $(this).empty().append($('<div>').addClass('dt-header-with-icon')",
            "                               .append(label).append(icon));",
            "    }",
            "  });",
            "}")))
      
      if (input$selectedPackage == "Rankprodpack") {
        dt <- dt %>% formatSignif(columns = intersect(c("RP_stat", "PFP", "pvalue", "p.adjust"), colnames(selected_data)), digits = 3)
      } else if (input$selectedPackage == "RRApack") {
        dt <- dt %>% formatSignif(columns = intersect(c("Score", "p.adjust"), colnames(selected_data)), digits = 3)}
      
      return(dt)})
    
    
    output$downloadButtonUI <- renderUI({
      req(combined_results())
      fluidRow(
        column(5, downloadButton(ns("downloadResultsCSV"), "Download CSV", style = "width: 100%;")),
        column(5, downloadButton(ns("downloadResultsTSV"), "Download TSV", style = "width: 100%;")),
        column(1, actionButton(ns("toggleExcludedGenes"), 
                               label = NULL, 
                               icon = icon("eye"), 
                               style = "width: 100%;", 
                               title = "Show/hide table of excluded genes with unique appearance")),
        column(1, downloadButton(ns("downloadExcludedGenesList"), 
                                 label = NULL, 
                                 icon = icon("download"), 
                                 style = "width: 100%;", 
                                 title = "Download excluded gene list (.txt)")))})
    
    # DOWNLOAD CSV
    output$downloadResultsCSV <- downloadHandler(
      filename = function() {
        paste0("MetaRank_Results_", Sys.Date(), ".csv")},
      content = function(file) {
        req(combined_results(), input$selectedPackage, input$Analysis_results_table_rows_all)
        selected_columns <- switch(input$selectedPackage,
                                   "Rankprodpack" = input$selectedColumnsResultsRP,
                                   "RRApack" = input$selectedColumnsResultsRRA, NULL)
        req(selected_columns)
        data <- combined_results()
        filtered_data <- data[input$Analysis_results_table_rows_all, selected_columns, drop = FALSE]
        write.csv(filtered_data, file, row.names = FALSE)})
    
    
    # DOWNLOAD TSV
    output$downloadResultsTSV <- downloadHandler(
      filename = function() {
        paste0("MetaRank_Results_", Sys.Date(), ".tsv")},
      content = function(file) {
        req(combined_results(), input$selectedPackage, input$Analysis_results_table_rows_all)
        selected_columns <- switch(input$selectedPackage,
                                   "Rankprodpack" = input$selectedColumnsResultsRP,
                                   "RRApack" = input$selectedColumnsResultsRRA, NULL)
        req(selected_columns)
        data <- combined_results()
        filtered_data <- data[input$Analysis_results_table_rows_all, selected_columns, drop = FALSE]
        write.table(filtered_data, file, sep = "\t", row.names = FALSE, quote = FALSE)})
    
    
    # DOWNLOAD EXCLUDED
    output$downloadExcludedGenesList <- downloadHandler(
      filename = function() paste0("Excluded_Genes_List_", Sys.Date(), ".txt"),
      content = function(file) {
        excluded <- prepared_data()$excluded
        file_names <- prepared_data()$files
        combined <- do.call(rbind, lapply(seq_along(excluded), function(i) {
          df <- excluded[[i]]
          if (nrow(df) > 0) df$Source <- file_names[i]
          df}))
        if (!"Gene" %in% colnames(combined)) return()
        
        combined <- combined |>
          group_by(Gene) |>
          summarise(Source = paste(unique(Source), collapse = ", "), .groups = "drop")
        
        writeLines(sort(unique(combined$Gene)), con = file)})
    
    
    # EXCLUDED TABLE
    output$excludedGenesTable <- DT::renderDT({
      prep <- prepared_data()
      req(prep$excluded, prep$files)
      
      combined <- do.call(rbind, lapply(seq_along(prep$excluded), function(i) {
        df <- prep$excluded[[i]]
        if (nrow(df) > 0) df$Source <- prep$files[i]
        df}))
      req(nrow(combined) > 0, "No excluded genes.")
      
      combined <- combined |>
        group_by(Gene) |>
        summarise(Source = paste(unique(Source), collapse = ", "), .groups = "drop")
      
      raw   <- gene_data()
      names <- prep$files
      all_genes <- if (input$selectedPackage=="RRApack") raw else lapply(raw, `[[`, "Gene")
      app_df <- count_gene_appearance_rra(all_genes, names)
      combined <- combined |>
        left_join(app_df[, c("Gene","Count")], by = "Gene")
      
      datatable(
        setNames(combined[, c("Gene", "Count", "Source")], c("GeneID", "FileCount", "FileNames")),
        options   = list(pageLength = 5, autoWidth = TRUE),
        rownames  = FALSE,
        selection = "multiple")})
    
    
    
    
    # UPSET PLOT (PLOT + 2 DOWNLOAD BUTTONS)
    #----------------------------------------------------------------------------
    upsetPlot <- reactive({
      req(gene_data())
      req(combined_results())
      create_upset_plot(
        gene_data = gene_data(),
        text_size = input$upsetTextSize,
        sets_color = input$upsetSetsColor,
        intersection_color = input$upsetIntersectionColor)})
    
    output$upsetPlot <- renderPlot({
      req(upsetPlot())
      upsetPlot()})
    
    output$downloadButtonUpset <- renderUI({
      req(upsetPlot())
      fluidRow(
        column(6, downloadButton(ns("downloadUpsetPNG"), "Download PNG", style = "width: 100%;")),
        column(6, downloadButton(ns("downloadUpsetJPG"), "Download JPG", style = "width: 100%;")))})
    
    # DOWNLOAD PNG
    output$downloadUpsetPNG <- downloadHandler(
      filename = function() {
        paste0("UpSetR_Plot_", Sys.Date(), ".png")},
      content = function(file) {
        save_upset_plot_transparent(upsetPlot(), file)},
      contentType = "image/png")
    
    # DOWNLOAD JPG
    output$downloadUpsetJPG <- downloadHandler(
      filename = function() {
        paste0("UpSetR_Plot_", Sys.Date(), ".jpg")},
      content = function(file) {
        jpeg(file, width = 14.5 * 96, height = 8.33 * 96, res = 96)  # Abrir dispositivo JPEG
        print(create_upset_plot(
          gene_data = gene_data(),
          text_size = input$upsetTextSize,
          sets_color = input$upsetSetsColor,
          intersection_color = input$upsetIntersectionColor))
        dev.off()},
      contentType = "image/jpg")
    
    
    
    # HEATMAP (PLOT + MULTIPLE DOWNLOAD BUTTON (PNG, JPG, HTML))
    #----------------------------------------------------------------------------
    heatmap_obj <- reactive({
      req(gene_data())
      req(combined_results())
      generate_heatmap(
        gene_data = gene_data(),
        text_size = input$heatmapTickSize,
        title_size = input$heatmapTitleSize,
        tick_size = input$heatmapTickSize,
        colorscale = input$heatmapColorscale)})
    
    heatmap_flat_obj <- reactive({
      req(gene_data())
      req(combined_results())
      generate_flat_heatmap(
        gene_data = gene_data(),
        text_size = input$heatmapTickSize,
        title_size = input$heatmapTitleSize,
        tick_size = input$heatmapTickSize,
        colorscale = input$heatmapColorscale)})
    
    output$heatmapPlot <- renderPlotly({
      req(heatmap_obj())
      ggplotly(heatmap_obj())})
    
    output$downloadHeatmapUI <- renderUI({
      req(combined_results())
      dropdownButton(
        label = "Download Heatmap",
        icon = icon("download"),
        circle = FALSE,
        status = "secondary",
        width = "300px",
        tags$h4("Download Format"),
        downloadButton(ns("downloadHeatmapHTML"), "HTML", class = "btn-block"),
        downloadButton(ns("downloadHeatmapPNG"), "PNG", class = "btn-block"),
        downloadButton(ns("downloadHeatmapJPG"), "JPG", class = "btn-block"))})
    
    # DOWNLOAD HTML 
    output$downloadHeatmapHTML <- downloadHandler(
      filename = function() {
        paste("MetaRank_Heatmap_", Sys.Date(), ".html", sep = "")},
      content = function(file) {
        htmlwidgets::saveWidget(ggplotly(heatmap_obj()), file, selfcontained = TRUE)})
    
    # DOWNLOAD PNG
    output$downloadHeatmapPNG <- downloadHandler(
      filename = function() {
        paste("MetaRank_Heatmap_", Sys.Date(), ".png", sep = "")},
      content = function(file) {
        ggsave(file, plot = heatmap_flat_obj(), device = "png", width = 10, height = 8, dpi = 300)})
    
    # DOWNLOAD JPG
    output$downloadHeatmapJPG <- downloadHandler(
      filename = function() {
        paste("MetaRank_Heatmap_", Sys.Date(), ".jpg", sep = "")},
      content = function(file) {
        ggsave(file, plot = heatmap_flat_obj(), device = "jpeg", width = 10, height = 8, dpi = 300)})
    
    
    
    #----------------------------------------------------------------------------
    # DATA REPRESENTATION OF ENRICHMENT ANALYSIS (CARD)
    #----------------------------------------------------------------------------
    
    
    observeEvent(input$oraDatabase, {
      if (input$oraDatabase == "Personalized") {
        new_choices <- c("Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "FeatureCount", "FeatureID")
      } else {
        new_choices <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "GeneCount", "GeneID")}
      new_selected <- new_choices
      updateCheckboxGroupInput(
        session,
        "selectedColumnsRank",
        choices = new_choices, 
        selected = new_selected )
      
    }, ignoreInit = FALSE)
    
    
    observeEvent(input$selectedColumnsRank, {
      req(analysis_params(), !is.null(analysis_params()$columns))
      
      params <- analysis_params()
      if (all(input$selectedColumnsRank %in% params$columns)) {
        selected_cols_proxy(input$selectedColumnsRank)}
    }, ignoreNULL = FALSE, ignoreInit = TRUE)
    
    
    #----------------------------------------------------------------------------
    # ENRICHMENT ANALYSIS (TABLE + PLOT + DESCARGAS)
    #----------------------------------------------------------------------------
    
    output$enrichmentTable <- DT::renderDataTable({
      req(enrichment_results(), selected_cols_proxy())
      
      data <- as.data.frame(enrichment_results(), stringsAsFactors = FALSE)
      params <- analysis_params()
      is_personalized <- (!is.null(params$database) && params$database == "Personalized")
      
      if (!is_personalized) {
        if ("Count" %in% colnames(data)) colnames(data)[colnames(data) == "Count"] <- "GeneCount"
        if ("geneID" %in% colnames(data)) colnames(data)[colnames(data) == "geneID"] <- "GeneID"}
      
      cols_to_display <- intersect(selected_cols_proxy(), colnames(data))
      req(length(cols_to_display) > 0)
      selected_data <- data[, cols_to_display, drop = FALSE]
      
      if (is_personalized) {
        user_cols <- gsub("^GeneID$", "FeatureID", selected_data)
        user_cols <- gsub("^GeneCount$", "FeatureCount", selected_data)
      } else {
        user_cols <- gsub("^FeatureID$", "GeneID", selected_data)
        user_cols <- gsub("^FeatureCount$", "GeneCount", selected_data)}
      
      colTooltips <- if (is_personalized) {
        list(
          Description  = "Descriptive name from your annotation file.",
          GeneRatio    = "Proportion of your input features associated with this term.",
          BgRatio      = "Proportion of background features associated with this term.",
          pvalue       = "Raw p-value from the hypergeometric test.",
          p.adjust     = "Adjusted p-value (Benjamini–Hochberg).",
          FeatureCount = "Number of your input features mapped to this term.",
          FeatureID    = "List of your input features mapping to this term.")
      } else {
        list(
          ID           = "Stable identifier of the enriched term.",
          Description  = "Descriptive name of the term.",
          GeneRatio    = "Proportion of your input genes associated with this term.",
          BgRatio      = "Proportion of background genes associated with this term.",
          pvalue       = "Raw p-value from the hypergeometric test.",
          p.adjust     = "Adjusted p-value (Benjamini–Hochberg).",
          GeneCount    = "Number of your input genes mapped to this term.",
          GeneID       = "List of your input genes mapping to this term.")}
      
      id_col_name <- if (is_personalized) "FeatureID" else "GeneID"
      gene_col_idx <- which(colnames(selected_data) == id_col_name) - 1
      column_defs <- list(
        list(targets = "_all", className = "dt-wrap", width = paste0(round(100 / max(1, ncol(selected_data))), "%")))
      
      if (length(gene_col_idx) > 0) {
        column_defs <- c(column_defs,
                         list(list(
                           targets = gene_col_idx,
                           className = "dt-gene-wrap",
                           width = "300px",
                           render = htmlwidgets::JS(
                             "function(data, type, row, meta) {",
                             "  if (type === 'display' && data != null && data.length > 40) {",
                             "    return '<span title=\"' + data + '\">' + data.substr(0, 40) + '…</span>';",
                             "  } else {",
                             "    return data;",
                             "  }",
                             "}"))))}
      
      order_idx   <- which(colnames(selected_data) == "p.adjust")
      order_list  <- if (length(order_idx) > 0) list(list(order_idx - 1, "asc")) else list()
      num_rows    <- input$enrichmentTable_length
      table_h     <- paste0(num_rows * 40, "px")
      
      DT::datatable(
        selected_data,
        rownames = FALSE,
        selection = "multiple",
        filter = "top",
        escape = FALSE,
        options = list(
          pageLength = 10,
          scrollX = FALSE,
          scrollY = table_h,
          autoWidth = FALSE,
          columnDefs = column_defs,
          order = order_list,
          initComplete = htmlwidgets::JS(
            "function(settings, json) {",
            "  $(this.api().table().node()).css({'table-layout':'fixed','width':'100%'});",
            "}"),
          headerCallback = htmlwidgets::JS(
            "function(thead, data, start, end, display){",
            "  var tooltips = ", jsonlite::toJSON(colTooltips), ";",
            "  $('th', thead).each(function(){",
            "    var txt = $(this).text().trim();",
            "    if(tooltips[txt]){",
            "      var label = $('<span>').addClass('dt-header-label').text(txt);",
            "      var icon  = $('<span>?</span>').addClass('dt-info-icon').attr('title', tooltips[txt]);",
            "      $(this).empty().append($('<div>').addClass('dt-header-with-icon').append(label).append(icon));",
            "    }",
            "  });",
            "}"))) %>%
        formatSignif(columns = intersect(c("pvalue", "p.adjust"), colnames(selected_data)), digits = 4)})
    
    output$downloadButtonUI2 <- renderUI({
      req(combined_results())
      fluidRow(
        column(6, downloadButton(ns("downloadEnrichmentCSV"), "Download CSV", style = "width: 100%;")),
        column(6, downloadButton(ns("downloadEnrichmentTSV"), "Download TSV", style = "width: 100%;")))})
    
    # DOWNLOAD CSV
    output$downloadEnrichmentCSV <- downloadHandler(
      filename = function() { paste0("MetaRank_results_Enrich_", Sys.Date(), ".csv") },
      content = function(file) {
        req(enrichment_results())
        data <- as.data.frame(enrichment_results(), stringsAsFactors = FALSE)
        params <- analysis_params()
        is_personalized <- (!is.null(params$database) && params$database == "Personalized")
        
        if (!is_personalized) {
          if ("Count" %in% colnames(data)) colnames(data)[colnames(data) == "Count"] <- "GeneCount"
          if ("geneID" %in% colnames(data)) colnames(data)[colnames(data) == "geneID"] <- "GeneID"}
        
        user_cols <- input$selectedColumnsRank
        if (is_personalized) {
          user_cols <- gsub("^GeneID$", "FeatureID", user_cols)
          user_cols <- gsub("^GeneCount$", "FeatureCount", user_cols)
        } else {
          user_cols <- gsub("^FeatureID$", "GeneID", user_cols)
          user_cols <- gsub("^FeatureCount$", "GeneCount", user_cols)}
        
        cols_to_use <- intersect(user_cols, colnames(data))
        write.csv(data[, cols_to_use, drop = FALSE], file, row.names = FALSE)})
    
    output$downloadEnrichmentTSV <- downloadHandler(
      filename = function() { paste0("MetaRank_results_Enrich_", Sys.Date(), ".tsv") },
      content = function(file) {
        req(enrichment_results())
        data <- as.data.frame(enrichment_results(), stringsAsFactors = FALSE)
        params <- analysis_params()
        is_personalized <- (!is.null(params$database) && params$database == "Personalized")
        
        if (!is_personalized) {
          if ("Count" %in% colnames(data)) colnames(data)[colnames(data) == "Count"] <- "GeneCount"
          if ("geneID" %in% colnames(data)) colnames(data)[colnames(data) == "geneID"] <- "GeneID"}
        
        user_cols <- input$selectedColumnsRank
        if (is_personalized) {
          user_cols <- gsub("^GeneID$", "FeatureID", user_cols)
          user_cols <- gsub("^GeneCount$", "FeatureCount", user_cols)
        } else {
          user_cols <- gsub("^FeatureID$", "GeneID", user_cols)
          user_cols <- gsub("^FeatureCount$", "GeneCount", user_cols)}
        
        cols_to_use <- intersect(user_cols, colnames(data))
        write.table(data[, cols_to_use, drop = FALSE], file, row.names = FALSE, sep = "\t", quote = FALSE)})
    
    
    filtered_enrichment <- reactive({
      req(enrichment_results())
      if (!is.null(input$enrichmentTable_rows_all)) {
        enrichment_results()[as.numeric(input$enrichmentTable_rows_all), ]
      } else {
        enrichment_results()}})
    
    
    plot_obj <- reactive({
      req(filtered_enrichment())
      enriched_data <- filtered_enrichment()
      
      params <- analysis_params()
      is_personalized <- (!is.null(params$database) && params$database == "Personalized")
      
      if (!is_personalized) {
        if ("Count" %in% colnames(enriched_data)) colnames(enriched_data)[colnames(enriched_data) == "Count"] <- "GeneCount"
        if ("geneID" %in% colnames(enriched_data)) colnames(enriched_data)[colnames(enriched_data) == "geneID"] <- "GeneID"}
      
      num_terms <- min(input$nPlotTerms, nrow(enriched_data))
      data_plot <- enriched_data[order(enriched_data$p.adjust), ][1:num_terms, ]
      text_size <- input$textSize
      
      if (is_personalized) {
        data_plot$label_y <- data_plot$Description
        x_var <- "FeatureCount"
        
        p <- ggplot(data_plot, aes(x = .data[[x_var]], y = reorder(label_y, -p.adjust),
                                   text = paste0("Description: ", Description,
                                                 "<br>FeatureCount: ", .data[[x_var]],
                                                 "<br>p.adjust: ", signif(p.adjust, 3))))
      } else {
        data_plot$label_y <- switch(input$plotLabel,
                                    "ID" = data_plot$ID,
                                    "Description" = data_plot$Description,
                                    "ID_Description" = paste(data_plot$ID, data_plot$Description, sep = " - "))
        x_var <- "GeneCount"
        
        p <- ggplot(data_plot, aes(x = .data[[x_var]], y = reorder(label_y, -p.adjust),
                                   text = paste0("ID: ", ID,
                                                 "<br>Description: ", Description,
                                                 "<br>GeneCount: ", .data[[x_var]],
                                                 "<br>p.adjust: ", signif(p.adjust, 3))))}
      
      if (input$plotType == "dotplot") {
        p <- p + geom_point(aes(color = p.adjust, size = .data[[x_var]])) +
          scale_color_gradient(low = input$colorLow, high = input$colorHigh) +
          labs(x = x_var, y = NULL, color = "p.adjust", size = x_var, title = "ORA Dot Plot")
      } else {
        p <- p + geom_bar(aes(fill = p.adjust), stat = "identity") +
          scale_fill_gradient(low = input$colorLow, high = input$colorHigh) +
          labs(x = x_var, y = NULL, fill = "p.adjust", title = "ORA Bar Plot")}
      
      p + theme_minimal() +
        theme(axis.text.y  = element_text(size = text_size),
              axis.text.x  = element_text(size = text_size),
              axis.title.x = element_text(size = text_size + 1),
              legend.text  = element_text(size = text_size - 3),
              plot.title   = element_text(size = text_size + 2, face = "bold", hjust = 0.5))})
    
    dynamic_plot_height <- reactive({
      req(filtered_enrichment())
      n_terms <- min(input$nPlotTerms, nrow(filtered_enrichment()))
      calculated_height <- 200 + (n_terms * 25)
      max(600, calculated_height)})
    
    output$enrichmentDotPlotUI <- renderUI({
      req(dynamic_plot_height())
      height_px <- paste0(dynamic_plot_height(), "px")
      plotlyOutput(ns("enrichmentDotPlot"), height = height_px)})
    
    output$enrichmentDotPlot <- renderPlotly({
      req(plot_obj())
      ggplotly(plot_obj(), tooltip = "text")})
    
    output$downloadPlotUI3 <- renderUI({
      req(enrichment_results())
      dropdownButton(
        label = "Download Plot",
        icon = icon("download"),
        circle = FALSE,
        status = "secondary",
        width = "300px",
        tags$h4("Download Format"),
        downloadButton(ns("downloadPlotHTML"), "HTML", class = "btn-block"),
        downloadButton(ns("downloadPlotPNG"), "PNG", class = "btn-block"),
        downloadButton(ns("downloadPlotJPG"), "JPG", class = "btn-block"))})
    
    # DOWNLOAD HTML
    output$downloadPlotHTML <- downloadHandler(
      filename = function() {
        paste("MetaRank_Enrich_plot_", Sys.Date(), ".html", sep = "")},
      content = function(file) {
        h_px <- dynamic_plot_height()
        g <- ggplotly(plot_obj(), height = h_px) 
        htmlwidgets::saveWidget(g, file, selfcontained = TRUE)})
    
    # DOWNLOAD PNG (Ajustado dinámicamente)
    output$downloadPlotPNG <- downloadHandler(
      filename = function() {
        paste("MetaRank_Enrich_plot_", Sys.Date(), ".png", sep = "")},
      content = function(file) {
        req(dynamic_plot_height())
        h_in <- dynamic_plot_height() / 96 
        ggsave(file, plot = plot_obj(), device = "png",
               width = 14.5, height = h_in, units = "in", dpi = 96)})
    
    # DOWNLOAD JPG (Ajustado dinámicamente)
    output$downloadPlotJPG <- downloadHandler(
      filename = function() {
        paste("MetaRank_Enrich_plot_", Sys.Date(), ".jpg", sep = "")},
      content = function(file) {
        req(dynamic_plot_height())
        h_in <- dynamic_plot_height() / 96
        ggsave(file, plot = plot_obj(), device = "jpeg",
               width = 14.5, height = h_in, units = "in", dpi = 96)})
    
    
    
    #----------------------------------------------------------------------------
    # NAVIGATION EXTRAS (INFORMATION BUTTONS + DATA UPDATE)
    #----------------------------------------------------------------------------
    
    # HIDE BUTTONS BEFORE THE ANALYSIS
    output$Hidebuttons <- reactive({
      !is.null(combined_results())})
    outputOptions(output, "Hidebuttons", suspendWhenHidden = FALSE)
    
    
    # --------------------------------------------------------------------------
    # 1. EXAMPLE OF FILE FORMAT INFORMATION (RP - Weighted)
    # --------------------------------------------------------------------------
    observeEvent(input$showExample, {
      showModal(
        modalDialog(
          title = tags$div(icon("file-invoice"), " Example File Contents (RP)", 
                           style = "text-align: center; font-weight: bold; font-size: 20px; color: #2980b9;"),
          easy_close = TRUE,
          size = "l",
          
          # Code Blocks for CSV and TSV
          fluidRow(
            column(6, 
                   tags$h5(icon("table"), " Tab-Separated (.tsv/.txt)", style="color: #34495e; font-weight:bold; margin-top:0;"),
                   HTML("<pre style='background-color:#f8f9fa; border: 1px solid #e9ecef; padding:12px; border-radius:5px; font-family: monospace; color: #2c3e50;'>COL10A1\t12.01\nCOL11A1\t11.74\nGREM1\t10.74\nMMP1\t10.11\nMMP12\t9.29\nSPINK1\t8.97\n...</pre>")),
            column(6, 
                   tags$h5(icon("file-csv"), " Comma-Separated (.csv)", style="color: #34495e; font-weight:bold; margin-top:0;"),
                   HTML("<pre style='background-color:#f8f9fa; border: 1px solid #e9ecef; padding:12px; border-radius:5px; font-family: monospace; color: #2c3e50;'>COL10A1,12.01\nCOL11A1,11.74\nGREM1,10.74\nMMP1,10.11\nMMP12,9.29\nSPINK1,8.97\n...</pre>"))),
          
          # Info Box
          tags$div(
            style = "background-color: #eef6fb; padding: 15px; border-radius: 5px; border-left: 5px solid #2980b9; margin-top: 15px;",
            tags$h5(icon("info-circle"), " Important Notes", style="color: #2980b9; font-weight:bold; margin-top:0;"),
            tags$ul(
              style = "margin-bottom: 0; padding-left: 20px;",
              tags$li(HTML("<b>Accepted formats:</b> <code>.tsv</code>, <code>.csv</code>, and <code>.txt</code>."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>Structure:</b> Two columns required (Gene Name & Statistic)."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>Decimal format:</b> Use a <code>dot (.)</code>. Commas are not allowed as decimal separators."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>No headers:</b> Do <u>not</u> include any header row. Input should start directly with data."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>File Size Limit:</b> Maximum <b>30 MB</b> per file.")))),
          footer = modalButton("Close")))})
    
    
    # --------------------------------------------------------------------------
    # 2. EXAMPLE OF FILE FORMAT INFORMATION RRA (Unweighted)
    # --------------------------------------------------------------------------
    observeEvent(input$showExampleRRA, {
      showModal(
        modalDialog(
          title = tags$div(icon("file-alt"), " Example File Contents (RRA)", 
                           style = "text-align: center; font-weight: bold; font-size: 20px; color: #2980b9;"),
          easy_close = TRUE,
          size = "m",
          
          # Code Block
          tags$div(
            style = "margin-bottom: 15px;",
            tags$h5(icon("code"), " Text Format (.txt)", style="color: #34495e; font-weight:bold; margin-top:0;"),
            HTML("<pre style='background-color:#f8f9fa; border: 1px solid #e9ecef; padding:12px; border-radius:5px; font-family: monospace; color: #2c3e50;'>BRCA1\nEGFR\nMYC\nPTEN\nTP53\n...</pre>")),
          
          # Info Box
          tags$div(
            style = "background-color: #eef6fb; padding: 15px; border-radius: 5px; border-left: 5px solid #2980b9;",
            tags$h5(icon("info-circle"), " Important Notes", style="color: #2980b9; font-weight:bold; margin-top:0;"),
            tags$ul(
              style = "margin-bottom: 0; padding-left: 20px;",
              tags$li(HTML("<b>Accepted format:</b> Only <code>.txt</code> files are supported."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>One gene per line:</b> Each row must contain exactly one gene name."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>Ranked input:</b> Lists must be <u>pre-ranked</u> (most significant genes at the top)."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>No headers:</b> Do <u>not</u> include any header row."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>File Size Limit:</b> Maximum <b>30 MB</b> per file.")))),
          footer = modalButton("Close")))})
    
    
    # --------------------------------------------------------------------------
    # 3. EXAMPLE OF ORIGIN MANAGEMENT (Origin Info)
    # --------------------------------------------------------------------------
    observeEvent(input$showOriginInfo, {
      showModal(
        modalDialog(
          title = tags$div(icon("sitemap"), " Origin Vector Information", 
                           style = "text-align: center; font-weight: bold; font-size: 20px; color: #2980b9;"),
          size = "l",
          easy_close = TRUE,
          
          # Image
          div(style = "text-align: center; margin-bottom: 20px;",
              img(src = "org.png", height = "auto", width = "100%", style="max-width: 500px; border: 1px solid #eee; padding: 5px; border-radius: 5px;")),
          
          # Info Box
          tags$div(
            style = "background-color: #eef6fb; padding: 15px; border-radius: 5px; border-left: 5px solid #2980b9;",
            tags$h5(icon("info-circle"), " Configuration Rules", style="color: #2980b9; font-weight:bold; margin-top:0;"),
            tags$ul(
              style = "margin-bottom: 0; padding-left: 20px;",
              tags$li(HTML("<b>Purpose:</b> Group your gene lists by metadata (e.g., <i>lab, technology, date</i>)."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>Format:</b> Comma-separated numeric vector (e.g., <code>1,2,1,2,2</code>)."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>Length:</b> Must strictly match the number of uploaded lists."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>Replication:</b> Each group must contain at least two lists."), style="margin-bottom: 5px;"),
              tags$li(HTML("<b>Example Data:</b> For the provided example datasets, use <code>1,2,1,2</code>.")))),
          footer = modalButton("Close")))})
    
    
    # --------------------------------------------------------------------------
    # 4. EXAMPLE OF PERSONALIZED ANNOTATIONS (Annotations Info)
    # --------------------------------------------------------------------------
    observeEvent(input$showPersonalizedInfo, {
      showModal(
        modalDialog(
          title = tags$div(icon("database"), " Personalized Annotations", 
                           style = "text-align: center; font-weight: bold; font-size: 20px; color: #2980b9;"),
          easy_close = TRUE,
          size = "l", 
          
          # Code Block
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
          
          # Info Box
          tags$div(
            style = "background-color: #eef6fb; padding: 15px; border-radius: 5px; border-left: 5px solid #2980b9; margin-bottom: 20px;",
            tags$h5(icon("clipboard-check"), " File Requirements", style="color: #2980b9; font-weight:bold; margin-top:0;"),
            tags$ul(
              style = "margin-bottom: 0; padding-left: 20px;",
              tags$li(HTML("<b>Single File:</b> Upload exactly one annotation file."), style="margin-bottom: 6px;"),
              tags$li(HTML("<b>Header Required:</b> The file <b>must</b> include a header row (e.g., <code>FeatureID\\tAnnotationValue</code>)."), style="margin-bottom: 6px;"),
              tags$li(HTML("<b>Columns:</b> 1st Column = Features (Genes/IDs), 2nd Column = Annotation (Location/Function)."), style="margin-bottom: 6px;"),
              tags$li(HTML("<b>Separators:</b> Tab-separated (<code>.tsv</code>, <code>.txt</code>) or Comma-separated (<code>.csv</code>).")))),
          
          # Critical Warning Box
          tags$div(
            style = "background-color: #fff5f5; padding: 15px; border-radius: 5px; border-left: 5px solid #c0392b;",
            tags$h5(icon("exclamation-triangle"), " CRITICAL WARNING: Nomenclature", style="color: #c0392b; font-weight:bold; margin-top:0;"),
            tags$p(HTML("Ensure the feature names (IDs) in this annotation file <b>exactly match</b> the nomenclature (case-sensitivity, version) used in your input lists to avoid recognition issues."), 
                   style = "color: #7f2a21; margin-bottom: 0;")),
          footer = modalButton("Close")))})
    
    
    
    # SHOW/HIDE EXCLUDED GENES TABLE 
    observeEvent(input$toggleExcludedGenes, {
      toggle("excludedGenesSection")})
    
    
    # SLIDER AUTO-UPDATE 
    observeEvent({
      input$files
      input$useExampleData
      input$pasteText
      input$pasteText1
      input$files_rra
      input$useExampleData_rra
      input$pasteText_rra}, 
      
      {total_lists <- NULL
      if (input$selectedPackage == "Rankprodpack") {
        if (input$useExampleData && !is.null(example_files)) {
          total_lists <- length(example_files)
        } else if (input$inputMethod == "files" && !is.null(input$files)) {
          total_lists <- length(input$files$datapath)
        } else if (input$inputMethod == "paste") {
          if (input$pasteFormat == "csv" && !is.null(input$pasteText) && nzchar(input$pasteText)) {
            total_lists <- length(strsplit(input$pasteText, "###")[[1]])
          } else if (input$pasteFormat == "tsv" && !is.null(input$pasteText1) && nzchar(input$pasteText1)) {
            total_lists <- length(strsplit(input$pasteText1, "###")[[1]])}}
        
        if (!is.null(total_lists) && total_lists > 0) {
          updateSliderInput(session, "gene_app_filter_rp",
                            min = 1,
                            max = total_lists,
                            step = 1,
                            value = total_lists)}
        
      } else if (input$selectedPackage == "RRApack") {
        if (input$useExampleData_rra && !is.null(example_files_rra)) {
          total_lists <- length(example_files_rra)
        } else if (input$inputMethod == "files" && !is.null(input$files_rra)) {
          total_lists <- length(input$files_rra$datapath)
        } else if (input$inputMethod == "paste" && !is.null(input$pasteText_rra) && nzchar(input$pasteText_rra)) {
          total_lists <- length(strsplit(input$pasteText_rra, "###")[[1]])}
        
        if (!is.null(total_lists) && total_lists > 0) {
          updateSliderInput(session, "gene_app_filter_rra",
                            min = 1,
                            max = total_lists,
                            step = 1,
                            value = total_lists)}}})
    
    selected_gene_app_filter <- reactive({
      if (input$selectedPackage == "Rankprodpack") {
        input$gene_app_filter_rp
      } else {
        input$gene_app_filter_rra}})
    
    
    # RANKPRODPACK EXAMPLE ZIP
    output$downloadExampleRank <- downloadHandler(
      filename = function() { "MetaRank_example_Rankprodpack.zip" },
      content = function(zipfile) {
        tmp <- tempdir()
        files_cp <- file.path(tmp, basename(example_files))
        file.copy(example_files, files_cp, overwrite = TRUE)
        owd <- setwd(tmp); on.exit(setwd(owd), add=TRUE)
        utils::zip(zipfile, basename(files_cp), extras = "-j")},
      contentType = "application/zip")
    
    # ROBUSTRANKAGGREGPACK EXAMPLE ZIP
    output$downloadExampleRRA <- downloadHandler(
      filename = function() { "MetaRank_example_RRApack.zip" },
      content = function(zipfile) {
        tmp <- tempdir()
        files_cp <- file.path(tmp, basename(example_files_rra))
        file.copy(example_files_rra, files_cp, overwrite = TRUE)
        owd <- setwd(tmp); on.exit(setwd(owd), add=TRUE)
        utils::zip(zipfile, basename(files_cp), extras = "-j")},
      contentType = "application/zip")
    
    prev_input_method <- reactiveVal("files")
    
    observeEvent(input$inputMethod, {
      prev <- prev_input_method()
      
      if (input$inputMethod == "paste" && prev == "files") {
        shinyjs::reset("files")
        shinyjs::reset("files_rra")
        
        updateTextInput(session, "pasteText", value = "")
        updateTextInput(session, "pasteText1", value = "")
        updateTextInput(session, "pasteText_rra", value = "")
        updateSwitchInput(session, "useExampleData", value = FALSE)
        updateSwitchInput(session, "useExampleData_rra", value = FALSE)
        
      } else if (input$inputMethod == "files" && prev == "paste") {
        updateTextInput(session, "pasteText", value = "")
        updateTextInput(session, "pasteText1", value = "")
        updateTextInput(session, "pasteText_rra", value = "")}
      
      prev_input_method(input$inputMethod)})
    
    
    # Descarga de archivo de anotación de ejemplo
    output$downloadPersonalizedExample <- downloadHandler(
      filename = function() {
        "protein_location_example.txt"},
      content = function(file) {
        src_path <- "./example_data/protein_location.txt"
        if (!file.exists(src_path)) {
          stop("Error: El archivo de ejemplo no se encontró en la ruta especificada: ", src_path)}
        
        file.copy(src_path, file)})
    
  })}