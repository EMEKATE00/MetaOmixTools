###############################################################################
#                                                                             #
#           MetaOmixTools: suite of applications designed to integrate        #
#           and analyze omics data.                                           #
#             - MetaRank: Meta-analysis of ranked gene lists                  #
#             - MetaEnrichGO: Functional enrichment analysis                  #
#                                                                             #
#           Authors:                                                          #
#             - Maksym Kupchyk Tiurin,                                        #
#             - Francisco Javier Cordero Felipe,                              #
#             - Francisco García García (Supervisor),                         #
#             - Rubén Grillo Risco                                            #
#                                                                             #
###############################################################################

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
library(bslib)
library(R.utils)
library(markdown)

#------------------------------------------------------------------------------- 
# The user interface (ui) object controls the layout and appearance of your app
#-------------------------------------------------------------------------------

source("./modules/metaenrichgo/mod_metaenrichgo.R")
source("./modules/metarank/mod_metarank.R")
addResourcePath("style.css", "modules/metarank/www")
addResourcePath("style.css", "modules/metaenrichgo/www")

#------------------------------------------------------------------------------
# Define the main UI of the application using a tag list.
#------------------------------------------------------------------------------

ui <- tagList(
  useShinyjs(),  # Enable shinyjs for advanced JavaScript functionality within Shiny

  # Head section: includes custom CSS styling and JavaScript functions.
  tags$head(
    includeCSS("www/style.css"),
    tags$style(HTML("
      /* Global styles: reset and basic configuration */
      html, body {
        font-size: 16px;
        margin: 0;
        padding: 0;
        height: 100%;
        overflow-x: hidden;
        scroll-behavior: smooth;
      }
      
      /* Description section styling */
      .description-section {
        padding: 30px;
        background-color: #fdfdfd;
      }
      .description-row {
        display: flex;
        flex-wrap: wrap;
        align-items: stretch; 
      }
      .info-box, .citation-box {
        padding: 1em;
        margin-bottom: 1em;
        flex: 1; 
        border: 2px solid #000000;
        background-color: #ffffff;
      }
      
      /* Divider section: horizontal line and explanatory text */
      .divider-section {
        text-align: center;
        padding: 20px;
        background-color: #f5f5f5;
      }
      .divider-section hr {
        margin: 0 auto;
        width: 80%;
        border-top: 1px solid #000;
      }
      .divider-section p {
        margin: 10px 0 0;
        font-weight: bold;
      }
      
      /* Hero section: tool logos with interactive hover effects */
      .hero-section {
        position: relative;
        height: 100vh; 
        display: flex;
        flex-direction: row;
        overflow: hidden;
      }
      .hero-col {
        flex: 1;
        position: relative;
        overflow: hidden;
      }
      .hero-img {
        width: 100%;
        height: 100%;
        object-fit: cover;
        transition: transform 0.2s ease-out;
      }
      /* Overlay styling for hover effect */
      .hero-col .overlay {
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background: rgba(0, 0, 0, 0.8);
        opacity: 0;
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        color: #fff;
        text-align: center;
        transition: opacity 0.3s ease;
        padding: 20px;
      }
      .hero-col:hover .overlay {
        opacity: 100%;
      }
      .overlay h3 {
        margin-bottom: 10px;
      }
      .overlay button {
        margin-top: 10px;
      }
    ")),
    

    tags$script(HTML("
      function navigateToPanel(panel) {
        Shiny.setInputValue('navigateToPanel', panel, {priority: 'event'});
      }
    ")),
    
    
    
    tags$head(
      # CSS de GLightbox
      tags$link(rel = "stylesheet", href = "https://cdn.jsdelivr.net/npm/glightbox/dist/css/glightbox.min.css"),
      # JS de GLightbox
      tags$script(src = "https://cdn.jsdelivr.net/npm/glightbox/dist/js/glightbox.min.js"),
      # Inicializador
      tags$script(HTML("
    document.addEventListener('DOMContentLoaded', function() {
      const lightbox = GLightbox({ selector: '.glightbox' });
    });
  "))
    )
    
  ),
  
  
  
 
  page_navbar(
    title = "MetaOmixTools",
    id = "mot_navbar",
    
    # Home panel ("Home") with description and citation.
    nav_panel("Home", card(HTML('
      <!-- DESCRIPTION SECTION: Title, description, and citation information -->
      <section class="description-section">
        <div class="container">
          <h1 class="text-left">Welcome to MetaOmixTools!</h1>
          <div class="row description-row" style="margin-top: 20px;">
            <div class="col-md-8">
              <div class="info-box">
                <h3>Suite Description</h3>
                <p style="text-align:justify;">
                  MetaOmixTools is a suite of applications designed to integrate and analyze omics data.
                  In research, thousands of data points generate ranked gene lists that often vary across studies.
                  To address this, MetaOmixTools includes:
                  <strong>MetaRank</strong>, a web tool developed in R with Shiny for the meta-analysis of ranked gene lists,
                  and <strong>MetaEnrichGO</strong>, which performs gene ontology enrichment analysis.
                  <br><br>
                  <strong>MetaRank</strong> accepts input as tables in TSV/CSV formats or via direct copy-paste and produces
                  consensus gene rankings. Its outputs include tables (TSV, CSV) and figures (PNG, HTML).
                  <br><br>
                  <strong>MetaEnrichGO</strong> accepts input as TXT files or pasted data and outputs include figures (PNG)
                  and tables (TSV, CSV).
                  This integrated approach facilitates the identification of key biomarkers and their functional roles in biological processes.
                </p>
              </div>
            </div>
            <div class="col-md-4">
              <div class="citation-box">
                <h4>How to cite:</h4>
                <p>
                  Maksym Kupchyk Tiurin, Francisco Javier Cordero Felipe, Francisco García García, and Rubén Grillo Risco. 
                  <strong>MetaOmixTools:</strong> A suite for integrative meta-omics analysis combining gene ranking and functional enrichment. 
                  <i>Bioinformatics Journal</i>, 2025. doi:10.1234/bioinformatics.2025.12345
                </p>
                <h5 style="margin-top: 38px;">Funding:</h5>
                <p>
                  <a href="https://www.cipf.es/" target="_blank">
                    <img src="images/5.png" style="width:200px; align-self: auto;">
                  </a>
                </p>
              </div>
            </div>
          </div>
        </div>
      </section>
      
      <!-- DIVIDER SECTION: Horizontal line and clarifying text -->
      <section class="divider-section">
        <hr>
        <p>Discover the tools offered by MetaOmixTools:</p>
      </section>
      
      <!-- HERO SECTION: Displays tool logos with hover effects and additional tool details -->
      <section class="hero-section">
        <!-- Column for MetaRank -->
        <div class="hero-col">
          <a href="metarank">
            <img src="images/1.png" alt="MetaRank" class="hero-img">
          </a>
          <div class="overlay">
            <h3>MetaRank</h3>
            <p>
              Performs advanced statistical analyses by combining ranked gene lists.<br>
              <strong>Inputs:</strong> TSV/CSV files or direct copy-paste.<br>
              <strong>Outputs:</strong> TSV/CSV tables and figures in PNG/HTML formats.
            </p>
            <button id="btnMetaRank" class="btn btn-info" onclick="navigateToPanel(\'MetaRank\')">Try it!</button>
          </div>
        </div>
        <!-- Column for MetaEnrichGO -->
        <div class="hero-col">
          <a href="metaenrichgo">
            <img src="images/2.png" alt="MetaEnrichGO" class="hero-img">
          </a>
          <div class="overlay">
            <h3>MetaEnrichGO</h3>
            <p>
              Conducts functional enrichment analysis using Gene Ontology.<br>
              <strong>Inputs:</strong> TXT files or direct copy-paste.<br>
              <strong>Outputs:</strong> PNG figures and tables in TSV/CSV formats.
            </p>
            <button id="btnMetaEnrichGO" class="btn btn-info" onclick="navigateToPanel(\'MetaEnrichGO\')">Try it!</button>
          </div>
        </div>
      </section>
    '))),
    
    # MetaRank panel: dynamically loads the MetaRank module's UI.
    tabPanel("MetaRank", mod_metaRank_ui("metarank1")),
    
    # MetaEnrichGO panel: dynamically loads the MetaEnrichGO module's UI.
    tabPanel("MetaEnrichGO", mod_metaEnrichGO_ui("mega1")),
    
    # "About" dropdown menu with additional information tabs.
    navbarMenu("About",
               tabPanel("Instructions",
                        navlistPanel(
                          id    = ("instructions_nav"),  
                          widths = c(3, 9),   
                          
                          "Instructions",
                          br(),
                          tabPanel(title = "Overview",
                                   tags$iframe(
                                     src    = "overview_readme.html",
                                     width  = "100%", 
                                     height = "900px",
                                     style  = "border:none;")),
                          
                          tabPanel(title = "MetaRank",
                                   tags$iframe(
                                     src    = "metarank_readme.html",
                                     width  = "100%", 
                                     height = "900px",
                                     style  = "border:none;")),
                          
                          tabPanel(title = "MetaRank Example Usage",
                                   tags$iframe(
                                     src    = "metarank_example.html",
                                     width  = "100%", 
                                     height = "900px",
                                     style  = "border:none;")),
                          
                          tabPanel(title = "MetaEnrichGO",
                                   tags$iframe(
                                     src    = "metaenrichgo_readme.html",
                                     width  = "100%", 
                                     height = "900px",
                                     style  = "border:none;")),
                          tabPanel(title = "MetaEnrichGO Example Usage",
                                   tags$iframe(
                                     src    = "metaenrichgo_example.html",
                                     width  = "100%", 
                                     height = "900px",
                                     style  = "border:none;")))),
               
               tabPanel("Contact",
                        navlistPanel(
                          id    = ("instructions_nav"),  
                          widths = c(3, 9), 
                          
                          "Help Section",
                          br(),
                          tabPanel(title = "Contact us",
                                   tags$iframe(
                                     src    = "contact.html",
                                     width  = "100%", 
                                     height = "900px",
                                     style  = "border:none;")),
                          
                          tabPanel(title = "Questions & Answers",
                                   tags$iframe(
                                     src    = "FAQ_readme.html",
                                     width  = "100%", 
                                     height = "900px",
                                     style  = "border:none;")))),
                        
               tabPanel("Gallery",
                        navlistPanel(
                          id = "instructions_nav",
                          widths = c(3, 9),
                          
                          "Instructions",
                          br(),
                          
                          tabPanel(title = "Overview",
                                   fluidRow(column(6, tags$a(href = "images/metaomixtools_overview.png", 
                                              class = "glightbox", 
                                              `data-gallery` = "gallery-overview",
                                              tags$img(src = "images/metaomixtools_overview.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                            tags$p("General overview of the MetaOmixTools suite.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                           column(6,tags$a(
                                                    href = "images/metaomixtools_apps.png", 
                                                    class = "glightbox", 
                                                    `data-gallery` = "gallery-overview",
                                                    tags$img(src = "images/metaomixtools_apps.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                                  tags$p("General overview of the apps that MetaOmixTools offers.", style = "text-align: center; font-style: italic; margin-top: 5px;")))),
                          
                        tabPanel(title = "MetaRank",
                                 fluidRow(
                                   column(6, tags$a(
                                            href = "images/metarank_overview.png", 
                                            class = "glightbox", 
                                            `data-gallery` = "gallery-metarank",
                                            tags$img(
                                              src = "images/metarank_overview.png", 
                                              width = "100%", 
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 250px; object-fit: contain; background-color: white; display: block; margin: auto;")),
                                          tags$p("General overview of MetaRank module.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                   
                                   column(6, tags$a(
                                            href = "images/metarank_rankprod_input.png",
                                            class = "glightbox",
                                            `data-gallery` = "gallery-metarank",
                                            tags$img(
                                              src = "images/metarank_rankprod_input.png",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 250px; object-fit: contain; background-color: white; display: block; margin: auto;")),
                                          tags$p("RankProd input panel.", style = "text-align: center; font-style: italic; margin-top: 5px;"))),
                                 
                                 fluidRow(
                                   column(6, tags$a(
                                            href = "images/metarank_rra_input.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                            tags$img(src = "images/metarank_rra_input.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                          tags$p("Input interface for RRA analysis.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                   
                                   column(6,
                                          tags$a(
                                            href = "images/metarank_rankprod_info.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                            tags$img(src = "images/metarank_rankprod_info.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                          tags$p("Information box for RankProd mode.", style = "text-align: center; font-style: italic; margin-top: 5px;"))),
                                 
                                 fluidRow(
                                   column(6,
                                          tags$a(
                                            href = "images/metarank_rra_info.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                            tags$img(src = "images/metarank_rra_info.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                          tags$p("Information box for RRA mode.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                   
                                   column(6,
                                          tags$a(
                                            href = "images/metarank_rankprod_parameters.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                            tags$img(src = "images/metarank_rankprod_parameters.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                          tags$p("RankProd parameter settings.", style = "text-align: center; font-style: italic; margin-top: 5px;"))),
                                 
                                 fluidRow(
                                   column(6,
                                          tags$a(
                                            href = "images/metarank_rra_parameters.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                            tags$img(src = "images/metarank_rra_parameters.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("RRA parameter settings.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(
                                            href = "images/metarank_rankprod_origin.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                            tags$img(src = "images/metarank_rankprod_origin.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Group origin definition for RankProd.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(
                                            href = "images/metarank_enrich.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                            tags$img(src = "images/metarank_enrich.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Enrichment options after ranking.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(
                                            href = "images/metarank_data_1.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                            tags$img(src = "images/metarank_data_1.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Example data configuration (1).", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/metarank_data_2.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_data_2.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Example data configuration (2).", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/metarank_data_3.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_data_3.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Example data configuration (3).", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/metarank_data_4.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_data_4.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Example data configuration (4).", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/metarank_data_5.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_data_5.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Example data configuration (5).", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/metarank_error_1.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_error_1.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Error example 1.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/metarank_error_2.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_error_2.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Error example 2.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/metarank_error_3.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_error_3.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Error example 3.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/metarank_error_4.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_error_4.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Error example 4.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/metarank_error_5.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_error_5.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Error example 5.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/metarank_error_6.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_error_6.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Error example 6.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/metarank_error_7.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_error_7.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Error example 7.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/metarank_error_8.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_error_8.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Error example 8.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/metarank_error_9.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_error_9.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Error example 9.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/metarank_workflow.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_workflow.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("MetaRank full workflow.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/metarank_upsetplot.jpg", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/metarank_upsetplot.jpg", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("UpSet plot of example dataset.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/example2_metarank_upset.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/example2_metarank_upset.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("UpSet plot in Example 2.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/example1_metarank_heatmap.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/example1_metarank_heatmap.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Heatmap from Example 1.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/example2_metarank_heatmap.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/example2_metarank_heatmap.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Heatmap from Example 2.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   )
                                 ),
                                 fluidRow(
                                   column(6,
                                          tags$a(href = "images/example3_metarank_enrich.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                          tags$img(src = "images/example3_metarank_enrich.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")
                                          ),
                                          tags$p("Enrichment plot from Example 3.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                   ),
                                   column(6,
                                          tags$a(href = "images/example2_metarank_enrich.png", class = "glightbox", `data-gallery` = "gallery-metarank", tags$img(src = "images/example2_metarank_enrich.png", 
                                          width = "100%", 
                                          style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                          tags$p("Enrichment plot from Example 2.", style = "text-align: center; font-style: italic; margin-top: 5px;")))),
                          
                        tabPanel(
                          title = "MetaEnrichGO",
                          fluidRow(
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_overview.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_overview.png", 
                                              width = "100%", 
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("General overview of MetaEnrichGO.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_input.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_input.png", 
                                              width = "100%", 
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Input via file upload or text paste.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_i.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_i.png",
                                              width = "100%", 
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Contextual information menu (i).", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_parameters.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_parameters.png",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Functional analysis parameters.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_resultstable.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_resultstable.png",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Combined results table.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_dot.jpg",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_dot.jpg",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Dotplot of enriched GO terms.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_bar.jpg",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_bar.jpg",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Barplot of enriched GO terms.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_workflow.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_workflow.png",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Workflow diagram of the analysis.", style = "text-align: center; font-style: italic; margin-top: 5px;")),

                            column(6, tags$a(
                                     href = "images/metaenrichgo_file_error.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_file_error.png",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Input file format error.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_paste_error.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_paste_error.png",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Text paste input error.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_organism_error.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_organism_error.png",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Organism selection error.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_geneid_error.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_geneid_error.png",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Gene identifier error.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                            
                            column(6, tags$a(
                                     href = "images/metaenrichgo_noinfo_error.png",
                                     class = "glightbox",
                                     `data-gallery` = "gallery-metaenrichgo",
                                     tags$img(src = "images/metaenrichgo_noinfo_error.png",
                                              width = "100%",
                                              style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                   tags$p("Error: no functional information available.", style = "text-align: center; font-style: italic; margin-top: 5px;")))))),
               
               tabPanel("GitHub",
                        tags$div(style = "padding: 25px; font-family: 'Segoe UI', sans-serif; line-height: 1.7; max-width: 900px; margin: auto;",
                                 
                                 tags$h3("MetaOmixTools Suite & Individual Tools", style = "margin-bottom: 20px;"),
                                 
                                 tags$h4("Online Tools", style = "margin-top: 30px;"),
                                 tags$p("You can launch each application individually:"),
                                 tags$ul(
                                   tags$li(
                                     tags$b("MetaOmixTools Suite: "),
                                     tags$a(href = "https://maksymkupchyktiurin00.shinyapps.io/MetaOmixTools/", 
                                            "Launch App", target = "_blank", style = "color: #007ACC; text-decoration: none;")
                                   ),
                                   tags$li(
                                     tags$b("MetaRank: "),
                                     tags$a(href = "https://maksymkupchyktiurin00.shinyapps.io/MetaRank/", 
                                            "Launch App", target = "_blank", style = "color: #007ACC; text-decoration: none;")
                                   ),
                                   tags$li(
                                     tags$b("MetaEnrichGO: "),
                                     tags$a(href = "https://maksymkupchyktiurin00.shinyapps.io/MetaEnrichGO/", 
                                            "Launch App", target = "_blank", style = "color: #007ACC; text-decoration: none;")
                                   )
                                 ),
                                 
                                 tags$h4("GitHub Repositories", style = "margin-top: 30px;"),
                                 tags$p("All apps are open-source and available on GitHub:"),
                                 tags$ul(
                                   tags$li(tags$b("MetaOmixTools Suite: "), 
                                           tags$a(href = "https://github.com/your_user/MetaOmixTools", 
                                                  "GitHub Repo", target = "_blank")),
                                   tags$li(tags$b("MetaRank: "), 
                                           tags$a(href = "https://github.com/your_user/MetaRank", 
                                                  "GitHub Repo", target = "_blank")),
                                   tags$li(tags$b("MetaEnrichGO: "), 
                                           tags$a(href = "https://github.com/your_user/MetaEnrichGO", 
                                                  "GitHub Repo", target = "_blank"))
                                 ),
                                 
                                 tags$h4("Suite Architecture", style = "margin-top: 30px;"),
                                 tags$p("The suite is built using modular Shiny components. Each tool is encapsulated as a reusable module that can be independently loaded, improving maintainability, flexibility, and integration into larger workflows."),
                                 
                                 tags$h4("Web Tool Design Principles", style = "margin-top: 30px;"),
                                 tags$p("The development of MetaOmixTools adheres to best practices in web-based scientific tools:"),
                                 tags$ul(
                                   tags$li(tags$b("Reproducibility:"), " All steps and computations are transparent and traceable."),
                                   tags$li(tags$b("Modularity:"), " Tools are implemented as independent Shiny modules."),
                                   tags$li(tags$b("Clarity:"), " The UI is designed to be minimalistic, responsive, and intuitive."),
                                   tags$li(tags$b("Accessibility:"), " No installation required. Hosted online for direct usage."),
                                   tags$li( tags$b("Maintainability:"), " Clear file structure and logical code separation for ease of updates."),
                                   tags$li(tags$b("Documentation:"), " Every module and function is documented for both users and developers.")
                                 ),
                                 
                                 tags$div(style = "margin-top: 40px;",
                                          tags$a(href = "https://github.com/your_user/MetaOmixTools", target = "_blank",
                                                 class = "btn btn-primary", 
                                                 style = "color: white; background-color: #24292e; border: none; padding: 10px 18px; font-size: 16px; border-radius: 6px;",
                                                 shiny::icon("github"), " View Suite on GitHub")
                                 )
                        )
               )
               
               
               ),
    
    bslib::nav_spacer(),
  
    
    header = tagList(
      tags$head(
        tags$style(HTML("
        .navbar-custom-right {
          position: absolute;
          right: 50px;
          top: 15px;
        }
        .navbar-custom-right a {
          color: white;
          margin-left: 20px;
          font-size: 20px;
        }
      "))),
      
      tags$div(class = "navbar-custom-right",
               actionLink("go_contact", label = NULL, icon = icon("envelope"), class = "shiny-link-icon"),
               tags$a(href = "https://github.com/tu_usuario/tu_repo", icon("github"), target = "_blank"))
      
      
      
      
      )))



#------------------------------------------------------------------------------
# Define the server logic for the main application.
#------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  mod_metaRank_server("metarank1")
  mod_metaEnrichGO_server("mega1")
  
  observeEvent(input$navigateToPanel, {
    updateNavbarPage(session, "mot_navbar", selected = input$navigateToPanel)})

  observeEvent(input$go_contact, {
    updateTabsetPanel(session, inputId = "mot_navbar", selected = "Contact")
  })
  

}

#------------------------------------------------------------------------------
# Launch the Shiny application with the defined UI and server logic.
#------------------------------------------------------------------------------

shinyApp(ui = ui, server = server)
