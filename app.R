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

library(shiny)                  # Main framework for building interactive web applications.
library(shinycssloaders)        # Adds loading animations (spinners) to outputs while calculating.
library(shinyWidgets)           # Provides advanced and custom input controls/UI components.
library(shinyBS)                # Adds extra Bootstrap components (tooltips, popovers, modals).
library(shinyjs)                # Enables common JavaScript operations from R (hide, show, reset).
library(bslib)                  # Tools for customizing the visual theme (Bootstrap versions, colors).
library(dplyr)                  # Fundamental tool for data manipulation and transformation (data frames).
library(DT)                     # Creates interactive data tables with search, sort, and pagination.
library(ggplot2)                # Standard system for declaratively creating high-quality static graphics.
library(processx)               # Allows running system processes in the background.
library(plotly)                 # Converts static graphs (ggplot2) into interactive ones (zoom, hover).
library(R.utils)                # diverse utility functions (file handling, paths, arguments).
library(markdown)               # Renders text written in Markdown format into HTML within the app.

#------------------------------------------------------------------------------- 
# Load Modules
#-------------------------------------------------------------------------------

options(shiny.maxRequestSize = 30 * 1024^2)

source("./modules/metaenrichgo/mod_metaenrichgo.R")
source("./modules/metarank/mod_metarank.R")

#------------------------------------------------------------------------------
# Define the main UI of the application
#------------------------------------------------------------------------------

ui <- tagList(
  useShinyjs(),
  
  # --- HEAD: CSS, JS, Fuentes y GLightbox ---
  tags$head(
    includeCSS("www/style.css"),
    
    # Google Font (Lato)
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = "anonymous"),
    tags$link(href = "https://fonts.googleapis.com/css2?family=Lato:wght@300;400;700&display=swap", rel = "stylesheet"),
    
    # CSS de GLightbox (PARA LA GALERÍA)
    tags$link(rel = "stylesheet", href = "https://cdn.jsdelivr.net/npm/glightbox/dist/css/glightbox.min.css"),
    
    # --- ESTILOS CSS ---
    tags$style(HTML("
      :root {
        --primary-color: #34495e;
        --primary-hover: #4a657a; 
        --footer-bg: #2c3e50; 
        --footer-text: #ecf0f1; 
        --bg-white: #ffffff;
        --bg-light: #f8f9fa;
        --text-dark: #333;
        --text-light: #777;
        --shadow-sm: 0 4px 12px rgba(0,0,0,0.05);
        --shadow-md: 0 6px 20px rgba(0,0,0,0.1);
        --icon-size: 100px;
      }

      html, body {
        font-family: 'Lato', sans-serif;
        font-size: 16px;
        scroll-behavior: smooth;
        margin: 0;
        padding: 0;
        height: 100%;
        overflow-x: hidden;
      }
      
      body {
        display: flex;
        flex-direction: column;
        min-height: 100vh;
      }
      .shiny-html-output { 
         /* This targets the main div(style='display: flex...') */
         flex-grow: 1;
         display: flex; /* Make sure it's a flex container */
         flex-direction: column; /* Aligned to match its internal style */
      }

      /* 1. Sección de Bienvenida */
      .welcome-section {
        padding: 3rem 2rem; 
        text-align: center;
        background-color: var(--bg-light);
        border-bottom: 1px solid #dee2e6;
      }
      .welcome-section h1 {
        color: var(--primary-color);
        font-weight: 700;
        font-size: 2.5rem; 
        margin-bottom: 0.5rem;
      }
      .welcome-section p.subtitle {
        font-size: 1.15rem; 
        color: var(--text-light);
        font-weight: 300;
        max-width: 700px;
        margin: 0 auto;
      }

      /* 2. Sección de Contenido Principal (Descripción + Herramientas) */
      .main-content-section {
        padding: 2.5rem 1rem; 
        background-color: var(--bg-white);
      }
      
      /* Tarjeta de Descripción */
      .description-card {
        border: 1px solid #dee2e6;
        box-shadow: none;
        border-radius: 8px;
        display: flex;
        flex-direction: column;
        height: 100%; 
      }
      .description-card .card-header {
        background-color: var(--primary-color);
        color: white;
        font-weight: 600;
        font-size: 1.1rem;
        border-top-left-radius: 8px;
        border-top-right-radius: 8px;
        text-align: center; 
      }
      .description-card .card-body {
        font-size: 0.95rem;
        line-height: 1.7;
        text-align: justify; 
        flex-grow: 1; 
        padding: 1.5rem; 
      }
      
      /* Contenedor para las tarjetas de herramientas apiladas */
      .stacked-tools-container {
        display: flex;
        flex-direction: column;
        height: 100%; 
        justify-content: center; 
        gap: 1.5rem; 
      }
      
      /* Tarjeta de Herramienta (Layout: Icon - Name - Button) */
      .tool-card-stacked {
        min-height: calc(var(--icon-size) + 2rem); 
        border: 0;
        box-shadow: var(--shadow-md); 
        border-radius: 10px;
        transition: all 0.3s ease;
        background-color: var(--bg-white);
        overflow: hidden; 

        display: flex;
        flex-direction: row; 
        align-items: center; 
        padding: 1rem;
        gap: 1.5rem;
      }
      .tool-card-stacked:hover {
        transform: scale(1.02);
        box-shadow: 0 8px 24px rgba(0,0,0,0.15);
      }
      
      /* Contenedor del icono redondo */
      .tool-card-stacked .card-img-container {
        flex: 0 0 var(--icon-size); 
        height: var(--icon-size);   
        border-radius: 10%;      
        overflow: hidden;         
        background-color: transparent; 
        margin: 0; 
      }
      
      .tool-card-stacked .card-img-container img {
        width: 100%;
        height: 100%;
        object-fit: contain; 
      }
      
      /* Cuerpo de la tarjeta (Nombre + Botón) */
      .tool-card-stacked .card-body {
        flex-grow: 1;   
        flex-shrink: 1; 
        min-width: 0;   
        padding: 0;     
        text-align: left;
      }
      .tool-card-stacked .card-title {
        font-weight: 700;
        font-size: 1.4rem; 
        color: var(--primary-color);
        white-space: normal;
        overflow-wrap: break-word; 
        margin: 0; 
      }
      .btn-tool-launch {
        background-color: var(--primary-color);
        color: white;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        border-radius: 50px;
        padding: 6px 18px; 
        font-size: 0.85rem; 
        border: 0;
        transition: all 0.3s ease;
        white-space: nowrap; 
        flex-shrink: 0; 
      }
      .btn-tool-launch:hover {
        background-color: var(--primary-hover);
      }

      .clickable-tool-card {
        cursor: pointer;
      }
      
      /* 3. Sección de Footer (Centrado) */
      .footer-section {
        padding: 1.5rem 1rem; 
        background-color: var(--footer-bg);
        color: var(--footer-text);
        width: 100%; 
        box-sizing: border-box; 
        margin-top: 0px; 
        display: flex; 
        justify-content: center;
      }
      .footer-content {
        width: 100%;
        max-width: 1140px; 
        display: flex;
        flex-wrap: wrap; 
        justify-content: space-between; 
        align-items: center;
        gap: 1rem;
      }
      .footer-left {
        font-size: 0.85rem;
        font-weight: 300;
        line-height: 1.5;
        text-align: left;
        flex-basis: 60%; 
        flex-grow: 1;
      }
      .footer-right {
        display: flex;
        align-items: center;
        gap: 10px; 
        text-align: right;
        font-size: 0.85rem;
        font-weight: 300;
        flex-basis: 35%; 
        flex-grow: 1;
        justify-content: flex-end; 
      }
      .footer-right .logo-container {
         background-color: var(--bg-white); 
         padding: 5px 8px; 
         border-radius: 4px; 
         display: inline-flex; 
         # align-items: center;
      }
      .footer-right img {
        height: 35px; 
        display: block; 
      }
      
      .navbar-toggler-icon {
        filter: invert(1) brightness(2);
      }
      .navbar-toggler {
         border-color: rgba(255,255,255,0.2) !important;
      }

      .navbar-page {
        flex-grow: 1;
        display: flex;
        flex-direction: column;
        min-height: 0; 
      }
      
      .navbar-page .tab-content {
        flex-grow: 1;
        min-height: 0; 
      }

      .navbar-page .tab-content > .tab-pane {
        height: auto;
        min-height: initial;
      }
      
      .navbar-page .tab-content > .tab-pane:has(.mot-fill-panel) {
        min-height: 0;
        padding: 1rem 0; 
      }
      
      .mot-fill-panel {
        height: calc(100vh - 80px); 
        margin: 0 auto;  
        width: 98%;     
        border-radius: 8px;        
        display: flex;
        flex-direction: column;
        overflow: hidden; 
      }
      
      .mot-fill-panel > * {
        flex-grow: 1;
        min-height: 0;
        overflow-y: hidden; 
      }

      @media (max-width: 992px) {
        
        html, body {
          height: auto !important;
          overflow-x: hidden;
          overflow-y: auto !important;
        }

        .mot-fill-panel {
          height: auto !important;
          max-height: none !important;
          display: block !important;
          overflow: visible !important;
          margin-bottom: 2rem;
        }

        .mot-fill-panel > * {
          height: auto !important;
          overflow-y: visible !important;
        }

        .bslib-page-sidebar {
           height: auto !important;
           min-height: 100vh !important; 
        }
        
        .bslib-page-sidebar > .main {
           overflow: visible !important;
           height: auto !important;
        }
        
        .bslib-sidebar-layout > .sidebar-toggle {
           z-index: 1000; /* Por encima de otros elementos */
        }}

    ")), # Fin de tags$style
    
    # --- SCRIPTS JS (CON GLIGHTBOX) ---
    tags$script(src = "https://cdn.jsdelivr.net/npm/glightbox/dist/js/glightbox.min.js"),
    tags$script(HTML("
      // Función para navegar a paneles desde botones
      function navigateToPanel(panel) {
        Shiny.setInputValue('navigateToPanel', panel, {priority: 'event'});
      }
      
      // Inicializador de GLightbox (para la Galería)
      // Se asegura de que se reinicie si las pestañas de Shiny interfieren
      document.addEventListener('DOMContentLoaded', function() {
        const lightbox = GLightbox({ selector: '.glightbox' });
        
        // Añadir un listener para clics en la navbar que puedan cambiar la pestaña
        // y necesitar reinicializar GLightbox
        document.querySelector('.navbar').addEventListener('click', function() {
           // Pequeña espera para que Shiny actualice el DOM
           setTimeout(() => {
              lightbox.reload();
           }, 200);
        });

        // Observador por si el contenido de la galería se carga dinámicamente
        const observer = new MutationObserver(function(mutations) {
            mutations.forEach(function(mutation) {
                if (mutation.type === 'childList') {
                    lightbox.reload();
                }
            });
        });
        
        const galleryElement = document.getElementById('instructions_nav_gallery');
        if (galleryElement) {
            observer.observe(galleryElement, { childList: true, subtree: true });
        }
      });
    "))
  ), # Fin de tags$head
  
  # --- Resto de la UI ---
  div(style="display: flex; flex-direction: column; min-height: 100vh;", 
      page_navbar(
        title = "MetaOmixTools",
        id = "mot_navbar",
        position = "fixed-top",
        bg = "#34495e",          
        inverse = TRUE, 
        
        # --- Panel "Home" ---
        nav_panel(
          title = "Home", br(),
          
          # --- 1. Sección de Bienvenida ---
          div(
            class = "welcome-section",
            h1("Welcome to MetaOmixTools"),
            p(class = "subtitle",
              "A modern suite for integrative meta-omics analysis combining gene ranking and functional enrichment.")
          ), # Fin div (welcome-section)
          
          # --- 2. Sección de Contenido Principal ---
          div(
            class = "main-content-section container", 
            bslib::layout_columns(
              col_widths = c(7, 5), 
              
              # Columna Izquierda: Descripción
              bslib::card(
                class = "description-card",
                bslib::card_header("Suite Description"),
                bslib::card_body(
                  p(style = "text-align:justify;",
                    "MetaOmixTools is a bioinformatics suite designed for the integration and meta-analysis of omics data. The suite is composed of two main tools: ", tags$strong("MetaRank"), " and ", tags$strong("MetaEnrichGO"), ".",
                    br(), br(),
                    tags$strong("MetaRank"), " specialises in the meta-analysis of pre-ranked (classified) gene lists. It accepts multiple studies (weighted or unweighted) and uses methods such as RankProd or RRA to generate a consensus gene ranking. Optionally, it performs an enrichment analysis of the top genes from the final consensus.",
                    br(), br(),
                    tags$strong("MetaEnrichGO"), " performs a meta-analysis at the functional enrichment level. It takes one or multiple unranked gene lists, performs an overrepresentation analysis (ORA) on each of them, allowing the use of known databases such as GO, KEGG, Reactome, as well as custom annotations, and then combines the p-values (e.g., Fisher's method) to identify consensus functional terms enriched across all studies."
                  )
                )
              ),
              
              # Columna Derecha: Herramientas Apiladas
              div(class = "stacked-tools-container",
                
                # Tarjeta MetaRank
                tags$div(class = "clickable-tool-card",
                         onclick = "navigateToPanel('MetaRank')",

                  div(class = "tool-card-stacked", 
                    
                    # Hijo 1: Contenedor de Imagen
                    div(class = "card-img-container",
                        tags$img(src = "images/1.png", alt = "MetaRank Visual")),
                    
                    # Hijo 2: Contenedor de Texto
                    div(class = "card-body", 
                        h3(class = "card-title", "MetaRank")))),
                
                # Tarjeta MetaEnrichGO (Ahora toda la tarjeta es clicable)
                tags$div(class = "clickable-tool-card",
                         onclick = "navigateToPanel('MetaEnrichGO')",
                  
                  div(class = "tool-card-stacked",
                    
                    # Hijo 1: Contenedor de Imagen
                    div(class = "card-img-container",
                        tags$img(src = "images/2.png", alt = "MetaEnrichGO Visual")),
                    
                    # Hijo 2: Contenedor de Texto
                    div(class = "card-body", # Usamos un DIV simple
                        h3(class = "card-title", "MetaEnrichGO"))))
                
              ) # Fin stacked-tools-container
            ) # Fin layout_columns
          ), # Fin main-content-section
          
          tags$div(style = "height: -200px;"),
          
        ), # Fin nav_panel("Home")
        
        
        # --- Pestañas de las herramientas ---
        tabPanel("MetaRank", tags$div(style = "height: 60px;"),
                 div(class = "mot-fill-panel", 
                     mod_metaRank_ui("metarank1"))),
        
        tabPanel("MetaEnrichGO", tags$div(style = "height: 60px;"),
                 div(class = "mot-fill-panel",
                     mod_metaEnrichGO_ui("mega1"))),
        
        # --- Menú "About" ---
        navbarMenu("About", 
                   tabPanel("Instructions", tags$div(style = "height: 25px;"),
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
                   
                   tabPanel("Contact", tags$div(style = "height: 25px;"),
                            navlistPanel(
                              id    = ("instructions_nav_contact"),
                              widths = c(3, 9), 
                              "Help Section",
                              br(),
                              tabPanel(title = "Contact us", value="Contact us",
                                       tags$iframe(
                                         src    = "contact.html",
                                         width  = "100%", 
                                         height = "900px",
                                         style  = "border:none;")),
                              
                              tabPanel(title = "Questions & Answers", value="Q&A",
                                       tags$iframe(
                                         src    = "FAQ_readme.html",
                                         width  = "100%", 
                                         height = "900px",
                                         style  = "border:none;")))),
                   
                   # --- PESTAÑA GALLERY (RESTAURADA) ---
                   tabPanel("Gallery", tags$div(style = "height: 25px;"),
                            navlistPanel(
                              id = "instructions_nav_gallery",
                              widths = c(3, 9),
                              "Gallery Sections",
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
                                                  tags$p("This section indicates the year of creation and the responsible entity, along with an image that redirects to the official website of that organisation.", style = "text-align: center; font-style: italic; margin-top: 5px;")))),
                              
                              tabPanel(title = "MetaRank",
                                       fluidRow(
                                         # --- SECCIÓN 1: OVERVIEW & INTERFACE (De README) ---
                                         column(6, tags$a(
                                           href = "images/metarank_overview.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_overview.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("General interface of the MetaRank module for gene ranking meta-analysis.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_rankprod_input.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_rankprod_input.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Configuration of the Weighted (RankProd) input method with file upload.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_rra_input.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_rra_input.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Configuration of the Unweighted (RRA) input method for ordered gene lists.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_rp_info.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_rp_info.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Information modal detailing the RankProd score-based data format.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_rra_info.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_rra_info.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Information modal detailing the required plain text format for RRA.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_origin_info.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_origin_info.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Detailed guidance on defining the 'Origin' vector for batch correction in Advanced RankProd.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         # --- SECCIÓN 2: PARAMETERS & FILTERS (De README) ---
                                         column(6, tags$a(
                                           href = "images/metarank_rankprod_parameters.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_rankprod_parameters.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Analysis settings for RankProd: ranking direction, NA management, and appearance filter.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_rra_parameters.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_rra_parameters.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("RRA specific parameters: choice of aggregation model and list frequency threshold.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_enrich.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_enrich.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Functional enrichment module options for post-ranking interpretation.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         # --- SECCIÓN 3: DATA VISUALIZATION TAB SETTINGS (De README) ---
                                         column(6, tags$a(
                                           href = "images/metarank_data_1.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_data_1.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Table configuration: selection of ranking metrics and metadata columns.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_data_2.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_data_2.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("UpSet Plot customization: text size and bar color selection.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_data_3.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_data_3.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Heatmap visualization settings: axis labels and color scale selection.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_data_4.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_data_4.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Enrichment table settings: toggling functional annotation columns.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_data_5.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_data_5.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Enrichment plot customization: type, colors, and number of top terms.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         # --- SECCIÓN 4: TROUBLESHOOTING & ERRORS (De README) ---
                                         column(6, tags$a(
                                           href = "images/metarank_error_1_no_input.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_1_no_input.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Validation error: Notification for missing input data lists.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_2_single_list.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_2_single_list.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Input error: Alert for providing fewer than two lists for meta-analysis.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_3_file_rra.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_3_file_rra.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Format error: Mismatch in RRA input file structure (headers or commas detected).", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_4_paste_rra.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_4_paste_rra.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Pasting error: Invalid formatting detected in RRA text-area input.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_5_file_rp.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_5_file_rp.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Format error: Mismatch in RankProd file (numeric scores missing).", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_6_paste_rp.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_6_paste_rp.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Pasting error: Invalid numeric data or alignment in RankProd input.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_7_origin.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_7_origin.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Validation error: Origin vector length does not match input file count.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_7_2_origin.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_7_2_origin.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Batch error: Insufficient replicates within the defined Origin groups.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_8_threshold.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_8_threshold.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Filter error: Appearance threshold set higher than any actual gene recurrence.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_9_no_geness.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_9_no_geness.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Zero overlap error: No genes shared among input datasets.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_9_no_genes.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_9_no_genes.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Enrichment failure: The analysis returned no significant biological terms.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_10_no_annotation.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_10_no_annotation.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Validation error: Personalized database selected without uploading an annotation file.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_11_annotation_format.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_11_annotation_format.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Format error: Custom annotation file does not follow the required two-column structure.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_12_geneid_format.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_12_geneid_format.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Blocking error: Most identifiers do not match the selected Gene ID type.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_12_2_geneid_format.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_12_2_geneid_format.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Non-blocking warning: Some invalid identifiers detected and will be excluded.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_13_organism_format.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_13_organism_format.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Blocking error: Mismatch between gene nomenclature and selected organism.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_error_13_2_organism_format.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_error_13_2_organism_format.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Non-blocking warning: Some gene symbols do not match the selected organism's casing.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_workflow.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_workflow.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Complete MetaRank workflow: from data ingestion to final visualization.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         # --- SECCIÓN 5: RESULTS & EXAMPLES (De EXAMPLES) ---
                                         column(6, tags$a(
                                           href = "images/example1_metarank_excluded.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example1_metarank_excluded.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Table of genes excluded due to recurrence filters in Example 1.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/example1_metarank_upset.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example1_metarank_upset.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Intersection analysis UpSet plot for Study 1 datasets.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/example1_metarank_heatmap.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example1_metarank_heatmap.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Pairwise similarity heatmap showing gene overlaps in Example 1.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/example2_metarank_upset.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example2_metarank_upset.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("UpSet plot for Advanced RankProd with yellow/purple custom styling.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/example2_metarank_heatmap.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example2_metarank_heatmap.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Portland-scale similarity heatmap for the Advanced RankProd analysis.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/example2_metarank_enrich.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example2_metarank_enrich.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Barplot of enriched GO biological processes from consensus genes.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/example3_metarank_excluded.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example3_metarank_excluded.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Table of genes excluded during RRA meta-analysis.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/example3_metarank_upset.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example3_metarank_upset.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("UpSet intersection analysis for the unweighted (RRA) datasets.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/example3_metarank_heatmap.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example3_metarank_heatmap.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Heatmap showing the proportion of shared genes in Example 3.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/example3_metarank_enrich.png", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/example3_metarank_enrich.png", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Personalized database enrichment barplot: protein location results.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, tags$a(
                                           href = "images/metarank_upsetplot.jpg", class = "glightbox", `data-gallery` = "gallery-metarank",
                                           tags$img(src = "images/metarank_upsetplot.jpg", width = "100%", style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                           tags$p("Default UpSet plot styling showing intersections across studies.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, 
                                                tags$div(style = "border: 1px solid #ccc; border-radius: 8px; overflow: hidden; height: 400px; background: white;",
                                                         tags$iframe(src = "images/metarank_heatmap.html", width = "100%", height = "100%", frameborder = "0")),
                                                tags$p("Interactive similarity heatmap: explore pairwise gene list overlaps.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                         
                                         column(6, 
                                                tags$div(style = "border: 1px solid #ccc; border-radius: 8px; overflow: hidden; height: 400px; background: white;",
                                                         tags$iframe(src = "images/metarank_bar_html.html", width = "100%", height = "100%", frameborder = "0")),
                                                tags$p("Interactive Enrichment Bar Plot: dynamic tooltips and filtering.", style = "text-align: center; font-style: italic; margin-top: 5px;"))
                                       )),
                              
                              # --- CONTENIDO DE METAENRICHGO RESTAURADO ---
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
                                    tags$p("General overview of the MetaEnrichGO module interface.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_input.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_input.png", 
                                             width = "100%", 
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Selection of input method: file upload (.txt/.tsv) or manual gene list pasting.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_i.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_i.png",
                                             width = "100%", 
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Interactive help tooltips providing details on format requirements and parameters.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_parameters.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_parameters.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Configuration of analysis parameters: database selection, organism, and meta-analysis method.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metarank_personalized.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metarank_personalized.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Configuration panel for personalized annotations using custom user-provided files.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_resultstable.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_resultstable.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Interactive results table showing combined p-values and term descriptions.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_excludedtable.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_excludedtable.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Table of terms excluded due to the minimum dataset frequency threshold.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_dot.jpg",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_dot.jpg",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Interactive Dotplot visualization of enriched biological processes.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_bar.jpg",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_bar.jpg",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Barplot representation of the top significant enriched terms.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_workflow.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_workflow.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Functional workflow: from gene list input to over-representation and meta-analysis.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_1_no_input_data.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_1_no_input_data.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Validation error: No input data detected before running analysis.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_2_no_annotation_data.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_2_no_annotation_data.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Validation error: Missing custom annotation file for Personalized database mode.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_3_file_format.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_3_file_format.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Validation error: Incorrect file format or forbidden characters detected in upload.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_4_paste_format.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_4_paste_format.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Validation error: Formatting issues in manually pasted gene lists (e.g., hidden tabs).", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_5_annotation_format.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_5_annotation_format.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Validation error: Annotation file does not match the required two-column structure.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_6_organism_format.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_6_organism_format.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Nomenclature mismatch: Gene symbols do not match the selected organism's case format.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_6_2_organism_format.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_6_2_organism_format.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Warning: High percentage of genes excluded due to organism nomenclature mismatch.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_7_geneid_format.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_7_geneid_format.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Validation error: Input identifiers do not match the selected ID type (SYMBOL/ENTREZ/ENSEMBL).", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_7_2_geneid_format.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_7_2_geneid_format.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Warning: Some identifiers were excluded as they didn't follow the expected ID pattern.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_8_no_mapping.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_8_no_mapping.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Critical error: No identifiers were successfully mapped to the selected database.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_error_9_to_strict.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_error_9_to_strict.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Filter error: No terms reached the minimum dataset frequency threshold.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_example1dot.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_example1dot.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Dot plot visualization for Basic Example: displaying the top 20 enriched KEGG pathways.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/combined.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/combined.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Combined meta-analysis results: consensus enrichment terms across all input datasets.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/genes1.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/genes1.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Individual enrichment analysis results for the first dataset (genes1.txt).", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/genes2.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/genes2.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Individual enrichment analysis results for the second dataset (genes2.txt).", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/genes3.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/genes3.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Individual enrichment analysis results for the third dataset (genes3.txt).", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/genes4.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/genes4.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Individual enrichment analysis results for the fourth dataset (genes4.txt).", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, tags$a(
                                    href = "images/metaenrichgo_example.png",
                                    class = "glightbox",
                                    `data-gallery` = "gallery-metaenrichgo",
                                    tags$img(src = "images/metaenrichgo_example.png",
                                             width = "100%",
                                             style = "border: 1px solid #ccc; border-radius: 8px; height: 400px; object-fit: contain; background: white;")),
                                    tags$p("Table of excluded terms: biological processes that did not meet the minimum recurrence threshold.", style = "text-align: center; font-style: italic; margin-top: 5px;")),
                                  
                                  column(6, 
                                         tags$div(
                                           style = "border: 1px solid #ccc; border-radius: 8px; overflow: hidden; height: 400px; background: white;",
                                           tags$iframe(
                                             src = "images/MetaEnrichGO_plot_2025-05-08.html",
                                             width = "100%",
                                             height = "100%",
                                             frameborder = "0"
                                           )
                                         ),
                                         tags$p("Interactive HTML Dot Plot: dynamic visualization with tooltips and zoom capabilities.", style = "text-align: center; font-style: italic; margin-top: 5px;")
                                  ))
                              ))),
                   
                   
                   tabPanel("GitHub", tags$div(style = "height: 25px;"),
                            tags$div(style = "padding: 25px; font-family: 'Segoe UI', sans-serif; line-height: 1.7; max-width: 900px; margin: auto;",
                                     
                                     tags$h3("MetaOmixTools Suite & Individual Tools", style = "margin-bottom: 20px;"),
                                     
                                     tags$h4("Online Tools", style = "margin-top: 30px;"),
                                     tags$p("You can launch each application individually:"),
                                     tags$ul(
                                       tags$li(
                                         tags$b("MetaOmixTools Suite: "),
                                         tags$a(href = "https://maksymkupchyktiurin00.shinyapps.io/MetaOmixTools/", 
                                                "Launch App", target = "_blank", style = "color: #007ACC; text-decoration: none;")),
                                       tags$li(
                                         tags$b("MetaRank: "),
                                         tags$a(href = "https://maksymkupchyktiurin00.shinyapps.io/MetaRank/", 
                                                "Launch App", target = "_blank", style = "color: #007ACC; text-decoration: none;")),
                                       tags$li(
                                         tags$b("MetaEnrichGO: "),
                                         tags$a(href = "https://maksymkupchyktiurin00.shinyapps.io/MetaEnrichGO/", 
                                                "Launch App", target = "_blank", style = "color: #007ACC; text-decoration: none;"))),
                                     
                                     
                                     tags$h4("GitHub Repositories", style = "margin-top: 30px;"),
                                     tags$p("All apps are open-source and available on GitHub:"),
                                     tags$ul(
                                       tags$li(tags$b("MetaOmixTools Suite: "), 
                                               tags$a(href = "https://github.com/EMEKATE00/MetaOmixTools.git", 
                                                      "GitHub Repo", target = "_blank")),
                                       tags$li(tags$b("MetaRank: "), 
                                               tags$a(href = "https://github.com/EMEKATE00/MetaRank.git", 
                                                      "GitHub Repo", target = "_blank")),
                                       tags$li(tags$b("MetaEnrichGO: "), 
                                               tags$a(href = "https://github.com/EMEKATE00/MetaEnrichGO.git", 
                                                      "GitHub Repo", target = "_blank"))),
                                     
                                     
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
                                       tags$li(tags$b("Documentation:"), " Every module and function is documented for both users and developers.")),
                                     
                                     tags$div(style = "margin-top: 40px;",
                                              tags$a(href = "https://github.com/EMEKATE00/MetaOmixTools.git", target = "_blank",
                                                     class = "btn btn-primary", 
                                                     style = "color: white; background-color: #24292e; border: none; padding: 10px 18px; font-size: 16px; border-radius: 6px;",
                                                     shiny::icon("github"), " View Suite on GitHub")
                                     )))
                   ) # Fin tabPanel Gallery
      ) # Fin NavbarMenu About
  ), # Fin de page_navbar
  
  
  # --- FOOTER ---
  div(class = "footer-section",
    div(class = "footer-content", 
      div(class = "footer-left",
          HTML("&copy; Copyright ") , paste(format(Sys.Date(), "%Y")), HTML(", MetaOmixTools."), br(), 
          "A modern suite for integrative meta-omics analysis combining gene ranking and functional enrichment."),
      div(class = "footer-right",
          span("Developed by The Centro de Investigación Príncipe Felipe (CIPF)"),
          div(class="logo-container", 
             tags$a(href = "https://www.cipf.es/", target = "_blank",
                    tags$img(src = "images/5.png", alt = "CIPF Logo"))))) 
  
    ) # Fin footer-section
) # Fin div contenedor principal flex



#------------------------------------------------------------------------------
# Define the server logic
#------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  mod_metaRank_server("metarank1")
  mod_metaEnrichGO_server("mega1")
  
  observeEvent(input$navigateToPanel, {
    req(input$navigateToPanel)
    updateNavbarPage(session, "mot_navbar", selected = input$navigateToPanel)
  })
  
}

#------------------------------------------------------------------------------
# Launch the Shiny application
#------------------------------------------------------------------------------

shinyApp(ui = ui, server = server)