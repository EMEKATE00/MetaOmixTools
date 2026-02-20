# MetaRank_Functions.R
library(dplyr)
library(ggplot2)
library(RankProd)
library(bslib)
library(DT)
library(shinycssloaders)
library(shinyWidgets)
library(shiny)
library(shinyBS)
library(plotly)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)
library(UpSetR)
library(ReactomePA)
library(reactome.db)
library(RobustRankAggreg)
library(viridis)
library(RColorBrewer)
library(magick)

###############################################################################

#-------------------------------- VALIDATIONS --------------------------------#

###############################################################################

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
#' @2.MIN_LISTS_CHECK
# ------------------------------------------------------------------------------
validate_min_lists_metarank <- function(data_lists) {
  if (length(data_lists) < 2) {
    
    showModal(modalDialog(
      title = tags$div(icon("exclamation-circle"), " Insufficient Data", 
                       style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "m",
      
      tags$div(
        class = "alert alert-danger", role = "alert",
        tags$h4(icon("layer-group"), " Single List Detected"),
        tags$p("The meta-analysis cannot be performed because only one gene list has been provided.")),
      
      tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #c0392b; margin-top: 15px;",
        tags$h5("Action Required:", style="font-weight:bold; margin-top:0;"),
        tags$ul(
          tags$li("Upload at least ", tags$b("two files"), "."),
          tags$li("Or paste multiple lists separated by ", tags$code("###"), "."))),
      
      footer = modalButton("Close")))
    return(FALSE)}
  return(TRUE)}



# ------------------------------------------------------------------------------
#' @3.FILE_FORMAT_CHECK_RRA
# ------------------------------------------------------------------------------
validate_file_format_rra <- function(input, session) {
  bad_files <- list()
  files <- if (input$inputMethod == "files") input$files_rra$datapath else NULL
  file_names <- if (input$inputMethod == "files") input$files_rra$name else NULL
  
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
      title = tags$div(icon("exclamation-circle"), " Format Error: RRA File Input", 
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
        tags$h5("Expected Format (RRA):", style="font-weight:bold; margin-top:0;"),
        tags$ul(tags$li("Plain text file (.txt) with one gene per line."),
                tags$li("No headers, comments (#), commas, or tabs."))),
      
      footer = modalButton("Close and Fix Files")))
    return(FALSE)}
  return(TRUE)}



# ------------------------------------------------------------------------------
#' @4.PASTE_FORMAT_CHECK_RRA
# ------------------------------------------------------------------------------
validate_paste_format_rra <- function(input, session) {
  if (!(input$inputMethod == "paste" && nzchar(input$pasteText_rra))) return(TRUE)
  
  txt <- input$pasteText_rra
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
      title = tags$div(icon("exclamation-circle"), " Format Error: RRA Paste Input", 
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
        tags$h5("Expected Format (RRA)", style="margin-top:0; font-weight:bold;"),
        tags$ul(tags$li("Single column of gene symbols."),
                tags$li("Separate multiple lists using ", tags$code("###"), "."),
                tags$li("No headers, hidden tabs, or extra columns."),
                tags$li("Tip: Paste into Notepad first to clean hidden formatting."))),
      
      footer = modalButton("Close and Fix")))
    
    shinyjs::disable("Run"); shinyjs::delay(1000, shinyjs::enable("Run"))
    return(FALSE)}
  return(TRUE)}



# ------------------------------------------------------------------------------
#' @5.FILE_FORMAT_CHECK_RP 
# ------------------------------------------------------------------------------
validate_file_format_rp <- function(input, session) {
  files <- if (input$inputMethod == "files") input$files$datapath else NULL
  file_names <- if (input$inputMethod == "files") input$files$name else NULL
  
  if (is.null(files)) return(TRUE)
  
  fatal_errors <- list()   # Errores que detienen el análisis (Rojo)
  warning_msgs <- list()   # Advertencias que permiten continuar (Amarillo)
  css_box <- "margin-top:4px; background-color:#fff3cd; padding:5px; border-radius:3px; font-family:monospace; font-size:0.85em; color:#856404;"
  
  for (i in seq_along(files)) {
    fname <- file_names[i]
    fpath <- files[i]
    
    lines_raw <- readLines(fpath, warn = FALSE, encoding = "UTF-8")
    
    bad_hash <- grep("^#", lines_raw)
    if (length(bad_hash) > 0) {
      preview <- htmltools::htmlEscape(lines_raw[bad_hash[1]])
      fatal_errors[[fname]] <- paste0(
        "<div style='margin-bottom: 10px; border-bottom: 1px solid #eee; padding-bottom: 5px;'>",
        "<strong>Header/Comment (#):</strong> Detected at line ", bad_hash[1], 
        "<div style='", css_box, "'>", preview, "</div>",
        "</div>")
      next 
    }
    
    ext <- tolower(tools::file_ext(fname))
    sep_char <- if (ext == "csv") "," else "\t"
    
    df <- tryCatch({
      read.table(fpath, header = FALSE, sep = sep_char, 
                 stringsAsFactors = FALSE, fill = TRUE, 
                 quote = "", comment.char = "")
    }, error = function(e) NULL)
    
    if (is.null(df) || nrow(df) == 0) {
      fatal_errors[[fname]] <- paste0(
        "<div style='margin-bottom: 10px; border-bottom: 1px solid #eee; padding-bottom: 5px;'>",
        "<strong>Read Error:</strong> File appears empty or unreadable.",
        "<div style='", css_box, "'>Could not parse table structure.</div></div>")
      next
    }
    
    if (ncol(df) != 2) {
      fatal_errors[[fname]] <- paste0(
        "<div style='margin-bottom: 10px; border-bottom: 1px solid #eee; padding-bottom: 5px;'>",
        "<strong>Column Count Error:</strong> Detected <b>", ncol(df), 
        "</b> columns. Required <b>2</b>.",
        "<div style='", css_box, "'>Check separator (", ifelse(sep_char==",","Comma","Tab"), 
        ") or extra delimiters.</div></div>")
      next
    }
    
    vals_raw <- as.character(df[[2]])
    vals_raw <- gsub(",", ".", vals_raw) 
    vals_num <- suppressWarnings(as.numeric(vals_raw))
    genes_raw <- trimws(as.character(df[[1]]))
    
    is_bad_row <- (genes_raw == "") | is.na(vals_num)
    bad_row_count <- sum(is_bad_row)
    total_rows <- nrow(df)
    
    if (bad_row_count > 0) {
      if (bad_row_count == total_rows) {
        fatal_errors[[fname]] <- paste0(
          "<div style='margin-bottom: 10px; border-bottom: 1px solid #eee; padding-bottom: 5px;'>",
          "<strong>No Valid Data:</strong> All ", total_rows, " rows are invalid.",
          "<div style='", css_box, "'>Verify decimal format (use dot or comma) and ensure no missing genes.</div></div>")
      } else {
        preview_bad_idx <- which(is_bad_row)[1]
        bad_gene <- genes_raw[preview_bad_idx]; if(is.na(bad_gene)) bad_gene <- ""
        bad_val  <- vals_raw[preview_bad_idx];  if(is.na(bad_val)) bad_val <- ""
        preview_str <- paste(bad_gene, bad_val, sep = ifelse(sep_char==",",",","\t"))
        
        warning_msgs[[fname]] <- paste0(
          "<div style='margin-bottom: 10px; border-bottom: 1px solid #eee; padding-bottom: 5px;'>",
          "<strong>Incomplete Rows Ignored:</strong> <b>", bad_row_count, 
          "</b> rows will be skipped.",
          "<div style='", css_box, "'>First instance at row ", preview_bad_idx, ": ", htmltools::htmlEscape(preview_str), "</div></div>")}}}
  
  if (length(fatal_errors) > 0) {
    showModal(modalDialog(
      title = tags$div(icon("times-circle"), " Analysis Stopped: RP File Format Error", 
                       style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      size = "l", easyClose = TRUE,
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("ban"), " Critical Format Issues"),
               tags$p("The analysis cannot proceed due to the following errors:")),
      tags$div(style="max-height: 350px; overflow-y: auto; border: 1px solid #e2e2e2; border-radius: 5px; padding: 10px;",
               lapply(names(fatal_errors), function(f) {
                 tags$div(style = "margin-bottom: 15px;",
                          tags$h5(icon("file-alt"), tags$b(f), style="color: #c0392b; margin-top:0; border-bottom: 2px solid #c0392b; padding-bottom:5px;"),
                          HTML(fatal_errors[[f]]))})),
      footer = modalButton("Close and Fix")))
    return(FALSE) # Detiene el proceso
    }
  
  if (length(warning_msgs) > 0) {
    showModal(modalDialog(
      title = tags$div(icon("exclamation-triangle"), " Data Warning", 
                       style = "color: #e67e22; font-weight: bold; font-size: 18px;"),
      size = "l", easyClose = TRUE,
      tags$div(class = "alert alert-warning", role = "alert",
               tags$h4(icon("info-circle"), " Analysis Continuing with Warnings"),
               tags$p("Some rows were incomplete and have been automatically ignored. The analysis will proceed with the valid data.")),
      tags$div(style="max-height: 350px; overflow-y: auto; border: 1px solid #e2e2e2; border-radius: 5px; padding: 10px;",
               lapply(names(warning_msgs), function(f) {
                 tags$div(style = "margin-bottom: 15px;",
                          tags$h5(icon("file-alt"), tags$b(f), style="color: #e67e22; margin-top:0; border-bottom: 2px solid #e67e22; padding-bottom:5px;"),
                          HTML(warning_msgs[[f]]))})),
      footer = modalButton("Ok, I understand")))
    return(TRUE) # Permite continuar
  }
  
  return(TRUE)}



# ------------------------------------------------------------------------------
#' @6.PASTE_FORMAT_CHECK_RP
# ------------------------------------------------------------------------------
validate_paste_format_rp <- function(input, session) {
  if (input$inputMethod != "paste") return(TRUE)
  
  if (input$pasteFormat == "csv") {
    txt <- input$pasteText; sep_char <- ","; req_txt <- input$pasteText
  } else {
    txt <- input$pasteText1; sep_char <- "\t"; req_txt <- input$pasteText1}
  
  if (!nzchar(req_txt)) return(TRUE)
  
  blocks <- strsplit(txt, "###", fixed = TRUE)[[1]]
  
  fatal_errors <- list()
  warning_msgs <- list()
  
  css_box <- "margin-top:2px; background-color:#fff3cd; padding:4px; border-radius:3px; font-family:monospace; font-size:0.85em; color:#856404;"
  
  for (i in seq_along(blocks)) {
    b <- trimws(blocks[i])
    list_name <- paste0("List ", i)
    
    if (!nzchar(b)) next
    
    lines <- unlist(strsplit(b, "\n"))
    list_errors_html <- character() 
    
    # -------------------------------------------------------
    # 1. CHEQUEO DE CARACTERES PROHIBIDOS (Punto y coma ;)
    # -------------------------------------------------------
    bad_semi_idx <- grep(";", lines, fixed = TRUE)
    if (length(bad_semi_idx) > 0) {
      preview_line <- lines[bad_semi_idx[1]]
      preview_html <- htmltools::htmlEscape(substr(preview_line, 1, 60))
      
      msg <- paste0(
        "<li><strong>Invalid Character (;):</strong> Detected semicolon at line ", bad_semi_idx[1], ".",
        "<div style='", css_box, "'>", preview_html, "</div></li>"
      )
      list_errors_html <- c(list_errors_html, msg)}
    
    # -------------------------------------------------------
    # 2. CHEQUEO DE CABECERAS (#)
    # -------------------------------------------------------
    bad_hash_idx <- grep("^#", lines)
    if (length(bad_hash_idx) > 0) {
      preview_line <- lines[bad_hash_idx[1]]
      preview_html <- htmltools::htmlEscape(substr(preview_line, 1, 60))
      
      msg <- paste0(
        "<li><strong>Header/Comment (#):</strong> Detected at line ", bad_hash_idx[1], ".",
        "<div style='", css_box, "'>", preview_html, "</div></li>"
      )
      list_errors_html <- c(list_errors_html, msg)}
    
    # -------------------------------------------------------
    # 3. CHEQUEO DE ESTRUCTURA Y COLUMNAS
    # -------------------------------------------------------
    df <- tryCatch({
      read.table(text = b, header = FALSE, sep = sep_char, 
                 stringsAsFactors = FALSE, fill = TRUE, quote = "")
    }, error = function(e) NULL)
    
    if (is.null(df)) {
      # Error crítico de parseo
      msg <- paste0(
        "<li><strong>Parse Error:</strong> Could not interpret text structure.",
        "<div style='", css_box, "'>Ensure you are using the correct separator (", 
        ifelse(sep_char == ",", "Comma", "Tab"), ").</div></li>")
      list_errors_html <- c(list_errors_html, msg)
      
    } else {
      if (ncol(df) != 2) {
        msg <- paste0(
          "<li><strong>Column Mismatch:</strong> Found <b>", ncol(df), "</b> columns. Required <b>2</b>.",
          "<div style='", css_box, "'>Check if your data has extra delimiters or missing values.</div></li>"
        )
        list_errors_html <- c(list_errors_html, msg)}
      
      # -------------------------------------------------------
      # 4. CHEQUEO DE FILAS INCOMPLETAS (Solo si estructura es válida)
      # -------------------------------------------------------
      if (ncol(df) == 2 && length(list_errors_html) == 0) {
        vals_raw <- as.character(df[[2]])
        vals_raw <- gsub(",", ".", vals_raw)
        vals_num <- suppressWarnings(as.numeric(vals_raw))
        genes_raw <- trimws(as.character(df[[1]]))
        
        is_bad_row <- (genes_raw == "") | is.na(vals_num)
        bad_count <- sum(is_bad_row)
        
        if (bad_count > 0) {
          if (bad_count == nrow(df)) {
            msg <- paste0(
              "<li><strong>No Valid Data:</strong> All rows are invalid (non-numeric values or empty genes).</li>"
            )
            list_errors_html <- c(list_errors_html, msg)
          } else {
            preview_idx <- which(is_bad_row)[1]
            bad_gene <- genes_raw[preview_idx]; if(is.na(bad_gene)) bad_gene <- ""
            bad_val  <- vals_raw[preview_idx];  if(is.na(bad_val)) bad_val <- ""
            preview_str <- paste(bad_gene, bad_val, sep=ifelse(sep_char==",", ",", "\t"))
            
            warning_msgs[[list_name]] <- paste0(
              "<div><strong>Incomplete Rows:</strong> <b>", bad_count, "</b> rows will be skipped.",
              "<div style='", css_box, "'>Example (row ", preview_idx, "): ", 
              htmltools::htmlEscape(preview_str), "</div></div>")}}}}
    
    if (length(list_errors_html) > 0) {
      fatal_errors[[list_name]] <- paste0("<ul style='padding-left: 20px; margin-bottom: 0;'>", 
                                          paste(list_errors_html, collapse=""), "</ul>")}}
  
  # --- MOSTRAR ERRORES FATALES (ROJO) ---
  if (length(fatal_errors) > 0) {
    showModal(modalDialog(
      title = tags$div(icon("times-circle"), " PR Paste Format Error", 
                       style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "l",
      
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("ban"), " Analysis Stopped"),
               tags$p("Please fix the following issues in your pasted text:")),
      
      tags$div(style="max-height: 400px; overflow-y: auto; border: 1px solid #e2e2e2; border-radius: 5px; padding: 15px;",
               lapply(names(fatal_errors), function(f) {
                 tags$div(style = "margin-bottom: 15px; border-bottom: 1px solid #eee; padding-bottom: 10px;",
                          tags$h5(icon("list-ul"), tags$b(f), style="color: #c0392b; margin-top:0;"),
                          HTML(fatal_errors[[f]]) )})),
      
      footer = modalButton("Close and Fix")))
    
    shinyjs::disable("Run"); shinyjs::delay(1000, shinyjs::enable("Run"))
    return(FALSE)}
  
  # --- MOSTRAR ADVERTENCIAS (AMARILLO) ---
  if (length(warning_msgs) > 0) {
    showModal(modalDialog(
      title = tags$div(icon("exclamation-triangle"), " Data Warning", 
                       style = "color: #e67e22; font-weight: bold; font-size: 18px;"),
      size = "l", easyClose = TRUE,
      tags$div(class = "alert alert-warning", role = "alert",
               tags$h4(icon("info-circle"), " Analysis Continuing with Warnings"),
               tags$p("Some rows contained errors and will be automatically skipped.")),
      
      tags$div(style="max-height: 350px; overflow-y: auto; border: 1px solid #e2e2e2; border-radius: 5px; padding: 10px;",
               lapply(names(warning_msgs), function(f) {
                 tags$div(style = "margin-bottom: 15px;",
                          tags$h5(tags$b(f), style="color: #e67e22; margin-top:0; border-bottom: 1px solid #ddd;"),
                          HTML(warning_msgs[[f]]) )})),
      
      footer = modalButton("Ok, I understand")))
    return(TRUE)}
  
  return(TRUE)}



# ------------------------------------------------------------------------------
#' @7.ORIGIN_VALIDATION 
# ------------------------------------------------------------------------------
validate_origin_rp <- function(input, session, pd) {
  if (!(input$selectedPackage == "Rankprodpack" && input$RankBasedMethod == "rp.adv")) {
    return(TRUE)}
  
  req(pd$valid)
  n_lists <- length(pd$valid)
  origin_raw <- input$origin
  
  css_box <- "margin-top:4px; background-color:#fff3cd; padding:5px; border-radius:3px; font-family:monospace; font-size:0.85em; color:#856404;"
  
  if (is.null(origin_raw) || !nzchar(origin_raw)) {
    showModal(modalDialog(
      title = tags$div(icon("times-circle"), " Origin Required", 
                       style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "m",
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("ban"), " Missing Parameter"),
               tags$p("You selected 'RankProd Advanced' but did not provide an Origin vector.")),
      tags$div(style = css_box, "Please enter a comma-separated numeric vector (e.g., 1,1,2,2)."),
      footer = modalButton("Close")))
    return(FALSE)}
  
  # 2. Caracteres inválidos (Permite solo números y comas)
  if (grepl("[^0-9,]", origin_raw)) {
    showModal(modalDialog(
      title = tags$div(icon("times-circle"), " Invalid Format", 
                       style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "m",
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("keyboard"), " Character Error"),
               tags$p("The origin vector must contain ONLY numbers and commas.")),
      tags$div(style = css_box, paste("Invalid Input:", htmltools::htmlEscape(origin_raw))),
      footer = modalButton("Close")))
    return(FALSE)}
  
  # Procesar vector
  origin_vals <- strsplit(origin_raw, ",")[[1]]
  origin_vals <- origin_vals[nzchar(origin_vals)] # Eliminar vacíos accidentales "1,,2"
  
  # 3. Longitud
  if (length(origin_vals) != n_lists) {
    showModal(modalDialog(
      title = tags$div(icon("times-circle"), " Length Mismatch", 
                       style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "m",
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("ruler-horizontal"), " Dimension Error"),
               tags$p(HTML(paste0("You have <b>", n_lists, "</b> valid gene lists, but provided <b>", length(origin_vals), "</b> origin values.")))),
      tags$div(style = css_box, "The number of origin values must exactly match the number of lists."),
      footer = modalButton("Close")))
    return(FALSE)}
  
  # 4. Replicados (Minimo 2 por grupo)
  counts <- table(origin_vals)
  singletons <- names(counts[counts < 2])
  
  if (length(singletons) > 0) {
    showModal(modalDialog(
      title = tags$div(icon("times-circle"), " Insufficient Replicates", 
                       style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "m",
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("layer-group"), " Grouping Error"),
               tags$p("Each origin group must contain at least 2 lists.")),
      tags$div(style = css_box, paste("Groups with only 1 list detected:", paste(singletons, collapse=", "))),
      footer = modalButton("Close")))
    return(FALSE)}
  
  return(TRUE)}



# ------------------------------------------------------------------------------
#' @8.THRESHOLD_VALIDATION 
# ------------------------------------------------------------------------------
validate_gene_threshold <- function(all_genes_list, file_names, input, session) {
  
  pkg <- input$selectedPackage
  slider_id <- if (pkg == "Rankprodpack") "gene_app_filter_rp" else "gene_app_filter_rra"
  user_threshold <- input[[slider_id]]
  
  app_df <- count_gene_appearance_rra(all_genes_list, file_names)
  max_freq_found <- if (nrow(app_df) > 0) max(app_df$Count) else 0
  
  css_box <- "margin-top:4px; background-color:#fff3cd; padding:5px; border-radius:3px; font-family:monospace; font-size:0.85em; color:#856404;"
  
  # --- CASO A: NO HAY GENES COMUNES (STOP) ---
  if (max_freq_found < 2) {
    showModal(modalDialog(
      title = tags$div(icon("times-circle"), " No Shared Genes Found", 
                       style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "m",
      
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("ban"), " Meta-Analysis Impossible"),
               tags$p("The analysis cannot proceed because there are no genes shared between any of the input lists.")),
      
      tags$div(style = "margin-bottom: 15px;",
               "<strong>Diagnostic:</strong>",
               tags$div(style = css_box, "Max overlap found: 1 (Unique lists). Required: >= 2.")),
      
      tags$p("Please check gene nomenclature (e.g., SYMBOL vs ENTREZ) or ensure the datasets are comparable."),
      
      footer = modalButton("Close")))
    
    return(NULL) # Devuelve NULL para detener
  }
  
  # --- CASO B: UMBRAL DEMASIADO ALTO (WARNING + AUTO-FIX) ---
  if (user_threshold > max_freq_found) {
    
    updateSliderInput(session, slider_id, value = max_freq_found)
    
    showModal(modalDialog(
      title = tags$div(icon("exclamation-triangle"), " Threshold Too Strict", 
                       style = "color: #e67e22; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "m",
      
      tags$div(class = "alert alert-warning", role = "alert",
               tags$h4(icon("wrench"), " Threshold Auto-Adjusted"),
               tags$p("No genes meet the selected appearance threshold. The analysis will continue with the maximum possible overlap.")),
      
      tags$div(style = "margin-bottom: 10px;",
               "<strong>Adjustment details:</strong>",
               tags$div(style = css_box, 
                        paste0("Requested: ", user_threshold, " datasets\n",
                               "Max Found: ", max_freq_found, " datasets\n",
                               "Action: Lowered threshold to ", max_freq_found, "."))),
      
      footer = modalButton("Continue with Adjusted Value")))
    
    return(max_freq_found) }
  
  # --- CASO C: TODO OK ---
  return(user_threshold)}



# ------------------------------------------------------------------------------
#' @9.ENRICHMENT_RESULT_CHECK 
# ------------------------------------------------------------------------------
check_empty_enrichment_results <- function(enrich_res, session) {
  
  if (is.null(enrich_res)) {
    has_rows <- FALSE
  } else {
    df <- tryCatch({
      as.data.frame(enrich_res)
    }, error = function(e) NULL)
    has_rows <- !is.null(df) && nrow(df) > 0}
  
  if (!has_rows) {
    css_box <- "margin-top:4px; background-color:#fff3cd; padding:5px; border-radius:3px; font-family:monospace; font-size:0.85em; color:#856404;"
    
    showModal(modalDialog(
      title = tags$div(icon("search-minus"), " No Significant Results", 
                       style = "color: #e67e22; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "m",
      
      tags$div(class = "alert alert-warning", role = "alert",
               tags$h4(icon("info-circle"), " No Terms Found"),
               tags$p("The enrichment analysis finished, but returned 0 significant terms.")),
      
      tags$div(style = "margin-bottom: 15px;",
               tags$strong("Possible reasons:"),
               tags$ul(style = "margin-top: 5px;",
                       tags$li("Input gene list is too small or contains no annotated genes."),
                       tags$li("Strict filtering (p-value/q-value)."),
                       tags$li("Database mismatch (e.g. wrong Organism selected)."))),
      
      tags$div(style = "margin-top: 10px;",
               tags$strong("Suggestion:"),
               tags$div(style = css_box, 
                        "Try increasing the 'Number of Top Genes' slider or changing the Ontology/Database.")),
      
      footer = modalButton("Close")))
    return(FALSE)}
  return(TRUE)}



# ------------------------------------------------------------------------------
#' @10.CHECK_ANNOTATION_FILE 
# ------------------------------------------------------------------------------
check_annotation_file_format <- function(file_path, session) {
  
  css_box <- "margin-top:4px; background-color:#fff3cd; padding:5px; border-radius:3px; font-family:monospace; font-size:0.85em; color:#856404;"
  
  show_error_modal <- function(title_text, error_msg, actual_preview = NULL) {
    preview_html <- ""
    if (!is.null(actual_preview)) {
      preview_text <- paste(as.character(actual_preview), collapse="\n")
      preview_html <- paste0(
        "<div style='margin-bottom:10px;'><strong>File Preview (First 3 lines):</strong>",
        "<div style='", css_box, "'>", htmltools::htmlEscape(preview_text), "</div></div>")}
    
    showModal(modalDialog(
      title = tags$div(icon("table"), " Annotation File Error", 
                       style = "text-align: center; font-weight: bold; font-size: 18px; color: #c0392b;"),
      size = "l", easyClose = TRUE,
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("exclamation-triangle"), " Format Mismatch"),
               tags$p(HTML(error_msg))),
      HTML(preview_html),
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


# ------------------------------------------------------------------------------
# AUXILIARY: Generate File-Grouped Error Report
# ------------------------------------------------------------------------------
generate_file_grouped_report <- function(invalid_genes, gene_file_map) {
  if (is.null(gene_file_map)) {
    return(paste(head(unique(invalid_genes), 20), collapse = ", "))}
  
  bad_map <- gene_file_map[gene_file_map$GeneID %in% invalid_genes, ]
  if (nrow(bad_map) == 0) return("No file mapping available.")
  file_errors <- list()
  files_split <- strsplit(as.character(bad_map$FileNames), ", ")
  genes_vec   <- bad_map$GeneID
  
  for (i in seq_along(genes_vec)) {
    gene <- genes_vec[i]
    f_vec <- files_split[[i]] 
    
    for (f in f_vec) {
      f <- trimws(f) 
      if (nzchar(f)) {
        file_errors[[f]] <- c(file_errors[[f]], gene)}}}
  
  html_output <- character()
  
  for (f_name in names(file_errors)) {
    genes <- unique(file_errors[[f_name]])
    n_genes <- length(genes)
    
    genes_display <- paste(head(genes, 10), collapse = ", ")
    if (n_genes > 10) genes_display <- paste0(genes_display, ", ... (+", n_genes - 10, ")")
    
    block <- paste0(
      "<div style='margin-bottom: 8px;'>",
      "<strong>", f_name, ":</strong> ",
      "<span style='color:#c0392b; font-family:monospace;'>", genes_display, "</span>",
      "</div>")
    html_output <- c(html_output, block)}
  
  return(paste(html_output, collapse = ""))}



# ------------------------------------------------------------------------------
#' @11.VALIDATE_GENE_IDS
# ------------------------------------------------------------------------------
validate_gene_ids <- function(gene_vector, ID_type, session, gene_file_map = NULL) {
  if (length(gene_vector) == 0) return(gene_vector)
  
  gene_list <- as.character(gene_vector)
  ensembl_pattern <- "^ENS[A-Z]{0,4}[0-9]{6,}$"
  
  pattern <- switch(ID_type, 
                    SYMBOL="^[A-Za-z][A-Za-z0-9._/\\-]*$", 
                    ENTREZID="^[0-9]+$", 
                    ENSEMBL=ensembl_pattern, 
                    stop("Error"))
  
  is_valid <- grepl(pattern, gene_list)
  if (ID_type == "SYMBOL") is_valid <- is_valid & !grepl(ensembl_pattern, gene_list)
  
  valid_genes <- gene_list[is_valid]
  invalid_genes <- gene_list[!is_valid]
  percent_invalid <- (length(invalid_genes) / length(gene_list)) * 100
  
  if (length(invalid_genes) == 0) return(gene_list)
  
  report_html <- generate_file_grouped_report(invalid_genes, gene_file_map)
  
  css_box <- "margin-top:4px; background-color:#fff3cd; padding:10px; border-radius:3px; font-size:0.9em; color:#856404; max-height: 250px; overflow-y: auto;"
  
  # --- CASO 1: ERROR CRITICO (>50%) ---
  if (percent_invalid >= 50) {
    showModal(modalDialog(
      title = tags$div(icon("ban"), " Analysis Stopped: ID Mismatch", style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "l",
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("times-circle"), paste0(round(percent_invalid, 1), "% Invalid IDs")),
               tags$p(paste0("Most genes (>50%) do not match the selected ID type (", ID_type, "). Analysis halted."))),
      
      tags$div(style="margin-bottom:15px;",
               tags$strong("Invalid Genes by File:"),
               HTML(paste0("<div style='", css_box, "'>", report_html, "</div>"))),
      
      tags$div(style="background-color: #f8f9fa; padding: 10px; border-radius: 5px;",
               tags$strong("Expected Formats:"),
               tags$ul(tags$li(tags$b("SYMBOL:"), " e.g., TP53, Brca1"), 
                       tags$li(tags$b("ENTREZID:"), " Digits only (e.g., 7157)"),
                       tags$li(tags$b("ENSEMBL:"), " Starts with 'ENS'"))),
      
      footer = modalButton("Close and Fix Data")))
    
    return(NULL) # Detener análisis
  }
  
  # --- CASO 2: ADVERTENCIA ---
  showModal(modalDialog(
    title = tags$div(icon("exclamation-circle"), " ID Format Warning", style = "color: #e67e22; font-weight: bold; font-size: 18px;"),
    easyClose = TRUE, size = "l",
    tags$div(class = "alert alert-warning", role = "alert",
             tags$h4(icon("info-circle"), paste(length(invalid_genes), "genes will be excluded")),
             tags$p("The analysis will continue, but the following genes were removed because they don't match the ID type.")),
    
    tags$div(style="margin-bottom:15px;",
             tags$strong("Excluded Genes by File:"),
             HTML(paste0("<div style='", css_box, "'>", report_html, "</div>"))),
    
    footer = modalButton("I understand, continue")))
  
  return(valid_genes)}



# ------------------------------------------------------------------------------
#' @12.VALIDATE_ORGANISM 
# ------------------------------------------------------------------------------
validate_organism <- function(gene_list, organism, ID_type, session, gene_file_map = NULL) {
  if (length(gene_list) == 0) return(gene_list)
  
  pattern <- switch(ID_type,
                    SYMBOL = switch(organism, 
                                    Hsa="^[A-Z][A-Z0-9._/\\-]*$", 
                                    Mmu="^[A-Z][a-z0-9._/\\-]*$", 
                                    Rno="^[A-Z][a-z0-9._/\\-]*$", NULL),
                    ENTREZID = "^[0-9]+$",
                    ENSEMBL = switch(organism, 
                                     Hsa="^ENSG[0-9]{11}$", 
                                     Mmu="^ENSMUSG[0-9]{11}$", 
                                     Rno="^ENSRNOG[0-9]{11}$", NULL), NULL)
  
  if (is.null(pattern)) return(gene_list)
  
  is_valid <- grepl(pattern, gene_list)
  valid_genes <- gene_list[is_valid]; invalid_genes <- gene_list[!is_valid]
  percent_invalid <- (length(invalid_genes) / length(gene_list)) * 100
  org_name <- switch(organism, Hsa="Homo sapiens", Mmu="Mus musculus", Rno="Rattus norvegicus")
  
  if (length(invalid_genes) == 0) return(gene_list)
  
  report_html <- generate_file_grouped_report(invalid_genes, gene_file_map)
  
  css_box <- "margin-top:4px; background-color:#fff3cd; padding:10px; border-radius:3px; font-size:0.9em; color:#856404; max-height: 250px; overflow-y: auto;"
  
  # --- CASO 1: ERROR CRITICO (>50%) ---
  if (percent_invalid >= 50) {
    showModal(modalDialog(
      title = tags$div(icon("ban"), " Analysis Stopped: Data Mismatch", style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "l",
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("exclamation-triangle"), paste0(round(percent_invalid, 1), "% invalid genes")),
               tags$p(paste0("Most genes (>50%) do not match the nomenclature for ", org_name, ". Analysis halted."))),
      
      tags$div(style="margin-bottom:15px;",
               tags$strong("Invalid Genes by File:"),
               HTML(paste0("<div style='", css_box, "'>", report_html, "</div>"))),
      
      tags$div(style="background-color: #f8f9fa; padding: 10px; border-radius: 5px;",
               tags$strong("Correct Nomenclature:"),
               tags$ul(tags$li(tags$b("Human:"), " All Uppercase (TP53)"), 
                       tags$li(tags$b("Mouse/Rat:"), " Title Case (Tp53)"))),
      
      footer = modalButton("Close and Fix Data")))
    return(NULL)} 
  
  # --- CASO 2: ADVERTENCIA ---
  showModal(modalDialog(
    title = tags$div(icon("exclamation-circle"), " Nomenclature Warning", style = "color: #e67e22; font-weight: bold; font-size: 18px;"),
    easyClose = TRUE, size = "l",
    tags$div(class = "alert alert-warning", role = "alert",
             tags$h4(icon("info-circle"), paste(length(invalid_genes), "genes will be excluded")),
             tags$p(paste0("The analysis will continue, but some genes were removed because they don't match ", org_name, " nomenclature."))),
    
    tags$div(style="margin-bottom:15px;",
             tags$strong("Excluded Genes by File:"),
             HTML(paste0("<div style='", css_box, "'>", report_html, "</div>"))),
    
    footer = modalButton("I understand, continue")))
  
  return(valid_genes)}



# ------------------------------------------------------------------------------
#' @13.CHECK_DATABASE_EXISTS 
# ------------------------------------------------------------------------------
check_database_exists_metarank <- function(db_type, annot_file) {
  if (db_type == "Personalized" && is.null(annot_file)) {
    
    showModal(modalDialog(
      title = tags$div(icon("times-circle"), " Annotation File Required", 
                       style = "color: #c0392b; font-weight: bold; font-size: 18px;"),
      easyClose = TRUE, size = "m",
      
      tags$div(class = "alert alert-danger", role = "alert",
               tags$h4(icon("file-upload"), " Missing Input"),
               tags$p("The 'Personalized' database option requires a custom annotation file.")),
      
      tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 5px solid #c0392b; margin-top: 15px;",
        tags$h5("Action Required:", style="font-weight:bold; margin-top:0;"),
        tags$ul(tags$li("Please upload a .txt, .tsv, or .csv file in the 'Choose annotation file' input field."))),
      
      footer = modalButton("Close")))
    
    return(FALSE)}
  return(TRUE)}






###############################################################################

#----------------------- AUXILIARY FUNCTIONS ---------------------------------#

###############################################################################

# 1
count_gene_appearance_with_files <- function(gene_data, file_names) {
  if (length(gene_data) != length(file_names)) {
    stop("Length of 'gene_data' and 'file_names' must match.")}
  
  gene_file_tracker <- lapply(seq_along(gene_data), function(i) {
    genes <- gene_data[[i]]$Gene
    
    if (length(genes) == 0) {
      return(NULL)}
    
    data.frame(
      Gene = genes,
      File = file_names[i],
      stringsAsFactors = FALSE)})
  
  gene_file_tracker <- Filter(Negate(is.null), gene_file_tracker)
  
  if (length(gene_file_tracker) == 0) {
    return(data.frame(Gene = character(), File = character(), Count = integer()))}
  
  combined_data <- do.call(rbind, gene_file_tracker)
  
  appearance_summary <- aggregate(
    File ~ Gene,
    data = combined_data,
    FUN = function(x) paste(unique(x), collapse = ", "))
  
  appearance_summary$Count <- sapply(
    strsplit(appearance_summary$File, ", "),
    length)
  
  return(appearance_summary)}



# 2
clean_gene_list <- function(lines) {
  lines <- as.character(lines)          
  lines <- lines[lines != "" & !is.na(lines)]
  lines <- trimws(lines)
  lines <- gsub("///.*|//.*", "", lines)
  unique(lines[lines != ""])}



# 3
read_rra_files <- function(paths) {
  lists <- lapply(paths, function(p) {
    lines <- readLines(p, warn = FALSE)
    clean_gene_list(as.character(lines))})
  names(lists) <- basename(paths)
  lists}



# 4
paste_rra_lists <- function(raw_text) {
  blocks <- unlist(strsplit(raw_text, "###", fixed = TRUE))
  blocks <- blocks[nzchar(trimws(blocks))]
  
  lists <- lapply(blocks, function(b) {
    lines <- unlist(strsplit(b, "\n"))
    clean_gene_list(as.character(lines))})
  
  lists <- Filter(function(x) length(x) > 0, lists)
  
  if (length(lists) == 0) {
    stop("No valid lists found after parsing. Check your ### separators.")}
  
  names(lists) <- paste0("List_", seq_along(lists))
  lists}



# 5
read_rankprod_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  df <- switch(ext,
               csv = read.csv(path, header = FALSE, stringsAsFactors = FALSE),
               tsv = read.delim(path, header = FALSE, stringsAsFactors = FALSE, sep = "\t"),
               txt = read.delim(path, header = FALSE, stringsAsFactors = FALSE, sep = "\t"),
               stop("Unsupported file format: ", path))
  if (ncol(df) < 2) {
    warning("Unexpected format: ", path)
    return(NULL)}
  colnames(df)[1:2] <- c("Gene", "Stat.data")
  df <- df %>%
    filter(!is.na(Gene), Gene != "") %>%
    mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1), Gene = (Gene)) %>%
    distinct(Gene, .keep_all = TRUE)
  buscar_redondeos_v2(df)}



# 6
read_rankprod_files <- function(paths) {
  dfs <- lapply(paths, read_rankprod_file)
  names(dfs) <- basename(paths)
  return(dfs)}



# 7
paste_rankprod_lists <- function(raw_text, format = c("csv", "tsv")) {
  fmt <- match.arg(format)
  blocks <- strsplit(raw_text, "###")[[1]]
  
  dfs <- lapply(blocks, function(b) {
    con <- textConnection(b)
    df <- if (fmt == "csv") read.csv(con, header = TRUE, stringsAsFactors = FALSE)
    else read.delim(con, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    close(con)
    
    df <- df %>%
      filter(!is.na(Gene), Gene != "") %>%
      mutate(
        Gene = as.character(Gene),
        Gene = ifelse(grepl("//", Gene), sapply(strsplit(Gene, "//"), `[`, 1), Gene)
      ) %>%
      distinct(Gene, .keep_all = TRUE)
    
    buscar_redondeos_v2(df)})
  
  names(dfs) <- paste0("List_", seq_along(dfs))
  return(dfs)}



# 8
buscar_redondeos_v2 <- function(tabla) {
  if ("Gene" %in% names(tabla)) {
    tabla$Gene <- (tabla$Gene)}
  for (col in names(tabla)) {
    if (is.numeric(tabla[[col]])) {
      valores_formateados <- format(tabla[[col]], scientific = FALSE)
      indices <- which(
        tabla[[col]] == floor(tabla[[col]]) &
          grepl("\\.0{2,}", valores_formateados))
      if (length(indices) > 0) {
        tabla[[col]][indices] <- NA}}}
  return(tabla)}



# 9
count_gene_appearance_with_files <- function(gene_data, file_names) {
  if (length(gene_data) != length(file_names)) {
    stop("Length of 'gene_data' and 'file_names' must match.")}
  
  gene_file_tracker <- lapply(seq_along(gene_data), function(i) {
    data.frame(
      Gene = gene_data[[i]]$Gene,
      File = file_names[i],
      stringsAsFactors = FALSE)})
  
  combined_data <- do.call(rbind, gene_file_tracker)
  appearance_summary <- aggregate(
    File ~ Gene,
    data = combined_data,
    FUN = function(x) paste(unique(x), collapse = ", "))
  
  appearance_summary$Count <- sapply(
    strsplit(appearance_summary$File, ", "),
    length)
  
  return(appearance_summary)}



# 9
count_gene_appearance_rra <- function(gene_data, file_names) {
  if (length(gene_data) != length(file_names)) {
    stop("Length of 'gene_data' and 'file_names' must be the same.")}
  
  valid_gene_data <- gene_data[sapply(gene_data, length) > 0]
  valid_file_names <- file_names[sapply(gene_data, length) > 0]
  
  if (length(valid_gene_data) == 0) {
    return(data.frame(Gene = character(0), File = character(0), Count = integer(0), stringsAsFactors = FALSE))}
  
  gene_file_tracker <- do.call(rbind, lapply(seq_along(valid_gene_data), function(i) {
    data.frame(
      Gene = valid_gene_data[[i]],
      File = valid_file_names[i],
      stringsAsFactors = FALSE)}))
  
  appearance_summary <- aggregate(
    File ~ Gene,
    data = gene_file_tracker,
    FUN = function(x) paste(unique(x), collapse = ", "))
  
  appearance_summary$Count <- sapply(
    strsplit(appearance_summary$File, ", "),
    length)
  
  return(appearance_summary)}



# 10
generate_statdata_matrix <- function(gene_data, na_management = "ignore", ranking_direction = "ascending") {
  
  all_genes <- unique(unlist(lapply(gene_data, function(dataset) dataset$Gene)))
  statdata_matrix <- sapply(gene_data, function(dataset) {
    stat_values <- setNames(as.numeric(dataset$Stat.data), dataset$Gene)
    stat_values[all_genes]})
  
  statdata_matrix <- as.data.frame(statdata_matrix)
  rownames(statdata_matrix) <- all_genes
  
  if (na_management == "penalize") {
    penalty_value <- if (ranking_direction == "ascending") {
      max(statdata_matrix, na.rm = TRUE)
    } else if (ranking_direction == "descending") {
      min(statdata_matrix, na.rm = TRUE)
    } else {
      stop("Unrecognized ranking_direction; use 'ascending' or 'descending'.")}
    statdata_matrix[is.na(statdata_matrix)] <- penalty_value}
  return(statdata_matrix)}



# 11
perform_RP <- function(logFC_matrix, cl, logged = TRUE, na.rm = TRUE, plot = FALSE,
                       ranking_direction = "ascending", progress = NULL, MinNumOfValidPairs = NA) {
  if (!is.null(progress)) progress(0.1, detail = "Starting RP computation...")
  
  result <- RankProducts(
    logFC_matrix,
    cl = cl,
    na.rm = na.rm,
    logged = logged,
    gene.names = rownames(logFC_matrix),
    plot = plot,
    MinNumOfValidPairs = MinNumOfValidPairs)
  
  if (!is.null(progress)) progress(0.3, detail = "RP computation finished. Processing results...")
  rank_col <- ifelse(ranking_direction == "ascending", 1, 2)
  
  ranked_genes <- data.frame(
    Gene    = rownames(logFC_matrix),
    Rank    = result$RPrank[, rank_col],
    RP_stat = result$RPs[, rank_col],
    PFP     = result$pfp[, rank_col],
    pval    = result$pval[, rank_col])
  
  ranked_genes$p.adjust <- p.adjust(ranked_genes$pval, method = "BH")
  
  ranked_genes <- ranked_genes[order(ranked_genes$Rank), ]
  if (!is.null(progress)) progress(0.8, detail = "Results processed.")
  return(ranked_genes)}



# 12
perform_RP_advance <- function(logFC_matrix, cl, origin, logged = TRUE, na.rm = TRUE, plot = FALSE,
                               ranking_direction = "ascending", progress = NULL) {
  if (!is.null(progress)) progress(0.1, detail = "Starting advanced RP computation...")
  
  result <- RP.advance(
    logFC_matrix,
    cl = cl,
    origin = origin,
    na.rm = na.rm,
    logged = logged,
    gene.names = rownames(logFC_matrix),
    plot = plot)
  
  if (!is.null(progress)) progress(0.6, detail = "Advanced RP computation finished. Processing results...")
  rank_col <- ifelse(ranking_direction == "ascending", 1, 2)
  
  ranked_genes <- data.frame(
    Gene    = rownames(logFC_matrix),
    Rank    = result$RPrank[, rank_col],
    RP_stat = result$RPs[, rank_col],
    PFP     = result$pfp[, rank_col],
    pval    = result$pval[, rank_col])
  
  ranked_genes$p.adjust <- p.adjust(ranked_genes$pval, method = "BH")
  
  ranked_genes <- ranked_genes[order(ranked_genes$Rank), ]
  if (!is.null(progress)) progress(1, detail = "Results processed.")
  return(ranked_genes)}



# 13
run_rra_analysis <- function(gene_data, method = "RRA", full = TRUE) {
  if (!is.list(gene_data) || length(gene_data) == 0) {
    stop("gene_data must be a non-empty list of gene vectors.")}
  
  rra_result <- RobustRankAggreg::aggregateRanks(
    glist   = gene_data,
    N       = NA,
    method  = method,
    full    = full,
    exact   = FALSE,
    topCutoff = NA)
  
  rra_result <- as.data.frame(rra_result)
  rra_result$Ranking <- seq_len(nrow(rra_result))
  colnames(rra_result)[1:2] <- c("Gene", "Score")
  
  # Ajuste BH para Score
  rra_result$p.adjust <- p.adjust(rra_result$Score, method = "BH")
  
  appearance_info <- count_gene_appearance_rra(gene_data, names(gene_data))
  
  final_table <- merge(rra_result, appearance_info, by = "Gene", all.x = TRUE)
  colnames(final_table)[colnames(final_table) == "Count"] <- "Appereance"
  colnames(final_table)[colnames(final_table) == "File"]  <- "Files"
  
  final_table <- final_table[order(final_table$Ranking), 
                             c("Gene", "Ranking", "Score", "p.adjust", "Appereance", "Files")]
  
  colnames(final_table) <- c("GeneID", "Rank", "Score", "p.adjust", "FileCount", "FileNames")
  
  return(final_table)}



# 14
create_combined_rank_table <- function(ranked_genes, gene_appearance) {
  if (!"Gene" %in% colnames(ranked_genes) ||
      !"Gene" %in% colnames(gene_appearance)) {
    stop("Both tables must contain a 'Gene' column.")}
  
  gene_appearance$File <- sapply(strsplit(gene_appearance$File, " "), function(files) {
    paste(basename(files), collapse = " ")})
  merged_data <- merge(ranked_genes, gene_appearance, by = "Gene", all.x = TRUE, sort = FALSE)
  
  merged_data$Rank <- ranked_genes$Rank
  merged_data <- merged_data[, c("Gene","Rank","Count","File","RP_stat","PFP","pval","p.adjust")]
  
  colnames(merged_data) <- c("GeneID", "Rank", "FileCount", "FileNames", 
                             "RP_stat", "PFP", "pvalue", "p.adjust")
  
  return(merged_data)}



# 15
create_combined_rank_table2 <- function(ranked_genes, gene_appearance) {
  if (!"Gene" %in% colnames(ranked_genes) ||
      !"Gene" %in% colnames(gene_appearance)) {
    stop("Both tables must contain a 'Gene' column.")}
  
  total_lists <- max(gene_appearance$Count, na.rm = TRUE)
  max_rank    <- max(ranked_genes$Rank, na.rm = TRUE)
  
  merged_data <- merge(ranked_genes, gene_appearance, by = "Gene", all.x = TRUE, sort = FALSE)
  merged_data$AdjustedRank <- merged_data$Rank +
    ((total_lists - merged_data$Count) * (max_rank / total_lists))
  merged_data <- merged_data[order(merged_data$AdjustedRank), ]
  merged_data$Rank <- seq_len(nrow(merged_data))
  
  merged_data <- merged_data[, c("Gene","Rank","Count","File","RP_stat","PFP","pval","p.adjust")]
  
  colnames(merged_data) <- c("GeneID", "Rank", "FileCount", "FileNames", 
                             "RP_stat", "PFP", "pvalue", "p.adjust")
  
  return(merged_data)}



# 16
create_upset_plot <- function(gene_data,
                              text_size = 16,
                              sets_color = "#56B4E9",
                              intersection_color = "#666666") {
  par(bg = "transparent")
  
  gene_lists <- lapply(gene_data, function(entry) {
    if (is.data.frame(entry)) {
      col <- entry[[1]]  
      if (is.null(col)) return(NULL)
      return(unique(na.omit(as.character(col))))
    } else if (is.vector(entry)) {
      return(unique(na.omit(as.character(entry))))
    } else {
      return(NULL)}})
  
  gene_lists <- Filter(Negate(is.null), gene_lists)
  
  if (is.null(names(gene_lists)) || any(names(gene_lists) == "")) {
    names(gene_lists) <- paste0("List_", seq_along(gene_lists))}
  
  if (length(gene_lists) == 0) stop("No hay listas válidas de genes para el UpSet plot.")
  
  all_genes <- unique(unlist(gene_lists))
  binary_data <- lapply(gene_lists, function(x) as.numeric(all_genes %in% x))
  binary_df <- as.data.frame(do.call(cbind, binary_data))
  
  colnames(binary_df) <- names(gene_lists)
  rownames(binary_df) <- all_genes
  
  binary_df[] <- lapply(binary_df, function(col) {
    suppressWarnings(as.numeric(as.character(col)))})
  
  if (!all(sapply(binary_df, is.numeric))) {
    print(str(binary_df))
    stop("Error: La tabla binaria contiene columnas no numéricas.")}
  
  plot <- UpSetR::upset(
    binary_df,
    sets = names(gene_lists),
    sets.bar.color = sets_color,
    main.bar.color = intersection_color,
    order.by = "freq",
    empty.intersections = "on",
    text.scale = c(
      text_size / 16, text_size / 16,
      text_size / 14, text_size / 14,
      text_size / 12, text_size / 12))
  
  return(plot)}



# 17
save_upset_plot_transparent <- function(plot, output_file, width = 14.5, height = 8.33, dpi = 96) {
  temp_file <- tempfile(fileext = ".png")
  png(temp_file, width = width * dpi, height = height * dpi, res = dpi, bg = "white")
  print(plot)
  dev.off()
  
  img <- image_read(temp_file)
  img_trans <- image_transparent(img, "white", fuzz = 10)
  image_write(img_trans, path = output_file, format = "png")
  message("Plot saved in: ", output_file)}



# 18
generate_heatmap <- function(gene_data, 
                             text_size = 12, 
                             title_size = 17, 
                             axis_title_size = 18, 
                             tick_size = 12, 
                             colorscale = "Viridis") {
  
  gene_lists <- lapply(gene_data, function(entry) {
    if (is.data.frame(entry) && "Gene" %in% colnames(entry)) {
      unique(na.omit(entry$Gene))
    } else if (is.vector(entry)) {
      unique(na.omit(entry))
    } else {
      NULL}})
  
  gene_lists <- Filter(Negate(is.null), gene_lists)
  names(gene_lists) <- names(gene_data)
  
  if (length(gene_lists) < 2) {
    stop("At least two gene lists are required to generate the heatmap.")}
  
  n <- length(gene_lists)
  intersection_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in seq_along(gene_lists)) {
    for (j in seq_along(gene_lists)) {
      inter <- intersect(gene_lists[[i]], gene_lists[[j]])
      uni   <- union(gene_lists[[i]], gene_lists[[j]])
      intersection_matrix[i, j] <- length(inter) / length(uni)}}
  
  rownames(intersection_matrix) <- names(gene_lists)
  colnames(intersection_matrix) <- names(gene_lists)
  
  heatmap_plot <- plot_ly(
    z = intersection_matrix,
    x = colnames(intersection_matrix),
    y = rownames(intersection_matrix),
    type = "heatmap",
    colorscale = colorscale,
    colorbar = list(
      title     = "Proportional<br>Intersection",
      titleside = "right",
      titlefont = list(size = tick_size))) %>%
    
    layout(
      title = list(
        text = "Intersection Proportions Between Gene Lists",
        font = list(size = title_size)),
      xaxis = list(
        title     = "",
        titlefont = list(size = axis_title_size),
        tickfont  = list(size = tick_size)),
      yaxis = list(
        title     = "",
        titlefont = list(size = axis_title_size),
        tickfont  = list(size = tick_size)),
      margin = list(t = 80))
  
  return(heatmap_plot)}



# 19
generate_flat_heatmap <- function(gene_data, 
                                  text_size = 12, 
                                  title_size = 17, 
                                  axis_title_size = 18, 
                                  tick_size = 12, 
                                  colorscale = "Viridis", 
                                  output_file = "heatmap") {
  
  gene_lists <- lapply(gene_data, function(entry) {
    if (is.data.frame(entry) && "Gene" %in% colnames(entry)) {
      unique(na.omit(entry$Gene))
    } else if (is.vector(entry)) {
      unique(na.omit(entry))
    } else {
      NULL}})
  
  gene_lists <- Filter(Negate(is.null), gene_lists)
  names(gene_lists) <- names(gene_data)
  
  if (length(gene_lists) < 2) {
    stop("At least two gene lists are required to generate the heatmap.")}
  
  n <- length(gene_lists)
  intersection_matrix <- matrix(0, ncol = n, nrow = n)
  for (i in seq_along(gene_lists)) {
    for (j in seq_along(gene_lists)) {
      inter <- intersect(gene_lists[[i]], gene_lists[[j]])
      uni   <- union(gene_lists[[i]], gene_lists[[j]])
      intersection_matrix[i, j] <- length(inter) / length(uni)}}
  
  rownames(intersection_matrix) <- names(gene_lists)
  colnames(intersection_matrix) <- names(gene_lists)
  
  df_matrix <- as.data.frame(as.table(intersection_matrix))
  colnames(df_matrix) <- c("Gene_List_1", "Gene_List_2", "Intersection_Score")
  
  if (colorscale == "Viridis") {
    palette <- viridis::viridis(256)
  } else if (colorscale == "Cividis") {
    palette <- viridis::cividis(256)
  } else if (colorscale == "Portland") {
    palette <- RColorBrewer::brewer.pal(9, "YlOrRd")
  } else {
    palette <- viridis::viridis(256)}
  
  heatmap_plot <- ggplot2::ggplot(df_matrix, 
                                  ggplot2::aes(x = Gene_List_1, 
                                               y = Gene_List_2, 
                                               fill = Intersection_Score)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = tick_size, angle = 45, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = tick_size),
      axis.title.x = ggplot2::element_text(size = axis_title_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size),
      plot.title   = ggplot2::element_text(size = title_size, hjust = 0.5)) +
    ggplot2::labs(
      title = "Intersection Proportions Between Gene Lists",
      x     = NULL,
      y     = NULL,
      fill  = "Proportion")
  
  return(heatmap_plot)}



# 20
get_top_genes <- function(combined_results, top_n) {
  req(combined_results)  
  req(top_n > 0)         
  
  top_genes <- head(combined_results$GeneID, n = top_n)
  return(top_genes)}



# 21
get_gene_positions <- function(gene_data_output, ranking_direction = "asc") {
  
  is_rra <- is.character(gene_data_output[[1]])
  processed_lists <- lapply(gene_data_output, function(x) {
    if (is_rra) {
      return(x) 
    } else {
      df <- x
      
      stat_vals <- as.numeric(df[[2]]) 
      genes     <- as.character(df[[1]])
      
      if (ranking_direction == "desc") {
        ord <- order(stat_vals, decreasing = TRUE, na.last = NA)
      } else {
        ord <- order(stat_vals, decreasing = FALSE, na.last = NA)}
      
      return(genes[ord])}})
  
  all_genes <- unique(unlist(processed_lists))
  gene_positions <- lapply(all_genes, function(gene) {
    sapply(processed_lists, function(sorted_list) {
      match(gene, sorted_list)
    })})
  
  result_df <- data.frame(GeneID = all_genes, stringsAsFactors = FALSE)
  result_df$GenePositions <- gene_positions
  result_df$GenePositions <- sapply(result_df$GenePositions, function(pos_vector) {
    paste(ifelse(is.na(pos_vector), "NA", pos_vector), collapse = ",")})
  
  return(result_df)}


###############################################################################

#-------------------------------- NOT USED -----------------------------------#

###############################################################################

Load_gene_file_tsv <- function(files) {
  gene_data <- lapply(files, function(file) {
    data <- read.delim(file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", row.names = NULL)
    
    if (ncol(data) < 2) {
      warning(paste("El archivo", file, "no tiene el formato esperado."))
      return(NULL)}
    
    colnames(data) <- c("Gene", "Stat.data")
    data <- data %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1)) %>%
      mutate(Gene = (Gene)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    data <- buscar_redondeos_v2(data)
    
    return(data)})
  file_names <- basename(files)
  names(gene_data) <- file_names
  return(gene_data = gene_data)}


Load_gene_file_csv <- function(files) {
  gene_data <- lapply(files, function(file) {
    data <- read.csv(file, header = FALSE, stringsAsFactors = FALSE, row.names = NULL)
    
    if (ncol(data) < 2) {
      warning(paste("El archivo", file, "no tiene el formato esperado."))
      return(NULL)}
    
    colnames(data) <- c("Gene", "Stat.data")
    
    data <- data %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1)) %>%
      mutate(Gene = (Gene)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    data <- buscar_redondeos_v2(data)
    
    return(data)})
  
  file_names <- basename(files)
  names(gene_data) <- file_names
  return(list(gene_data = gene_data, file_names = file_names))}


Load_gene_paste_tsv <- function(raw_text) {
  lists <- strsplit(raw_text, split = "###")[[1]]
  gene_data <- lapply(lists, function(list) {
    con <- textConnection(list)
    data <- read.delim(con, header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = NULL)
    close(con)
    data <- data %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1)) %>%
      mutate(Gene = (Gene)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    data <- buscar_redondeos_v2(data)
    
    return(data)})
  return(gene_data)}


Load_gene_paste_csv <- function(raw_text) {
  lists <- strsplit(raw_text, split = "###")[[1]]
  gene_data <- lapply(lists, function(list) {
    con <- textConnection(list)
    data <- read.csv(con, header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
    close(con)
    data <- data %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(Gene = sapply(strsplit(Gene, "//"), `[`, 1)) %>%
      mutate(Gene = (Gene)) %>%
      distinct(Gene, .keep_all = TRUE)
    
    data <- buscar_redondeos_v2(data)
    
    return(data)})
  return(gene_data)}


normalize_symbols <- function(genes, species) {
  switch(species,
         "Hsa" = toupper(genes),
         "Mmu" = sub("^(.)", "\\U\\1", tolower(genes), perl = TRUE),
         "Rno" = sub("^(.)", "\\U\\1", tolower(genes), perl = TRUE),
         genes)}


perform_RS <- function(logFC_matrix, cl, logged = TRUE, na.rm = TRUE, plot = FALSE) {
  result <- RankProducts(
    logFC_matrix,
    cl = cl,
    na.rm = na.rm,
    logged = logged,
    gene.names = rownames(logFC_matrix),
    plot = plot,
    calculateProduct = FALSE)
  ranked_genes <- data.frame(
    Gene = rownames(logFC_matrix),
    Rank = result$RPrank[, 2],
    RS_stat = result$RSs[, 1],
    PFP = result$pfp[, 1],
    pval = result$pval[, 1])
  ranked_genes <- ranked_genes[order(ranked_genes$Rank), ]
  return(ranked_genes)}


perform_RSadvance <- function(logFC_matrix, cl, origin, logged = TRUE, na.rm = TRUE, plot = FALSE) {
  result <- RSadvance(
    logFC_matrix,
    cl = cl,
    origin = origin,
    na.rm = na.rm,
    logged = logged,
    gene.names = rownames(logFC_matrix),
    plot = plot)
  ranked_genes <- data.frame(
    Gene = rownames(logFC_matrix),
    Rank = result$RPrank[, 2],
    RS_stat = result$RSs[, 1],
    PFP = result$pfp[, 1],
    pval = result$pval[, 1])
  ranked_genes <- ranked_genes[order(ranked_genes$Rank), ]
  return(ranked_genes)}


save_combined_rank_table <- function(merged_data, output_file) {
  tsv_file <- paste0(output_file, ".tsv")
  csv_file <- paste0(output_file, ".csv")
  
  write.table(
    merged_data, 
    file = tsv_file, 
    sep = "\t", 
    row.names = FALSE, 
    quote = FALSE)
  
  write.csv(
    merged_data, 
    file = csv_file, 
    row.names = FALSE)}


save_heatmap <- function(gene_data, output_file) {
  heatmap <- generate_heatmap(gene_data)
  htmlwidgets::saveWidget(heatmap, file = output_file)
  message("Heatmap saved in: ", output_file)}


perform_ora_enrichment <- function(gene_list, 
                                   ont = "BP", 
                                   pvalue_cutoff = 0.05, 
                                   qvalue_cutoff = 0.2, 
                                   OrgDb = org.Hs.eg.db) {
  if (length(gene_list) == 0) {
    stop("La lista de genes está vacía. Proporcione genes válidos para el análisis.")}
  
  enrich_result <- enrichGO(
    gene          = gene_list,
    OrgDb         = OrgDb,  
    keyType       = "SYMBOL",  
    ont           = ont,       
    pvalueCutoff  = pvalue_cutoff,
    qvalueCutoff  = qvalue_cutoff,
    readable      = TRUE)
  
  if (!is.null(enrich_result) && nrow(enrich_result@result) > 0) {
    enrich_result@result$pvalue <- round(enrich_result@result$pvalue, digits = 8)
    enrich_result@result$p.adjust <- round(enrich_result@result$p.adjust, digits = 8)
    enrich_result@result$qvalue <- round(enrich_result@result$qvalue, digits = 8)}
  
  return(enrich_result)}


ora_result_tsv <- function(ora_result, file_path) {
  ora_df <- as.data.frame(ora_result)
  write.table(ora_df, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)}


ora_result_csv <- function(ora_result, file_path) {
  ora_df <- as.data.frame(ora_result)
  write.csv(ora_df, file = file_path, row.names = FALSE)}
