# Set the output directory for annotation files
output_dir <- "./database_annotations"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load required libraries
library(org.Hs.eg.db)     # Human gene annotation database
library(org.Mm.eg.db)     # Mouse gene annotation database
library(org.Rn.eg.db)     # Rat gene annotation database
library(org.Dm.eg.db)     # (Optional) Drosophila gene annotation database
library(GO.db)            # Gene Ontology database
library(AnnotationDbi)    # Tools for annotation databases
library(ReactomePA)       # Reactome pathway analysis
library(biomaRt)          # Interface to Ensembl BioMart
library(KEGGREST)         # Interface to KEGG REST API
library(reactome.db)



################################################################################
#                  FUNCTIONS TO GENERATE GO ANNOTATIONS                        #
################################################################################


# Function to select the appropriate OrgDb object based on the organism code
#----------------------------------------------------------------------------
getOrgDb <- function(organism) {
  if (organism == "Hsa") {
    return(org.Hs.eg.db)
  } else if (organism == "Mmu") {
    return(org.Mm.eg.db)
  } else if (organism == "Rno") {
    return(org.Rn.eg.db)
  } else {
    stop("Unsupported organism code")}}



# Function to generate and save GO annotation files for a given organism and ID type
#------------------------------------------------------------------------------------
generate_GO_files <- function(organism, ID) {
  
  orgdb <- getOrgDb(organism)
  keys_list <- keys(orgdb, keytype = ID)                                          # Retrieve the list of keys (gene identifiers) for the given ID type
  go_data <- AnnotationDbi::select(orgdb,
                                   keys = keys_list,
                                   columns = c(ID, "GO", "ONTOLOGY"),             # Extract GO annotations: gene ID, GO term, and ontology (BP/MF/CC)
                                   keytype = ID)
  go_data <- na.exclude(go_data)                                                  # Remove rows without GO annotation
  
  annotations <- list(
    BP = unique(go_data[go_data$ONTOLOGY == "BP", c("GO", ID)]),
    MF = unique(go_data[go_data$ONTOLOGY == "MF", c("GO", ID)]),                  # Split GO annotations by ontology category and remove duplicates
    CC = unique(go_data[go_data$ONTOLOGY == "CC", c("GO", ID)]))
  
  filename_BP <- file.path(output_dir, paste0(organism, "_GO_BP_", ID, ".txt"))
  filename_MF <- file.path(output_dir, paste0(organism, "_GO_MF_", ID, ".txt"))   # Define output filenames for each ontology category
  filename_CC <- file.path(output_dir, paste0(organism, "_GO_CC_", ID, ".txt"))
  
  write.table(annotations$BP, filename_BP, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(annotations$MF, filename_MF, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)     # Write the annotation files
  write.table(annotations$CC, filename_CC, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  cat("Files generated for", organism, "with ID type", ID, "\n")}



# Function to generate and save the GO term descriptions (TERM2NAME file)
#------------------------------------------------------------------------
generate_term2name <- function() {
  go_terms <- Term(GOTERM)                                                        # Get GO term descriptions
  term2name <- data.frame("GO" = names(go_terms), "NAME" = go_terms)
  out_fn <- file.path(output_dir, "TERM2NAME_GO.txt")
  write.table(term2name, out_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat("TERM2NAME_GO.txt file generated\n")}



#-----------------------------------------------------------------------------#
#                              GO SCRIPT EXECUTION                            #
#-----------------------------------------------------------------------------#

# Define the list of organisms (by code) and the types of gene IDs to process
organisms <- c("Hsa", "Mmu", "Rno")
ID_types <- c("ENTREZID", "SYMBOL", "ENSEMBL")

# Loop through each organism and each ID type to generate GO annotation files
for (org in organisms) {
  orgdb <- getOrgDb(org)
  available_IDs <- keytypes(orgdb)
  for (id in ID_types) {
    if (id %in% available_IDs) {
      generate_GO_files(org, id)
    } else {
      cat("ID type", id, "is not available for organism", org, "\n")}}}

# Generate the TERM2NAME file for GO terms
generate_term2name()














################################################################################
#                FUNCTIONS TO GENERATE KEGG ANNOTATIONS                        #
################################################################################

# Function to select the appropriate OrgDb object based on the organism code
#----------------------------------------------------------------------------
getGeneNameDb <- function(organism) {
  if (organism == "Hsa") {
    return(org.Hs.egGENENAME)
  } else if (organism == "Mmu") {
    return(org.Mm.egGENENAME)
  } else if (organism == "Rno") {
    return(org.Rn.egGENENAME)
  } else {
    stop("Unsupported organism code")}}


# Function to generate and save KEGG pathway annotation files for a given organism and ID type
#----------------------------------------------------------------------------------------------
generate_KEGG_files <- function(organism, ID) {
  orgdb <- getOrgDb(organism)
  
  if (!(ID %in% keytypes(orgdb))) {
    cat("ID type", ID, "is not available for", organism, "\n")                    # Check if the requested ID type is available for the organism
    return(NULL)}

  gene_keys <- keys(orgdb, keytype = ID)                                          # Retrieve gene identifiers of the selected type
  kegg_data <- AnnotationDbi::select(orgdb,
                                     keys = gene_keys,
                                     columns = c(ID, "PATH"),                     # Extract KEGG pathway annotations (ID and PATH)
                                     keytype = ID)
  
  kegg_data <- kegg_data[, c("PATH", ID)]                                         # Keep only the columns of interest: pathway and gene ID
  kegg_data <- na.exclude(kegg_data)                                              # Remove entries with missing annotations
  filename <- file.path(output_dir, paste0(organism, "_KEGG_", ID, ".txt"))       # Define the output filename
  write.table(kegg_data, filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)    # Write the KEGG annotation file
  
  cat("File generated:", filename, "\n")}


# Function to generate and save gene ID to name mapping files (ID2NAME)
#-----------------------------------------------------------------------
generate_geneName_file <- function(organism) {
  geneNameDb <- getGeneNameDb(organism)
  id2name <- toTable(geneNameDb)
  filename <- file.path(output_dir, paste0("ID2NAME_", organism, ".txt"))                          # Define the output filename (one file per organism)
  write.table(id2name, filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)   # Write the mapping file
  cat("File generated:", filename, "\n")}


generate_KEGG_pathway_desc <- function(organism) {
  # e.g. organism = "Hsa" → code_kegg = "hsa"
  code_kegg <- tolower(organism)
  # obtenemos todas las vías KEGG para el organismo
  pw_list <- KEGGREST::keggList("pathway", code_kegg)
  # keggList devuelve un named vector: names son "path: hsa00010", valores "Glycolysis / Gluconeogenesis – Homo sapiens (human)"
  df <- data.frame(
    pathway_id   = sub("^path:", "", names(pw_list)),
    description  = pw_list,
    stringsAsFactors = FALSE
  )
  # Quitamos el sufijo " – Homo sapiens (human)" (u otros) si queremos solo la descripción pura
  df$description <- sub(" –.*$", "", df$description)
  # Guardamos
  filename <- file.path(output_dir, paste0(organism, "_KEGG_pathway_desc.txt"))
  write.table(df, filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  cat("File generated:", filename, "\n")
}

generate_general_Kegg_map <- function() {
  # 1) traemos SOLO las vías de referencia
  pw_map <- KEGGREST::keggList("pathway", "map")
  
  # 2) construimos el data.frame
  df <- data.frame(
    # quitamos el prefijo "path:" si existiera, luego el "map"
    pathway_id = sub("^map", "",
                     sub("^path:", "", names(pw_map))),
    description = sub(" –.*$", "", pw_map),
    stringsAsFactors = FALSE
  )
  
  # 3) escribimos en TERM2NAME_KEGG.txt
  write.table(df,
              file.path(output_dir, "TERM2NAME_KEGG.txt"),
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  cat("File generado:", file.path(output_dir, "TERM2NAME_KEGG.txt"), "\n")
}



#-----------------------------------------------------------------------------#
#                            KEGG SCRIPT EXECUTION                            #
#-----------------------------------------------------------------------------#

organisms <- c("Hsa", "Mmu", "Rno")
ID_types <- c("ENTREZID", "SYMBOL", "ENSEMBL")

# Loop through each organism and ID type to generate KEGG annotation files
for (org in organisms) {
  for (id in ID_types) {
    generate_KEGG_files(org, id)}
  
  # Generate gene name mapping file
  generate_geneName_file(org)
  generate_KEGG_pathway_desc(org)}

generate_general_Kegg_map()












################################################################################
#                FUNCTIONS TO GENERATE REACTOME ANNOTATIONS                    #
################################################################################

# Mapping of organism codes to Ensembl dataset names
species_mart <- list(
  Hsa = "hsapiens_gene_ensembl",
  Mmu = "mmusculus_gene_ensembl",
  Rno = "rnorvegicus_gene_ensembl")

# ID types to process and their corresponding biomaRt attributes
ID_types <- c("ENTREZID", "SYMBOL", "ENSEMBL")
attr_map <- list(
  ENTREZID = "entrezgene_id",
  SYMBOL   = "external_gene_name",
  ENSEMBL  = "ensembl_gene_id")


# Function to generate Reactome annotation files for a given organism and ID type
#---------------------------------------------------------------------------------
generate_Reactome_files <- function(org_code, mart_dataset, ID_type) {
  message("Processing ", org_code, " - ", ID_type, " ...")

  mart <- useEnsembl(biomart = "genes", dataset = mart_dataset)                              # Connect to Ensembl BioMart
  df <- getBM(
    attributes = c(attr_map[[ID_type]], "reactome"),                                         # 1) Retrieve gene-to-pathway associations
    mart       = mart)
  df <- df[!is.na(df[,1]) & !is.na(df[,2]), , drop=FALSE]                                    # Remove rows with missing gene ID or Reactome pathway ID
  
  df2 <- data.frame(
    term = df[, 2],
    gene = df[, 1],                                                                          # 2) Reformat columns: 'term' = pathway ID, 'gene' = gene ID
    stringsAsFactors = FALSE)
  
  out_fn <- file.path(output_dir, paste0(org_code, "_Reactome_", ID_type, ".txt"))           # 3) Write to file with appropriate headers
  write.table(df2, out_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message("  -> File generated: ", out_fn)}



# Function to generate TERM2NAME_REACTOME.txt: mapping from Reactome pathway ID to name
#---------------------------------------------------------------------------------------
generate_term2name_reactome <- function() {
  message("Generating TERM2NAME_REACTOME.txt ...")
  
  if (!requireNamespace("reactome.db", quietly = TRUE))                                          # Ensure the reactome.db package is installed
    BiocManager::install("reactome.db")

  all_ids <- keys(reactome.db, keytype = "PATHID")                                               # Retrieve all Reactome pathway IDs
  t2n <- AnnotationDbi::select(
    reactome.db,
    keys    = all_ids,
    columns = c("PATHID", "PATHNAME"),                                                           # Get pathway names corresponding to each ID
    keytype = "PATHID")
  colnames(t2n) <- c("term", "name")                                                             # Rename columns to match expected format
  
  out_fn <- file.path(output_dir, "TERM2NAME_REACTOME.txt")
  write.table(t2n, out_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message("  -> File generated: ", out_fn)}


#-----------------------------------------------------------------------------#
#                        REACTOME SCRIPT EXECUTION                            #
#-----------------------------------------------------------------------------#

# Loop through all organisms and ID types to generate Reactome annotation files
for (org in names(species_mart)) {
  ds <- species_mart[[org]]
  for (id in ID_types) {
    tryCatch(generate_Reactome_files(org, ds, id),
      error = function(e) message("  ! Error with ", org, "/", id, ": ", e$message))}}

# Generate Reactome term-to-name mapping file
generate_term2name_reactome()

