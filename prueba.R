generar_dummy_valido <- function(nombre_archivo, megas_objetivo) {
  message(paste("Generando", nombre_archivo, "con formato válido..."))
  
  # Abrir conexión de escritura
  con <- file(nombre_archivo, "w")
  
  # Definir una línea base que parezca un dato real: "Gene_NAME \t Pvalue \n"
  # Hacemos la línea un poco larga para llegar rápido al peso deseado
  # Ejemplo: Gene_00000001    0.00001 ... (relleno)
  relleno <- paste(rep("A", 50), collapse = "") 
  linea_base <- paste0("Gene_%d_", relleno, "\t0.05\n")
  
  # Calcular bytes por línea aproximados
  bytes_por_linea <- nchar(sprintf(linea_base, 1))
  total_lineas <- ceiling((megas_objetivo * 1024^2) / bytes_por_linea)
  
  # Escribir en bloques para velocidad (10,000 líneas por bloque)
  bloque_size <- 10000
  num_bloques <- floor(total_lineas / bloque_size)
  resto <- total_lineas %% bloque_size
  
  # Escribir bloques
  for (i in 1:num_bloques) {
    # Creamos un vector de índices para este bloque
    indices <- ((i-1)*bloque_size + 1) : (i*bloque_size)
    # Generamos el bloque de texto
    texto_bloque <- paste0(sprintf(linea_base, indices), collapse = "")
    cat(texto_bloque, file = con)
    
    # Barra de progreso simple en consola
    if(i %% 10 == 0) cat(".")
  }
  
  # Escribir el resto
  if (resto > 0) {
    indices <- (num_bloques*bloque_size + 1) : total_lineas
    texto_resto <- paste0(sprintf(linea_base, indices), collapse = "")
    cat(texto_resto, file = con)
  }
  
  close(con)
  cat("\n")
  message(paste("✅ Archivo creado:", nombre_archivo, 
                "| Tamaño final:", round(file.size(nombre_archivo)/1024^2, 2), "MB"))
}

# --- EJECUTAR ---
# 1. Archivo de ~29 MB (Debería subir bien)
generar_dummy_valido("valid_29mb.txt", 26)

# 2. Archivo de ~31 MB (Debería ser rechazado por la barra de carga de Shiny)
generar_dummy_valido("valid_31mb.txt", 31)




# Función para crear archivos de texto dummy de un tamaño específico (en MB)
crear_archivo_dummy <- function(nombre, megas) {
  message(paste("Generando", nombre, "..."))
  
  bytes_objetivo <- megas * 1024^2
  linea <- paste0(paste(rep("X", 1023), collapse = ""), "\n") # 1 KB por línea aprox
  con <- file(nombre, "w")
  
  iteraciones <- ceiling(bytes_objetivo / nchar(linea))
  
  # Escribimos en bloques para ser eficientes
  bloque <- paste(rep(linea, 1000), collapse = "") # Bloque de ~1MB
  bloques_completos <- floor(iteraciones / 1000)
  resto <- iteraciones %% 1000
  
  for(i in 1:bloques_completos) {
    cat(bloque, file = con)
  }
  
  if(resto > 0) {
    cat(paste(rep(linea, resto), collapse = ""), file = con)
  }
  
  close(con)
  print(paste("Creado:", nombre, "| Tamaño:", file.size(nombre)/1024^2, "MB"))
}

# Generar los archivos
crear_archivo_dummy("test_29mb.txt", 27)
crear_archivo_dummy("test_31mb.txt", 31)


