################################################################################
#                            GUÍA DE TESTING (QA)                              #
#                          MetaRank - Validaciones                             #
#                                                                              #
#  Este documento describe cómo detonar cada una de las 13 validaciones de la  #
#  aplicación MetaRank para asegurar su correcto funcionamiento.               #
################################################################################

================================================================================
1. PRUEBA DE EXISTENCIA DE INPUT (validate_input_exists)
================================================================================
OBJETIVO: Verificar que la app impide ejecutar el análisis si no hay datos.

PASOS:
1. Abre la aplicación MetaRank.
2. Asegúrate de que las opciones "Use example data?" estén en "No".
3. Selecciona "Input Method" -> "Upload Files".
4. No subas ningún archivo.
5. Deja el resto de los parámetros por defecto/sin cambiar.
6. Haz clic en el botón [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "No Input Data Detected".
  Informa que no se han proporcionado listas de genes.

================================================================================
2. PRUEBA DE MÍNIMO DE LISTAS (validate_min_lists_metarank)
================================================================================
OBJETIVO: El meta-análisis requiere al menos 2 listas para funcionar.

PASOS:
1. Selecciona "Meta-Analysis type" -> "Weighted (Rankprodpack)".
2. En "Input Method" -> "Upload Files".
3. Sube SOLAMENTE UN archivo válido (ej. data1_1.tsv).
4. Deja el resto de los parámetros por defecto/sin cambiar.
5. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: ""Single List Detected".
  Indica que se requiere subir al menos dos archivos o pegar dos listas.

================================================================================
3. PRUEBA DE FORMATO DE ARCHIVO - RRA (validate_file_format_rra)
================================================================================
OBJETIVO: Verificar que los archivos para análisis "Unweighted" tengan el formato correcto.

PASOS:
1. Selecciona "Meta-Analysis type" -> "Unweighted (RRApack)".
2. Sube al menos dos archivos que contengan carácteres como '#' o comas (ej. 
   gene1.txt (válido) y rra_file_error.txt (con errores)).
3. Deja el resto de los parámetros por defecto/sin cambiar.
4. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Format Error: RRA File Input".
  Debe indicar caracteres prohibidos detectados (Header/Comment o Comma).

================================================================================
4. PRUEBA DE FORMATO DE PEGADO - RRA (validate_paste_format_rra)
================================================================================
OBJETIVO: Detectar errores al pegar listas simples (RRA).

PASOS:
1. Selecciona "Meta-Analysis type" -> "Unweighted (RRApack)".
2. Cambia a "Input Method" -> "Paste Genes".
3. Pega el siguiente texto (contiene carácteres no válidos):

FABP4
FMO2;
SLIT3
LDB2
###
FABP4
CLEC3B
STX11
FAM107A
LDB2,
SLIT3

4. Deja el resto de los parámetros por defecto/sin cambiar.
5. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Format Error: RRA Paste Input".
  Indica que se han detectado Hashtags/Headers.

================================================================================
5. PRUEBA DE FORMATO DE ARCHIVO - RANKPROD (validate_file_format_rp)
================================================================================
OBJETIVO: Verificar que los archivos para análisis "Weighted" tengan 2 columnas.
          y n contenga carácteres extraños.

PASOS:
1. Selecciona "Meta-Analysis type" -> "Weighted (Rankprodpack)".
2. Sube al menos dos archivos que NO tenga 2 columnas (ej. data1_1.tsv (válido) 
   y gene1.txt (no válido con una columna)).
3. Deja el resto de los parámetros por defecto/sin cambiar.
4. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Analysis Stopped: Format Error".
  Debe indicar "Column Count Error: Detected 1 columns. Required 2".

================================================================================
6. PRUEBA DE FORMATO DE PEGADO - RANKPROD (validate_paste_format_rp)
================================================================================
OBJETIVO: Detectar errores al pegar datos numéricos (RP).

PASOS:
1. Selecciona "Meta-Analysis type" -> "Weighted (Rankprodpack)".
2. Cambia a "Input Method" -> "Paste Genes".
3. Selecciona formato "CSV".
4. Pega el siguiente texto:
   
FABP4,0.2
FMO2,0.1
SLIT3,0.6
LDB2,0.21
###
FABP4,0.2
CLEC3B,0.25
STX11,0.3;
FAM107A,0.85
LDB2,,0.1
SLIT3,0.3

5. Deja el resto de los parámetros por defecto/sin cambiar.
6. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Amarillo o Rojo: "Data Warning" o "Analysis Stopped".
  Debe advertir sobre "Incomplete Rows".

================================================================================
7. PRUEBA DE ORIGEN - RANKPROD ADVANCED (validate_origin_rp)
================================================================================
OBJETIVO: Verificar consistencia del vector de origen (longitud y replicados).

PASOS:
1. Activa la opción de use "Use Example Data"
2. En "Parameters", en "Rank-Based Meta-Analysis Methods"" cambia a "RankProd advance".
3. En el campo "Origin", escribe: 1,1
   (Error: Tenemos 4 archivos, pero solo 2 valores de origen).
4. Deja el resto de los parámetros por defecto/sin cambiar.
5. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Length Mismatch".
  Indica discrepancia entre número de archivos y valores de origen.

*Nota adicional: Probar también "1,1,2,1" para detonar "Insufficient Replicates".*
o dejar en blanco el campo. 

================================================================================
8. PRUEBA DE THRESHOLD SIN INTERSECCIÓN (validate_gene_threshold)
================================================================================
OBJETIVO: Verificar qué pasa si las listas no tienen NINGÚN gen en común.

PASOS:
1. Selecciona "Unweighted (RRApack)" y "Paste Genes".
2. Pega el siguiente conjunto de listas:

TP53
EGFR
MYC
PTEN
###
TP53
EGFR
MYC
###
TP53
GAPDH
BDNF
INS
###
VEGFA
SOD1
NOS1
CREB1

3. Deja el resto de los parámetros por defecto salvo "Minimum number of 
   datasets", dejándolo en 4
4. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Threshold Too Strict".
  Informa que el análisis es imposible (Max overlap found: 3).
  
  *Nota adicional: Probar también dos listas sin ningún gen común, saltando 
  el error de "No Shared Genes Found".

================================================================================
9. PRUEBA DE RESULTADOS DE ENRIQUECIMIENTO VACÍOS (check_empty_enrichment_results)
================================================================================
OBJETIVO: Detectar fallo de mapeo o falta de resultados significativos tras el análisis.

PASOS:
1. Selecciona "Meta-Analysis type" -> "Weighted (Rankprodpack)".
2. Cambia a "Input Method" -> "Paste Genes".
3. Selecciona formato "CSV".
4. Pega el siguiente texto:

FAKEGENE1,0.01
FAKEGENE2,0.05
###
FAKEGENE1,0.02
FAKEGENE2,0.04

5. Activa "Perform Enrichment Analysis".
6. Deja el resto de los parámetros por defecto
7. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Naranja: "No Significant Results".
  Indica que el análisis finalizó pero devolvió 0 términos ("No Terms Found").


================================================================================
10. PRUEBA DE EXISTENCIA DE BASE DE DATOS (check_database_exists_metarank)
================================================================================
OBJETIVO: Verificar que se exige archivo para la base de datos personalizada.

PASOS:
1. Activa la opción de use "Use Example Data"
2. Activa "Perform Enrichment Analysis".
3. Selecciona "ORA Database" -> "Personalized".
4. NO subas ningún archivo en "Choose annotation file".
5. Deja el resto de los parámetros por defecto
6. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Annotation File Required".

================================================================================
11. PRUEBA DE FORMATO DE ANOTACIÓN (check_annotation_file_format)
================================================================================
OBJETIVO: El archivo de anotación personalizada debe tener 2 columnas.

PASOS:
1. Activa la opción de use "Use Example Data"
2. Activa "Perform Enrichment Analysis".
3. Activa Enriquecimiento -> Database "Personalized".
4. Sube un archivo que tenga 3 columnas o solo 1 (ej. gene1.txt).
5. Deja el resto de los parámetros por defecto
6. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Annotation File Error".
  Indica "Format Mismatch".

================================================================================
12. PRUEBA DE FORMATO DE ID DE GENES (validate_gene_ids)
================================================================================
OBJETIVO: Verificar que los IDs coinciden con el tipo seleccionado (Entrez vs Symbol).

PASOS:
1. Activa la opción de use "Use Example Data"
2. Activa "Perform Enrichment Analysis".
3. Deja el resto de los parámetros por defecto salvo "Gene ID".
4. En "Gene ID", selecciona "ENTREZID" (que espera números, mientras que nuetros
   datos tienen SYMBOL).
5. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Analysis Stopped: ID Mismatch".
  
*Nota importante*: Si se hay menos de 50% de genes en formato incorrecto, el
análisis continua ignorando esos genes. En caso de que se detecte >50% de error 
El análisis se detendrá por completo.

================================================================================
13. PRUEBA DE NOMENCLATURA DE ORGANISMO (validate_organism)
================================================================================
OBJETIVO: Verificar mayúsculas/minúsculas según la especie.

PASOS:
1. Activa la opción de use "Use Example Data"
2. Activa "Perform Enrichment Analysis".
3. Deja el resto de los parámetros por defecto salvo "Organism".
4. En "Organism", selecciona "Mus musculus" (nuestros datos son para Homo sapiens).
5. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Analysis Stopped: Data Mismatch".

*Nota importante*: Si se hay menos de 50% de genes en formato incorrecto, el
análisis continua ignorando esos genes. En caso de que se detecte >50% de error 
El análisis se detendrá por completo.

