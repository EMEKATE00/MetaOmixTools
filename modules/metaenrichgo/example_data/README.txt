################################################################################
#                            GUÍA DE TESTING (QA)                              #
#                       MetaEnrichGO - Validaciones                            #
#                                                                              #
#  Este documento describe cómo detonar cada una de las validaciones de la     #
#  aplicación para asegurar su correcto funcionamiento y visualización.        #
################################################################################

================================================================================
1. PRUEBA DE EXISTENCIA DE INPUT (validate_input_exists)
================================================================================
OBJETIVO: Verificar que la app impide ejecutar el análisis si no hay datos de input

PASOS:
1. Abre la aplicación (o recarga la página).
2. Asegúrate de que la opción "Use example data?" esté en "No".
3. No subas ningún archivo en "Upload Files".
4. No pegues nada si cambias a "Paste Genes".
5. No modifiques nada de opciónes (aunque modificandolas el efecto es el mismo)
6. Haz clic en el botón [Run Analysis].

RESULTADO ESPERADO:
- Modal sencillo: "No Input Data Detected"
                  Informa sobre la detención del análisis por la falta del archivo input 
                  y la recomendación de como proceder a continuación.

================================================================================
2. PRUEBA DE INPUT PERSONALIZADO FALTANTE (check_database_exists)
================================================================================
OBJETIVO: Verificar que, al elegir base de datos personalizada, se exige un archivo.

PASOS:
1. Sube cualquier archivo o pega un conjunto VÁLIDO de genes (o activa "Use example data").
2. En "Parameters" > "Database", selecciona "Personalized".
3. NO subas ningún archivo en el campo "Choose annotation file" que aparece.
4. No tocar el resto de los parámemtros (aunque no tienen efecto su cambio)
5. Haz clic en el botón [Run Analysis].

RESULTADO ESPERADO:
- Modal sencillo: "Annotation File Required".
                  Informa sobre la detención del análisis por la falta de archivo en
                  en el campo de Database -> Personalized + recomendación de como proceder

================================================================================
3. PRUEBA DE FORMATO DE ARCHIVO INCORRECTO (validate_file_format_metaenrichgo)
================================================================================
OBJETIVO: Detectar caracteres prohibidos (comas, tabulaciones, headers) en archivos .txt.

PASOS:
1. Localizar el archivo 'v3.txt' dentro de la carpeta example_data.
2. Abrir el archivo con el editor de preferencia y comprobar la presencia de comas (,),
   tabuladores (\t), punto y coma (;) o hastags (#).
3. Subir el archivo en el parametro Input -> Input method -> Upload Files.
4. Dejar el resto de los parametros por defecto (Database GO -> BP, Organism -> Homo 
   sapiens, GeneID -> SYMBOL, el resto no importa).
5. Haz clic en el botón [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Format Error: File Input".
              Informa sobre todos errores encontrados indicando el archivo y el tipo de error
              además del formato recomendado.

================================================================================
4. PRUEBA DE FORMATO DE PEGADO INCORRECTO (validate_paste_format)
================================================================================
OBJETIVO: Detectar copiado/pegado desde Excel (con columnas ocultas) o TXT.

PASOS:
1. Localizar el archivo 'v3.txt' dentro de la carpeta example_data.
2. Abrir el archivo con el editor de preferencia y comprobar la presencia de comas (,),
   tabuladores (\t), punto y coma (;) o hastags (#).
3. Copiar y pegar el contenido del archivo en el parametro Input -> Input method -> Paste Genes.
4. Dejar el resto de los parametros por defecto (Database GO -> BP, Organism -> Homo 
   sapiens, GeneID -> SYMBOL, el resto no importa).
5. Haz clic en el botón [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Format Error: Paste Input".
              Informa sobre todos errores encontrados indicando el archivo y el tipo de error
              además del formato recomendado.

================================================================================
5. PRUEBA DE FORMATO DE ANOTACIÓN PERSONALIZADA (check_annotation_file_format)
================================================================================
OBJETIVO: Verificar que el archivo de anotación tenga exactamente 2 columnas.

PASOS:
1. Sube una lista de genes válida (o usa example data).
2. Selecciona Database: "Personalized".
3. Sube 'v3.txt' en el campo de anotación personalizada.
4. Deja el resto de los parámetros por defecto
5. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Annotation File Error".
              Informa sobre el error específico, hace una previsualización del archivo
              subido e informa sobre el formato recomendado

================================================================================
6. PRUEBA DE FORMATO DE ID (validate_gene_ids)
================================================================================
OBJETIVO: Verificar que los identificadores coinciden con el tipo seleccionado.

CASO A: ERROR MASIVO (Bloqueante)
1. En "Input method", selecciona "Paste Genes".
2. Pega esta lista de IDs de Ensembl:

    Actb
    Gapdh
    Pten
    Tp53
    ###
    Foxp2
    Drd2
    Cyp2e1
    Nos1
    ENSG00000139618

3. En "Parameters" > "Gene ID", selecciona "ENTREZID".
4. En "Parameters" > "Organism", selecciona "Mus musculus" o "Rattus norvegicus".
5. Haz clic en [Run Analysis].
-> RESULTADO: Modal Rojo "Analysis Stopped". Muestra % de error (>50%).

CASO B: ADVERTENCIA (No bloqueante)
1. Pega una lista mixta (9 genes Entrez correctos y 1 incorrecto):

    Actb
    Gapdh
    Pten
    Tp53
    ###
    Foxp2
    Drd2
    Cyp2e1
    Nos1
    ENSG00000139618

2. Selecciona Gene ID: "SYMBOL".
3. En "Parameters" > "Organism", selecciona "Mus musculus" o "Rattus norvegicus".
4. Haz clic en [Run Analysis].
-> RESULTADO: Modal Naranja "ID Format Warning". Permite continuar ("I understand, continue analysis").

================================================================================
7. PRUEBA DE NOMENCLATURA DE ORGANISMO (validate_organism)
================================================================================
OBJETIVO: Verificar mayúsculas/minúsculas según el organismo (Humano vs Ratón).

CASO A: ERROR MASIVO (Bloqueante)
1. En "Input method", selecciona "Paste Genes".
2. Pega esta lista de IDs de Ensembl:

    Actb
    Gapdh
    Pten
    Tp53
    ###
    Foxp2
    Drd2
    Cyp2e1

3. En "Parameters" > "Gene ID", selecciona "SYMBOL".
4. En "Parameters" > "Organism", selecciona "Homo sapiens".
5. Haz clic en [Run Analysis].
-> RESULTADO: Modal Rojo "Analysis Stopped". Muestra % de error (>50%).

CASO B: ADVERTENCIA (No bloqueante)
1. Pega una lista mixta (9 genes Entrez correctos y 1 incorrecto):

    Actb
    Gapdh
    Pten
    Tp53
    ###
    Foxp2
    Drd2
    Cyp2e1
    NOS1

2. En "Parameters" > "Gene ID", selecciona "SYMBOL".
3. En "Parameters" > "Organism", selecciona "Mus musculus" o "Rattus norvegicus".
4. Haz clic en [Run Analysis].
-> RESULTADO: Modal Naranja "ID Format Warning". Permite continuar ("I understand, continue analysis").


================================================================================
8. PRUEBA DE FALLO DE MAPEO (check_mapping_failure)
================================================================================
OBJETIVO: Como resultado del análisis, no se han registrado resultados mapedos contra la base de datos.
          El formato es válido (texto limpio), pero los genes no existen en la DB.

CASO A ---
1. Selecciona Organism: "Homo sapiens", ID: "SYMBOL", DB: "GO".
2. Pega genes inventados que cumplan el formato de texto (letras):
    FAKEGENEAA
    NOMAPGENEBB
    RANDOMGENECC
3. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Gene Mapping Failure".

CASO B ---
1. Activa los datos de ejemplo.
2. Selecciona "Parameters" > "Database" > "Personalized"
3. Selecciona el archivo v4.txt disponible en la carpeta example_data
4. Deja el resto de los parametros por defecto
5. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "Gene Mapping Failure".

================================================================================
9. PRUEBA FILTRO MUY ESTRICTO (check_strict_threshold)
================================================================================
OBJETIVO: Informar sobre la ausencia de datos obtenidos trás el análisis debido
          un filtro de minimum treshold demasiado estrico.

PASOS:
1. Sube una lista de genes válida (o usa example data).
2. Selecciona Database: "Personalized".
3. Sube 'Chronic.tsv' en el campo de anotación personalizada.
4. Deja el resto de los parámetros por defecto
5. Haz clic en [Run Analysis].

RESULTADO ESPERADO:
- Modal Rojo: "No Significant Results Found".
              Informa sobre el error específico y hace una recomendación de como actuar