#!/bin/bash

# 01_download_muonia2011A_AOD_rootfiles.sh
#
# Purpose:
#   Interactive downloader for CMS Run2011A MuOnia AOD ROOT files using XRootD.
#   The script:
#     - Reads file-index lists from the directory `muonia_links/`
#       (each index contains one XRootD URL per line).
#     - Lets the user select which index file to use and which line range to process.
#     - Downloads the selected files into `root_files/`, renaming them as
#       `muonia_<slice>_<index>.root`.
#     - Creates a timestamped log in `root_files/logs_downloads/` summarising
#       successful, failed, and skipped downloads.
#
# Directories:
#   LINK_DIR : location of the file-index .txt files (XRootD URLs).
#   ROOT_DIR : destination for downloaded ROOT files.
#   LOG_DIR  : destination for download logs.

# Path to the directory containing the txt files
LINK_DIR="/code/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor/muonia_links"
ROOT_DIR="/code/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor/root_files"
LOG_DIR="/code/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor/root_files/logs_downloads"
mkdir -p "$ROOT_DIR"
mkdir -p "$LOG_DIR"

# List of txt files (in a manually defined order)
LINK_FILES=(
    "CMS_Run2011A_MuOnia_AOD_12Oct2013-v1_00000_file_index.txt"
    "CMS_Run2011A_MuOnia_AOD_12Oct2013-v1_00001_file_index.txt"
    "CMS_Run2011A_MuOnia_AOD_12Oct2013-v1_20000_file_index.txt"
    "CMS_Run2011A_MuOnia_AOD_12Oct2013-v1_210000_file_index.txt"
)

echo "📁 Archivos disponibles para descarga:"
for i in "${!LINK_FILES[@]}"; do
    full_path="${LINK_DIR}/${LINK_FILES[$i]}"
    if [ -f "$full_path" ]; then
        num_lines=$(wc -l < "$full_path")
        echo "  $((i+1)). ${LINK_FILES[$i]} ($num_lines líneas)"
    else
        echo "  $((i+1)). ${LINK_FILES[$i]} (❌ archivo no encontrado)"
    fi
done


# File selection
read -p "Selecciona un número (1-4) según el archivo deseado: " FILE_CHOICE

if ! [[ "$FILE_CHOICE" =~ ^[1-4]$ ]]; then
    echo "❌ Selección inválida. Debe ser un número del 1 al 4."
    exit 1
fi

# Determine txt file and prefix for ROOT file name
TXT_FILE="${LINK_DIR}/${LINK_FILES[$((FILE_CHOICE-1))]}"
PREFIX="muonia_${FILE_CHOICE}_"

# Check existence
if [ ! -f "$TXT_FILE" ]; then
    echo "❌ Error: No se encontró el archivo '$TXT_FILE'"
    exit 1
fi

# Ask for the range
echo "🔢 ¿Desde qué número de archivo deseas comenzar la descarga?"
read -p "Inicio: " START

echo "🔢 ¿Hasta qué número de archivo deseas terminar la descarga?"
read -p "Fin: " END

if ! [[ "$START" =~ ^[0-9]+$ && "$END" =~ ^[0-9]+$ && "$START" -le "$END" ]]; then
    echo "❌ Error: Ingresa dos números válidos (y que el inicio sea menor o igual que el fin)."
    exit 1
fi

# Start timestamp
START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
echo "📋 Registro de descargas del $START al $END del archivo ${LINK_FILES[$((FILE_CHOICE-1))]}" | tee "$LOG_DIR/temp_log.txt"
echo "🕒 Inicio: $START_TIME" | tee -a "$LOG_DIR/temp_log.txt"
echo "" | tee -a "$LOG_DIR/temp_log.txt"

# Prepare final log with timestamp in the name
LOG_TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/download_log_${PREFIX}${START}_to_${END}_$LOG_TIMESTAMP.txt"
mv "$LOG_DIR/temp_log.txt" "$LOG_FILE"

# Initialise counters
SUCCESS_COUNT=0
FAIL_COUNT=0
SKIPPED_COUNT=0

# Read and process the selected line range
CURRENT=1
while IFS= read -r url; do
    if [ "$CURRENT" -ge "$START" ] && [ "$CURRENT" -le "$END" ]; then
        custom_name="${PREFIX}${CURRENT}.root"
        target_path="$ROOT_DIR/$custom_name"

        if [ ! -f "$target_path" ]; then
            echo "🔽 Descargando $custom_name desde $url..." | tee -a "$LOG_FILE"
            xrdcp "$url" "$target_path"
            if [ $? -eq 0 ]; then
                echo "✅ Guardado como: $target_path" | tee -a "$LOG_FILE"
                echo "$custom_name : $url" >> "$LOG_FILE"
                ((SUCCESS_COUNT++))
            else
                echo "❌ Error al descargar $custom_name (línea $CURRENT)" | tee -a "$LOG_FILE"
                ((FAIL_COUNT++))
            fi
        else
            echo "ℹ️ El archivo $custom_name ya existe. Se omite descarga." | tee -a "$LOG_FILE"
            ((SKIPPED_COUNT++))
        fi
    fi
    ((CURRENT++))
done < "$TXT_FILE"

# End timestamp
END_TIME=$(date "+%Y-%m-%d %H:%M:%S")

# Summary
{
    echo ""
    echo "📁 Archivos almacenados en: $ROOT_DIR"
    echo "🗒 Log guardado en: $LOG_FILE"
    echo ""
    echo "📊 Resumen de la operación:"
    echo "   ✅ Descargados correctamente : $SUCCESS_COUNT"
    echo "   ❌ Fallidos                 : $FAIL_COUNT"
    echo "   ⚠️  Omitidos (ya existen)    : $SKIPPED_COUNT"
    echo ""
    echo "🕓 Fin: $END_TIME"
} | tee -a "$LOG_FILE"
