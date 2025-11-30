#!/bin/bash

# 04_merge_jpsitrktrktrk_outputs.sh
#
# Purpose:
#   Single-pass merge of all *_jpsitrktrktrk.root files into one ROOT file
#   using ROOT 5.32 (CMSSW_5_3_32 environment).
#   The script:
#     - Scans the JPsiTrkTrkTrk output directory for *_jpsitrktrktrk.root files.
#     - Excludes any file that already contains 'merged' in its name.
#     - Writes a manifest with the full list of input files.
#     - Runs a single 'hadd -f' call to create the final merged ROOT file.
#     - Prints a short preview of the inputs and a quick size check of the output.
#
# Directories and files:
#   INPUT_DIR : location of the per-file JPsiTrkTrkTrk output ROOT files.
#   OUT_DIR   : directory where the merged file and manifest will be stored.
#   OUT_FILE  : path to the merged ROOT file.
#   MANIFEST  : text file listing all input ROOT files used in the merge.
#
# Notes:
#   - Compatible with older bash (no 'mapfile -d' is used).
#   - Assumes ROOT's 'hadd' is available in the current CMSSW_5_3_32 environment.

set -euo pipefail

# === CONFIGURATION ===
INPUT_DIR="/code/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor/output_files_jpsitrktrktrk"
OUT_DIR="${INPUT_DIR}/merged"
OUT_FILE="${OUT_DIR}/muonia_all_jpsitrktrktrk_merged.root"
MANIFEST="${OUT_DIR}/merged_input_list.txt"

mkdir -p "$OUT_DIR"

echo "🔎 Buscando archivos en: $INPUT_DIR"

# === COLLECT FILES TO MERGE ===
# Collect file names (without path), sort them "naturally",
# and exclude any file that already has 'merged' in its name.
# Use newline separators (no NUL), which is sufficient for normal names.
files=()
while IFS= read -r fname; do
  # Skip empty lines for safety
  [[ -z "$fname" ]] && continue
  files+=( "${INPUT_DIR}/${fname}" )
done < <(
  find "$INPUT_DIR" -maxdepth 1 -type f -name "*_jpsitrktrktrk.root" ! -name "*merged*.root" -printf "%f\n" \
  | sort -V
)

if ((${#files[@]} == 0)); then
  echo "❌ No se encontraron archivos *_jpsitrktrktrk.root en: $INPUT_DIR"
  exit 1
fi

# === WRITE FULL MANIFEST ===
: > "$MANIFEST"
for f in "${files[@]}"; do
  printf "%s\n" "$f" >> "$MANIFEST"
done

# === PREVIEW INPUT LIST ===
echo "📦 Archivos a unir: ${#files[@]}"
# Preview (first 15)
preview_count=$(( ${#files[@]} < 15 ? ${#files[@]} : 15 ))
for ((i=0; i<preview_count; i++)); do
  printf '  • %s\n' "${files[i]}"
done
(( ${#files[@]} > 15 )) && echo "  … (y $((${#files[@]}-15)) más)"
echo "🗒️  Listado completo guardado en: $MANIFEST"

# Remove any previous incomplete merged file
[[ -f "$OUT_FILE" ]] && rm -f "$OUT_FILE"

# === RUN ROOT HADD ===
echo "🔧 Ejecutando un único hadd (ROOT 5.32, usando -f)…"
hadd -f "$OUT_FILE" "${files[@]}"

# === QUICK CHECK ===
if [[ -f "$OUT_FILE" ]]; then
  echo "✅ Merge completo:"
  echo "   → $OUT_FILE"
  if command -v du >/dev/null 2>&1; then
    echo -n "   Tamaño del archivo final: "
    du -h "$OUT_FILE" | awk '{print $1}'
  fi
  echo "   Entradas combinadas: ${#files[@]} (ver $MANIFEST)"
else
  echo "❌ Algo salió mal: no se generó $OUT_FILE"
  exit 1
fi
