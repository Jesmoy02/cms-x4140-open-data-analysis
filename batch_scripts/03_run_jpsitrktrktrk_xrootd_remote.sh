#!/bin/bash

# 03_run_jpsitrktrktrk_xrootd_remote.sh
#
# Purpose:
#   Batch runner for the JPsiTrkTrkTrk CMSSW configuration using remote AOD files
#   accessed via XRootD.
#   The script:
#     - Reads XRootD URLs from index files in `muonia_links/*_file_index.txt`.
#     - Lets the user select which index file to use and which line range to process
#       (or uses the full range when SKIP_PROMPT=1).
#     - Runs cmsRun with `poet_cfg_jpsitrktrktrk.py` for each selected URL.
#     - Writes skimmed outputs to `output_files_jpsitrktrktrk/` with a
#       `_jpsitrktrktrk.root` suffix.
#     - Stores per-job logs in `output_files_jpsitrktrktrk/logs/` and prints an
#       OK/FAIL/SKIP summary, including a suggestion to retry failed jobs.
#
# Directories:
#   BASE_DIR  : base path of the PhysObjectExtractorTool area.
#   LINK_DIR  : location of the file-index .txt files with XRootD URLs.
#   OUTPUT_DIR: destination for the JPsiTrkTrkTrk skimmed output ROOT files.
#   LOG_DIR   : destination for per-job cmsRun logs.
#
# Notes:
#   - OVERWRITE=1 (default) overwrites existing outputs; set OVERWRITE=0 to keep them.
#   - SKIP_PROMPT=1 selects the first index file and processes the full line range.
#   - cmsRun signature: cmsRun poet_cfg_jpsitrktrktrk.py <input> <output>
#
# ===============================
# JPsiTrkTrkTrk ONLINE batch runner (local-style)
# - Reads XRootD links from muonia_links/*_file_index.txt
# - Colours, per-file timing, OK/FAIL/SKIP summary
# - OVERWRITE=1 by default (overwrite existing outputs)
# - SKIP_PROMPT=1 to process without prompts
# - Output suffix: *_jpsitrktrktrk_preselection.root
# - cmsRun call: cmsRun poet_cfg_jpsitrktrktrk.py <input> <output>
# ===============================

set -u  # do not use 'set -e' so the loop is not stopped on a single file failure

# ---------- Colours ----------
if command -v tput >/dev/null 2>&1; then
  t_red="$(tput setaf 1)"; t_green="$(tput setaf 2)"; t_yellow="$(tput setaf 3)"
  t_blue="$(tput setaf 4)"; t_cyan="$(tput setaf 6)"; t_bold="$(tput bold)"
  t_dim="$(tput dim)"; t_reset="$(tput sgr0)"
else
  t_red=""; t_green=""; t_yellow=""; t_blue=""; t_cyan=""
  t_bold=""; t_dim=""; t_reset=""
fi
hr() { printf '%s\n' "------------------------------------------------------------"; }

# ===============================
# 1) PATHS AND CONFIGURATION
# ===============================
BASE_DIR="/code/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor"
LINK_DIR="$BASE_DIR/muonia_links"

OUTPUT_DIR="$BASE_DIR/output_files_jpsitrktrktrk"
LOG_DIR="${OUTPUT_DIR}/logs"

CFG_JPSITR="$BASE_DIR/python/cmssw_cfg_jpsitrktrktrk.py"
TAG="jpsitrktrktrk"

OVERWRITE="${OVERWRITE:-1}"      # 1=sobrescribe, 0=no
SKIP_PROMPT="${SKIP_PROMPT:-0}"  # 1=saltar prompts

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo -e "${t_cyan}🌐 Iniciando análisis JPsiTrkTrkTrk ONLINE${t_reset}"
echo "📚 LINK_DIR  : $LINK_DIR"
echo "📤 OUTPUT_DIR: $OUTPUT_DIR"
echo "🧾 LOG_DIR   : $LOG_DIR"
echo "⚙️  CFG       : $CFG_JPSITR"
echo "🏷️  TAG       : $TAG"
echo "♻️  OVERWRITE : $OVERWRITE"
hr

# Minimal checks
command -v cmsRun >/dev/null 2>&1 || { echo -e "${t_red}❌ Falta cmsRun en PATH${t_reset}"; exit 1; }
[ -f "$CFG_JPSITR" ] || { echo -e "${t_red}❌ No existe cfg: $CFG_JPSITR${t_reset}"; exit 1; }
[ -d "$LINK_DIR" ] || { echo -e "${t_red}❌ No existe LINK_DIR: $LINK_DIR${t_reset}"; exit 1; }

# ===============================
# 2) AVAILABLE INDEX FILES
# ===============================
LINK_FILES=(
  "CMS_Run2011A_MuOnia_AOD_12Oct2013-v1_00000_file_index.txt"
  "CMS_Run2011A_MuOnia_AOD_12Oct2013-v1_00001_file_index.txt"
  "CMS_Run2011A_MuOnia_AOD_12Oct2013-v1_20000_file_index.txt"
  "CMS_Run2011A_MuOnia_AOD_12Oct2013-v1_210000_file_index.txt"
)

echo "📁 Archivos índice en $LINK_DIR:"
for i in "${!LINK_FILES[@]}"; do
  file="$LINK_DIR/${LINK_FILES[$i]}"
  if [ -f "$file" ]; then
    num_lines=$(wc -l < "$file")
    printf "  [%d] %s (%d líneas)\n" $((i+1)) "${LINK_FILES[$i]}" "$num_lines"
  else
    printf "  [%d] %s %s\n" $((i+1)) "${LINK_FILES[$i]}" "❌ (no encontrado)"
  fi
done
hr

# ===============================
# 3) INDEX AND RANGE SELECTION
# ===============================
if [[ "$SKIP_PROMPT" == "1" ]]; then
  FILE_CHOICE=1
  INDEX_PATH="$LINK_DIR/${LINK_FILES[$((FILE_CHOICE-1))]}"
  [ -f "$INDEX_PATH" ] || { echo -e "${t_red}❌ Falta índice: $INDEX_PATH${t_reset}"; exit 1; }
  total_lines=$(wc -l < "$INDEX_PATH")
  START=1; END=$total_lines
  echo -e "➡️  ${t_bold}Índice auto:${t_reset} ${LINK_FILES[$((FILE_CHOICE-1))]}"
  echo -e "➡️  ${t_bold}Rango auto:${t_reset} $START .. $END  (total líneas: $total_lines)"
else
  read -p "➡️  Selecciona un índice (1-${#LINK_FILES[@]}): " FILE_CHOICE
  [[ "$FILE_CHOICE" =~ ^[1-4]$ ]] || { echo -e "${t_red}❌ Selección inválida.${t_reset}"; exit 1; }
  INDEX_PATH="$LINK_DIR/${LINK_FILES[$((FILE_CHOICE-1))]}"
  [ -f "$INDEX_PATH" ] || { echo -e "${t_red}❌ No existe: $INDEX_PATH${t_reset}"; exit 1; }

  total_lines=$(wc -l < "$INDEX_PATH")
  read -p "➡️  Ingresa el número de línea de INICIO (ENTER = 1): " START
  read -p "➡️  Ingresa el número de línea de FIN (ENTER = ${total_lines}): " END
  START=${START:-1}; END=${END:-$total_lines}

  if ! [[ "$START" =~ ^[0-9]+$ && "$END" =~ ^[0-9]+$ ]]; then
    echo -e "${t_red}❌ Rango inválido: use números enteros.${t_reset}"; exit 1
  fi
  if [[ "$START" -lt 1 || "$END" -gt "$total_lines" || "$START" -gt "$END" ]]; then
    echo -e "${t_red}❌ Rango fuera de límites (1..${total_lines}) o invertido.${t_reset}"; exit 1
  fi
fi

PREFIX="muonia_${FILE_CHOICE}"

# ===============================
# 4) READ LINKS FROM RANGE
# ===============================
mapfile -s $((START - 1)) -n $((END - START + 1)) -t LINKS < "$INDEX_PATH"
if [ "${#LINKS[@]}" -eq 0 ]; then
  echo -e "${t_yellow}⚠️  No hay enlaces en el rango especificado.${t_reset}"
  exit 1
fi

sel_total=${#LINKS[@]}
echo -e "🧮 Se procesarán ${t_bold}$sel_total${t_reset} enlaces."
hr

# ===============================
# 5) PROCESSING
# ===============================
ok=0; fail=0; skipped=0
ok_list=(); fail_list=(); skipped_list=()

sec_to_hms() {
  local s=$1 h m; ((h=s/3600)); ((m=(s%3600)/60)); ((s=s%60))
  printf "%02d:%02d:%02d" "$h" "$m" "$s"
}

idx=0
for link in "${LINKS[@]}"; do
  idx=$((idx+1))
  link="$(echo "$link" | tr -d '\r' | xargs)"
  if [[ -z "$link" ]]; then
    echo -e "${t_yellow}⏭️  Línea vacía, se omite.${t_reset}"
    skipped=$((skipped+1)); skipped_list+=("linea_${idx}_vacia"); hr; continue
  fi

  # Use the absolute line number in the index file
  line_abs=$((START + idx - 1))

  # Name "muonia" based on the actual line number
  muonia_name="${PREFIX}_${line_abs}"

  # Base name for output and log
  base="${muonia_name}"
  # (To keep output with sub-range index instead, you could use: base="${PREFIX}_${idx}")

  output_file="${OUTPUT_DIR}/${base}_${TAG}.root"
  log_file="${LOG_DIR}/${base}.log"

  # Main message (with muonia_X_Y)
  echo -e "${t_blue}[$idx/$sel_total]${t_reset} 🔧 Procesando ${t_bold}$link${t_reset} (${muonia_name})"

  if [[ -f "$output_file" && "$OVERWRITE" != "1" ]]; then
    echo "⏭️  Saltando (ya existe): $(basename "$output_file")"
    skipped=$((skipped+1)); skipped_list+=("$base"); hr; continue
  elif [[ -f "$output_file" && "$OVERWRITE" == "1" ]]; then
    echo "♻️  Sobrescribiendo → $(basename "$output_file")"
  fi

  t0=$(date +%s)

  # cmsRun with post-processing: append " (muonia_X_Y)" ONLY to the line "[CFG] Archivo de entrada:"
  cmsRun "$CFG_JPSITR" "$link" "$output_file" 2>&1 \
  | awk -v tag=" (${muonia_name})" '
      BEGIN { pat = "\\[CFG\\][[:space:]]+Archivo de entrada:" }
      { gsub(/\r/, "") }  # clean CR in case they come from eos/xrd
      $0 ~ pat {
        if (index($0, tag) == 0) { print $0 tag; fflush(); next }
      }
      { print; fflush() }
    ' \
  | tee "$log_file"
  code=${PIPESTATUS[0]}

  t1=$(date +%s); dt=$((t1 - t0)); t_hms=$(sec_to_hms "$dt")

  if [ $code -eq 0 ]; then
    if [ ! -s "$output_file" ] || [ "$(stat -c%s "$output_file" 2>/dev/null || echo 0)" -lt 10240 ]; then
      echo -e "⚠️  ${t_yellow}Salida muy pequeña, puede estar vacía.${t_reset} (${t_dim}$t_hms${t_reset})"
      fail=$((fail+1)); fail_list+=("$base (archivo pequeño)")
    else
      echo -e "✅ ${t_green}OK${t_reset} (${t_dim}$t_hms${t_reset})  →  ${t_dim}$(basename "$output_file")${t_reset}"
      ok=$((ok+1)); ok_list+=("$base")
    fi
  else
    echo -e "❌ ${t_red}ERROR=$code${t_reset} (${t_dim}$t_hms${t_reset}) — revisa log: ${t_dim}$log_file${t_reset}"
    echo -e "${t_yellow}─── Últimas líneas del log ───${t_reset}"
    tail -n 8 "$log_file" || true
    fail=$((fail+1)); fail_list+=("$base")
  fi
  hr
done

# ===============================
# 6) FINAL SUMMARY
# ===============================
echo -e "${t_bold}🏁 Análisis ONLINE completado.${t_reset}"
echo -e "   ✔️  Éxitos : ${t_green}$ok${t_reset}"
echo -e "   ⏭️  Skips  : ${t_yellow}$skipped${t_reset}"
echo -e "   ❌ Fallos  : ${t_red}$fail${t_reset}"
echo -e "📂 Resultados en: ${t_cyan}$OUTPUT_DIR${t_reset}"
echo -e "🧾 Logs en     : ${t_cyan}$LOG_DIR${t_reset}"

if [[ $ok -gt 0 ]]; then
  echo -e "\n${t_green}✔️  OK:${t_reset}"; for f in "${ok_list[@]}"; do echo "  - $f"; done
fi
if [[ $skipped -gt 0 ]]; then
  echo -e "\n${t_yellow}⏭️  Skips:${t_reset}"; for f in "${skipped_list[@]}"; do echo "  - $f"; done
fi
if [[ $fail -gt 0 ]]; then
  echo -e "\n${t_red}❌ Fallidos:${t_reset}"; for f in "${fail_list[@]}"; do echo "  - $f"; done
  echo -e "\n${t_yellow}Sugerencia para reintentar fallidos:${t_reset}"
  echo "for base in \\"
  for f in "${fail_list[@]}"; do echo "  \"$f\" \\"; done
  echo "; do"
  echo "  # cmsRun \"$CFG_JPSITR\" \"<link_remoto>\" \"$OUTPUT_DIR/\${base}_${TAG}.root\" 2>&1 | tee \"$LOG_DIR/\${base}.log\""
  echo "done"
fi
