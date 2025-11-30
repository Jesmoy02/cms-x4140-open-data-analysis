#!/bin/bash

# 02_run_jpsitrktrktrk_local_files.sh
#
# Purpose:
#   Batch runner for the JPsiTrkTrkTrk CMSSW configuration on locally stored ROOT files.
#   The script:
#     - Scans the directory `root_files/` and lists all available input ROOT files.
#     - Lets the user select a range of files to process (or processes all of them when SKIP_PROMPT=1).
#     - Runs cmsRun with `poet_cfg_jpsitrktrktrk.py` for each selected input file.
#     - Writes the skimmed outputs to `output_files_jpsitrktrktrk/` with a `_jpsitrktrktrk.root` suffix.
#     - Stores a per-file log in `output_files_jpsitrktrktrk/logs/` and prints a final OK/FAIL/SKIP summary.
#
# Directories:
#   INPUT_DIR  : location of the input ROOT files to be analysed.
#   OUTPUT_DIR : destination for the JPsiTrkTrkTrk skimmed output ROOT files.
#   LOG_DIR    : destination for per-file cmsRun logs.
#
# Notes:
#   - OVERWRITE=1 (default) overwrites existing outputs; set OVERWRITE=0 to keep them.
#   - SKIP_PROMPT=1 processes all files without asking for a range.

# ===============================
# JPsiTrkTrkTrk batch runner (improved)
# - Bash-compatible with SLC6
# - Colours, per-file timing, OK/FAIL/SKIP summary
# - OVERWRITE=1 by default (overwrite existing outputs)
# - SKIP_PROMPT=1 to process everything without prompts
# - Output suffix: *_jpsitrktrktrk_preselection.root
# ===============================

# We do not use 'set -e' so that the loop is not stopped if cmsRun fails
set -u

# ---------- Colours (if tput is available) ----------
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
INPUT_DIR="/code/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor/root_files"
OUTPUT_DIR="/code/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor/output_files_jpsitrktrktrk"
LOG_DIR="${OUTPUT_DIR}/logs"

# Current cfg:
CFG_JPSITR="/code/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor/python/cmssw_cfg_jpsitrktrktrk.py"

# Name tag for preselection outputs
TAG="jpsitrktrktrk"

# Default overwrite behaviour (1 = overwrite, 0 = do not overwrite)
OVERWRITE="${OVERWRITE:-1}"
# Process everything without prompts (1 = skip range prompts)
SKIP_PROMPT="${SKIP_PROMPT:-0}"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo -e "${t_cyan}🚀 Iniciando análisis JPsiTrkTrkTrk${t_reset}"
echo "📥 INPUT_DIR : $INPUT_DIR"
echo "📤 OUTPUT_DIR: $OUTPUT_DIR"
echo "🧾 LOG_DIR   : $LOG_DIR"
echo "⚙️  CFG       : $CFG_JPSITR"
echo "🏷️  TAG       : $TAG"
echo "♻️  OVERWRITE : $OVERWRITE"
hr

# Minimal checks
command -v cmsRun >/dev/null 2>&1 || { echo -e "${t_red}❌ Falta cmsRun en PATH${t_reset}"; exit 1; }
[ -f "$CFG_JPSITR" ] || { echo -e "${t_red}❌ No existe cfg: $CFG_JPSITR${t_reset}"; exit 1; }
[ -d "$INPUT_DIR" ] || { echo -e "${t_red}❌ No existe INPUT_DIR: $INPUT_DIR${t_reset}"; exit 1; }

# ===============================
# 2) DETECT AND SORT ROOT FILES
# ===============================
all_files=($(ls "$INPUT_DIR"/*.root 2>/dev/null | sort -V))
total_files=${#all_files[@]}

if [[ $total_files -eq 0 ]]; then
  echo -e "${t_red}❌ No se encontraron archivos .root en $INPUT_DIR${t_reset}"
  exit 1
fi

echo "📊 Archivos encontrados ($total_files):"
for i in "${!all_files[@]}"; do
  printf "  [%2d] %s\n" $((i+1)) "$(basename "${all_files[$i]}")"
done
hr

# ===============================
# 3) SELECT ANALYSIS RANGE
# ===============================
if [[ "$SKIP_PROMPT" == "1" ]]; then
  start=1
  end=$total_files
  echo -e "➡️  ${t_bold}Rango auto:${t_reset} $start .. $end  (total: $total_files)"
else
  echo ""
  read -p "➡️  Ingresa el número de archivo de INICIO (ENTER = 1): " start
  read -p "➡️  Ingresa el número de archivo de FIN (ENTER = $total_files): " end
  start=${start:-1}
  end=${end:-$total_files}
fi

if [[ $start -lt 1 || $end -lt $start || $end -gt $total_files ]]; then
  echo -e "${t_red}❌ Rango inválido. Cancelando.${t_reset}"
  exit 1
fi

selected_files=("${all_files[@]:$((start-1)):$((end-start+1))}")
sel_total=${#selected_files[@]}
echo -e "🧮 Se procesarán ${t_bold}$sel_total${t_reset} archivos."
hr

# ===============================
# 4) PROCESSING LOOP
# ===============================
ok=0
fail=0
skipped=0
ok_list=()
fail_list=()
skipped_list=()

sec_to_hms() {
  local s=$1 h m
  ((h=s/3600))
  ((m=(s%3600)/60))
  ((s=s%60))
  printf "%02d:%02d:%02d" "$h" "$m" "$s"
}

idx=0
for filepath in "${selected_files[@]}"; do
  idx=$((idx+1))
  filename=$(basename "$filepath")
  name="${filename%.*}"

  # Output name with preselection suffix
  output_file="${OUTPUT_DIR}/${name}_${TAG}.root"
  log_file="${LOG_DIR}/${name}.log"

  echo -e "${t_blue}[$idx/$sel_total]${t_reset} 🔧 Procesando ${t_bold}$filename${t_reset}"

  if [[ -f "$output_file" && "$OVERWRITE" != "1" ]]; then
    echo "⏭️  Saltando (ya existe): $(basename "$output_file")"
    skipped=$((skipped+1))
    skipped_list+=("$filename")
    hr
    continue
  elif [[ -f "$output_file" && "$OVERWRITE" == "1" ]]; then
    echo "♻️  Sobrescribiendo → $(basename "$output_file")"
  fi

  t0=$(date +%s)

  # Call consistent with your cfg: input.root output.root
  cmsRun "$CFG_JPSITR" "$filepath" "$output_file" 2>&1 | tee "$log_file"
  code=${PIPESTATUS[0]}

  t1=$(date +%s)
  dt=$((t1 - t0))
  t_hms=$(sec_to_hms "$dt")

  if [ $code -eq 0 ]; then
    echo -e "✅ ${t_green}OK${t_reset} (${t_dim}$t_hms${t_reset})  →  ${t_dim}$(basename "$output_file")${t_reset}"
    ok=$((ok+1))
    ok_list+=("$filename")
  else
    echo -e "❌ ${t_red}ERROR=$code${t_reset} (${t_dim}$t_hms${t_reset}) — revisa log: ${t_dim}$log_file${t_reset}"
    echo -e "${t_yellow}─── Últimas líneas del log ───${t_reset}"
    tail -n 8 "$log_file" || true
    fail=$((fail+1))
    fail_list+=("$filename")
  fi
  hr
done

# ===============================
# 5) FINAL SUMMARY
# ===============================
echo -e "${t_bold}🏁 Análisis completado.${t_reset}"
echo -e "   ✔️  Éxitos : ${t_green}$ok${t_reset}"
echo -e "   ⏭️  Skips  : ${t_yellow}$skipped${t_reset}"
echo -e "   ❌ Fallos  : ${t_red}$fail${t_reset}"
echo -e "📂 Resultados en: ${t_cyan}$OUTPUT_DIR${t_reset}"
echo -e "🧾 Logs en     : ${t_cyan}$LOG_DIR${t_reset}"

if [[ $ok -gt 0 ]]; then
  echo -e "\n${t_green}✔️  OK:${t_reset}"
  for f in "${ok_list[@]}"; do echo "  - $f"; done
fi

if [[ $skipped -gt 0 ]]; then
  echo -e "\n${t_yellow}⏭️  Skips:${t_reset}"
  for f in "${skipped_list[@]}"; do echo "  - $f"; done
fi

if [[ $fail -gt 0 ]]; then
  echo -e "\n${t_red}❌ Fallidos:${t_reset}"
  for f in "${fail_list[@]}"; do echo "  - $f"; done

  echo -e "\n${t_yellow}Sugerencia para reintentar fallidos:${t_reset}"
  echo "for f in \\"
  for f in "${fail_list[@]}"; do
    echo "  \"$INPUT_DIR/$f\" \\"
  done
  echo "; do"
  echo "  base=\$(basename \"\$f\" .root)"
  echo "  cmsRun \"$CFG_JPSITR\" \"\$f\" \"$OUTPUT_DIR/\${base}_${TAG}.root\" 2>&1 | tee \"$LOG_DIR/\${base}.log\""
  echo "done"
fi
