#!/bin/bash

set -e

export SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPTS_DIR/init.sh"

# Defaults
HOME_DIR=${SCRIPTS_DIR%/scripts}
settingsFile="${HOME_DIR}/json/ev-settings.json"
noCompile=false
cleanCompile=false
runNumber=""
smin=2
smax=15
sstep=1

print_help() {
  echo "*****************************************************************************"
  echo "Usage: $0 -r <run_number> [-j <settings_file>] [--smin N --smax N --sstep N] [--no-compile|--clean-compile] [-h]"
  echo "  -p | --home-path      path to repo home (default: autodetected)"
  echo "  -r | --run-number     run number (required)"
  echo "  -j | --json-settings  settings JSON (default: json/ev-settings.json)"
  echo "       --smin           min sigma (default: 2)"
  echo "       --smax           max sigma (default: 15)"
  echo "       --sstep          sigma step (default: 1)"
  echo "       --no-compile     skip compilation"
  echo "       --clean-compile  clean and compile"
  echo "  -h | --help           print this help"
  echo "*****************************************************************************"
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--home-path)      HOME_DIR="$2"; shift 2 ;;
    -r|--run-number)     runNumber="$2"; shift 2 ;;
    -j|--json-settings)  settingsFile="$2"; shift 2 ;;
    --smin)              smin="$2"; shift 2 ;;
    --smax)              smax="$2"; shift 2 ;;
    --sstep)             sstep="$2"; shift 2 ;;
    --no-compile)        noCompile=true; shift ;;
    --clean-compile)     cleanCompile=true; shift ;;
    -h|--help)           print_help; exit 0 ;;
    *)                   shift ;;
  esac
done

if [ -z "$runNumber" ]; then echo "Please specify -r <run_number>"; echo ""; print_help; exit 1; fi

inputDirectory=$(awk -F'"' '/inputDirectory/{print $4}' "$settingsFile")
outputDirectory=$(awk -F'"' '/outputDirectory/{print $4}' "$settingsFile")

echo "Using settings: $settingsFile"
echo "Input:  $inputDirectory"
echo "Output: $outputDirectory"

mkdir -p "$outputDirectory"

# Build as needed
BUILD_DIR="${HOME_DIR}/build"
if [ "$noCompile" = false ]; then
  echo "Preparing build in $BUILD_DIR"
  current_dir=$(pwd)
  mkdir -p "$BUILD_DIR"; cd "$BUILD_DIR"
  if [ "$cleanCompile" = true ]; then rm -rf "$BUILD_DIR"/*; fi
  . "$SCRIPTS_DIR/compile.sh"
else
  cd "$BUILD_DIR"
fi

# Locate .dat
filePath=$($SCRIPTS_DIR/findRun.sh -r "$runNumber" -i "$inputDirectory" | tail -n 1)
if [ -z "$filePath" ]; then echo "Run $runNumber not found"; exit 1; fi

fileName=$(basename "$filePath"); fileName=${fileName%.*}

# Convert if needed
if [ ! -f "${outputDirectory}/${fileName}.root" ]; then
  echo "Converting $filePath"
  ./PAPERO_convert "$filePath" "${outputDirectory}/${fileName}.root" --dune
fi

# Resolve calibration file: strictly use previous CAL run
resolve_prev_cal() {
  local runit=$1
  local calPath=""
  local prev=$((runit-1))
  while [ $prev -ge 0 ]; do
    local prevPaths=$($SCRIPTS_DIR/findRun.sh -r "$prev" -i "$inputDirectory" 2>/dev/null | tail -n 1)
    if [ -n "$prevPaths" ]; then
      for p in $prevPaths; do
        if basename "$p" | grep -qi "CAL"; then calPath="$p"; break; fi
      done
    fi
    [ -n "$calPath" ] && break
    prev=$((prev-1))
  done
  echo "$calPath"
}

calRunPath=$(resolve_prev_cal "$runNumber")
if [ -z "$calRunPath" ]; then echo "Previous CAL run not found for run $runNumber"; exit 1; fi

calBase=$(basename "$calRunPath"); calBase=${calBase%.*}
if [ ! -f "${outputDirectory}/${calBase}.root" ]; then
  echo "Converting CAL run: $calRunPath"
  ./PAPERO_convert "$calRunPath" "${outputDirectory}/${calBase}.root" --dune
fi
if [ ! -f "${outputDirectory}/${calBase}.cal" ]; then
  echo "Extracting CAL: ${calBase}.root"
  ./calibration "${outputDirectory}/${calBase}.root" --output "${outputDirectory}/${calBase}" --dune --fast
fi
calFile="${outputDirectory}/${calBase}.cal"

outPdf="${outputDirectory}/hits_vs_sigma_run${runNumber}.pdf"
echo "Running hits_vs_sigma..."
./hits_vs_sigma -r "${outputDirectory}/${fileName}.root" -c "$calFile" -o "$outPdf" -n "$runNumber" --smin "$smin" --smax "$smax" --sstep "$sstep"

# go back to previous location
cd $current_dir

echo "Done. PDF: $outPdf"
