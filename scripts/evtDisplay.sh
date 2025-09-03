#!/bin/bash

export SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPTS_DIR/init.sh

# Defaults
settingsFile=""
fileName=""
noCompile=false
cleanCompile=false

print_help() {
  echo "*****************************************************************************"
  echo "Usage: $0 -r <run_name|path> -j <settings_file> [-s <nSigma>] [--no-compile] [--clean-compile] [-h]"
  echo "  -r | --run-name       full run name without extension or absolute path to .dat/.root"
  echo "  -s | --n-sigma        number of sigma for threshold (overrides settings)"
  echo "  -j | --json-settings  json settings file"
  echo "  --no-compile          skip building"
  echo "  --clean-compile       clean build folder before building"
  echo "  -h | --help           print this help"
  echo "*****************************************************************************"
  exit 0
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--run-name)      fileName="$2"; shift 2 ;;
    -s|--n-sigma)       nsigma="$2"; shift 2 ;;
    -j|--json-settings) settingsFile="$2"; shift 2 ;;
    --no-compile)       noCompile=true; shift ;;
    --clean-compile)    cleanCompile=true; shift ;;
    -h|--help)          print_help ;;
    *) shift ;;
  esac
done

if [ -z "$settingsFile" ]; then
  echo "Please specify -j <settings_file>"; exit 1;
fi

cd "$HOME_DIR"

# Read defaults from settings
inputDirectory=$(awk -F'"' '/inputDirectory/{print $4}' "$settingsFile")
outputDirectory=$(awk -F'"' '/outputDirectory/{print $4}' "$settingsFile")
verbose=$(awk -F'"' '/verboseMode/{print $4}' "$settingsFile")
debug=$(awk -F'"' '/debugMode/{print $4}' "$settingsFile")
if [ -z "$nsigma" ]; then
  nsigma=$(awk -F'"' '/nSigma/{print $4}' "$settingsFile")
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

if [ "$noCompile" = false ]; then
  if [ "$cleanCompile" = true ]; then
    rm -rf "$BUILD_DIR"/*
  fi
  . ${HOME_DIR}/scripts/compile.sh
fi

cd "$BUILD_DIR"

# Resolve file base name and paths
if [[ -z "$fileName" ]]; then
  echo "Please provide -r <run_name|path>"; exit 1;
fi

if [[ "$fileName" == *.root ]]; then
  rootFile="$fileName"
  base=$(basename "$rootFile" .root)
elif [[ "$fileName" == *.dat ]]; then
  datPath="$fileName"
  base=$(basename "$datPath" .dat)
  rootFile="${outputDirectory}/${base}.root"
else
  # treat as run base name
  # Use finder to locate .dat, then derive base
  datPath=$($SCRIPTS_DIR/findRun.sh -r "$fileName" -i "$inputDirectory" | tail -n 1)
  if [ -z "$datPath" ]; then echo "Cannot locate run $fileName under $inputDirectory"; exit 1; fi
  base=$(basename "$datPath" .dat)
  rootFile="${outputDirectory}/${base}.root"
fi

# If needed, convert and extract calibration
if [ ! -f "$rootFile" ]; then
  # If .dat present in input dir, convert
  datCandidate=$(find "$inputDirectory" -name "${base}.dat" | head -n1)
  if [ -z "$datCandidate" ]; then
    echo "Cannot find ${base}.dat to convert, and ${rootFile} doesn't exist."; exit 1;
  fi
  ./PAPERO_convert "$datCandidate" "$rootFile" --dune || exit 1
fi

calFile="${outputDirectory}/${base}.cal"
if [ ! -f "$calFile" ]; then
  ./calibration "$rootFile" --output "${outputDirectory}/${base}" --dune --fast || exit 1
fi

# Launch interactive event display
cmd=("./event_display" -r "$rootFile" -c "$calFile" -s "$nsigma")
if [ -n "$RUN_NUMBER" ]; then
  cmd+=( -n "$RUN_NUMBER" )
fi
if [ "$verbose" = true ]; then
  cmd+=( -v )
fi
if [ "$debug" = true ]; then
  cmd+=( -d )
fi

echo "Executing: ${cmd[@]}"
"${cmd[@]}"
