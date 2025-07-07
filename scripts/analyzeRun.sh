#!/bin/bash

REPO_HOME=$(git rev-parse --show-toplevel) # does not have / at the end. 

# Set the default input and output directories, to be overwritten by the settings
inputDirectory="" 
outputDirectory=$REPO_HOME"/output/"

# Initialize variables
settingsFile=""
fileName=""
noCompile=false
cleanCompile=false

# Function to print help message
print_help() {
    echo "*****************************************************************************"
    echo "Usage: $0 -f <first_run> [-l <last_run>] -j <settings_file> [-h]"
    echo "  -p | --home-path      path to the code home directory (/your/path/tof-reco, no / at the end). There is no default."
    echo "  -r | --run-name       full run name, with ot without .dat"
    echo "  -f | --first-run      first run number"
    echo "  -l | --last-run       last run number" 
    echo "  -j | --json-settings  json settings file to run this app"
    echo "  --no-compile          do not recompile the code. Default is to recompile the code"
    echo "  --clean-compile       clean and recompile the code. Default is to recompile the code without cleaning"
    echo "  --source-lxplus       source the lxplus environment. Default is to not source the lxplus environment"
    echo "  -h | --help           print this help message"
    echo "*****************************************************************************"
    exit 0
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
      -p|--home-path)
        REPO_HOME="$2"
        shift 2
        ;;
      -r|--run-name)
        fileName="$2"
        shift 2
        ;;
      -f|--first-run)
        firstRun="$2"
        shift 2
        ;;
      -l|--last-run)
        lastRun="$2"
        shift 2
        ;;
      -j|--json-settings)
        settingsFile="$2"
        shift 2
        ;;
      --no-compile)
        noCompile=true
        shift
        ;;
      --clean-compile)
        cleanCompile=true
        shift
        ;;
      --source-lxplus)
        sourceLxplus=true
        shift
        ;;
      -h|--help)
        print_help
        ;;
      *)
        shift
        ;;
    esac
done

if [ -z "$settingsFile" ]
then
  echo "Please specify the settings file, using flag -j <settings_file>. Stopping execution."
  exit 0
fi

# THIS MIGHT BECOME NECESSARY LATER
# if sourceLxplus is true, source the lxplus environment.
# if [ "$sourceLxplus" = true ]; then  
#   echo "Sourcing the lxplus environment"
#   source $REPO_HOME/scripts/source-lxplus.sh
# fi

echo "Home is ${REPO_HOME}."
echo "Currently in $(pwd), moving to ${REPO_HOME} to run the script." # might not be necessary if using always absolute paths
cd $REPO_HOME

# read from settings file the input and output paths
inputDirectory=$(awk -F'"' '/inputDirectory/{print $4}' "$settingsFile")
outputDirectory=$(awk -F'"' '/outputDirectory/{print $4}' "$settingsFile")
verbose=$(awk -F'"' '/verboseMode/{print $4}' "$settingsFile")
debug=$(awk -F'"' '/debugMode/{print $4}' "$settingsFile")
nsigma=$(awk -F'"' '/nSigma/{print $4}' "$settingsFile") # TODO add fallback if not defined
showPlots=$(awk -F'"' '/showPlots/{print $4}' "$settingsFile") # TODO add fallback if not defined

# Check if the user has selected a run name or number(s)
if [ -z "$fileName" ] 
then
    if [ -z "$firstRun" ]
    then
    echo "Please select a run number (or an interval) from command line to convert the data into ROOT format. 
    You can select only the first run, or both the first and the last run.
    Usage: 
    $ ./analyzeRuns.sh -f <first_run> [-l <last_run>] -j <json_settings_file> [-h]"
    exit 0
    fi
    if [ -z "$lastRun" ]
    then
    echo "I have received only the number of the first run. Reading only run $firstRun"
    lastRun=$firstRun
    fi
else
    echo "I have received the run name: $fileName"
    # if it has the extension, remove it
    fileName=$(echo $fileName | awk -F'.' '{print $1}')
    # firstRun=$(echo $fileName | awk -F'_' '{print $2}')
fi

################################################
# Compile the code, with checks
if [ "$noCompile" = false ]
then
  echo "Currently we are in"
  pwd
  mkdir -p $REPO_HOME/build
  echo "Moving into build directory"
  cd $REPO_HOME/build;

  if [ "$cleanCompile" = true ]
  then
    echo "Cleaning the build directory"
    rm -r $REPO_HOME/build/*
  fi  
  # Compile the code
  . ${REPO_HOME}/scripts/compile.sh
else 
  echo "Skipping compilation."
fi

################################################
# Execute the code
cd $REPO_HOME/build # executable is here, for now

if [ -n "$fileName" ]
then
  # just to get in the loop below
  firstRun=0
  lastRun=0
fi

# Iterate over all the selected runs
for ((runit = $firstRun; runit <= $lastRun; runit++ ))
do
  if [ -z $fileName ]
  then
    echo "Looking for files by run number"
    echo "runit: $runit"
    runit_5d=$(printf "%05d" "$runit")
    filePath=$(find "$inputDirectory" -name "SCD_RUN${runit_5d}_*.dat")
    
    # Stop execution if the selected run is not present in the input directory
    if [ -z "$filePath" ]
    then
        echo "Run "$runit " not found. Stopping execution."
        continue
    fi

    echo "For this run number, found the following file(s):"
    echo $filePath

    echo "Dealing with file: "$filePath
    echo "outputDirectory: "$outputDirectory
  

    # in case it doesn't exist, creating output directory
    mkdir -p $outputDirectory

    echo ""
    echo "Currently we are in"
    pwd
    echo ""

    fileName=$(basename $filePath)
    fileName=${fileName%.*}
    echo "File name is: "$fileName
  else
    echo "Received a single file path ${fileName}"
    filePath="${inputDirectory}/${fileName}.dat"
  fi

  if [ -f "${outputDirectory}/${fileName}.root" ]
  then
      echo "File ${outputDirectory}/${fileName}.root already exists. Skipping conversion."
  else
      convert_data="./PAPERO_convert ${filePath} ${outputDirectory}/${fileName}.root --dune"
      echo "Executing command: "$convert_data
      $convert_data
  fi

  if [ -f "${outputDirectory}/${fileName}.cal" ]
  then
      echo "File ${outputDirectory}/${fileName}.cal already exists. Skipping calibration extraction."
  else
      extract_calibration="./calibration ${outputDirectory}/${fileName}.root --output ${outputDirectory}/${fileName} --dune --fast"
      echo "Executing command: "$extract_calibration
      $extract_calibration
  fi

  analyze_data="./dataAnalyzer -r ${outputDirectory}/${fileName}.root -c ${outputDirectory}/${fileName}.cal -o ${outputDirectory} -s ${nsigma}"
  
  if [ "$verbose" = true ]
  then
      echo "Verbose mode is on."
      analyze_data+=" -v"
  fi

  if [ "$debug" = true ]
  then
      echo "Debug mode is on."
      analyze_data+=" -d"
  fi

  if [ "$showPlots" = true ]
  then
      analyze_data+=" --show-plots"
  fi

  echo "Executing command: "$analyze_data
  $analyze_data
 
done

echo "All runs have been analyzed. Exiting."