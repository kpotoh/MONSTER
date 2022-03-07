#!/bin/bash

DOCSTRING="
USAGE: sh ./run_pipeline.sh [OPTION...] [FILE...]\n
\n
Positiosal arguments:\n
FILE - one or many tsv files with a few records (5-10)\n
\n
Options:\n
-t, --threads [INT]\t         number of threads (default: 0 (all))\n
-o, --outdir [STR] \t         output directory (default: ./data)\n
-p, --tmpdir [STR] \t         temporary directory (default: ./tmp)\n
-s, --show-browser \t         indicate if the browser should be visible\n
-c, --context      \t\t       indicate if context of mutations should be taken into account when calculating mutational spectrum\n
-r, --testrun      \t\t       first test run (required on new machine)\n
-d, --dry-run      \t\t       print the job to run on stdout (need for debugging)\n
-v, --verbose      \t\t       print all input filenames\n
-h, --help         \t\t       print this help\n
"
VERBOSE=NO
TESTRUN=NO
CONTEXT=""
SHOW_BROWSER=""
DRY=""
THREADS=0  # which will run as many jobs in parallel as possible
PYTHON=/usr/bin/python3  # change here to run locally
PARALLEL=/usr/bin/parallel
SCRIPT=/opt/MONSTER/scripts/pipeline/MONSTER.py
DEFAULT_DATA=/opt/MONSTER/data/share/ControlFile1.txt
OUTDIR=data
TMPDIR=tmp
POSITIONAL_ARGS=()

hit_size_treshold=50
hits_nb=500


if [ $# -eq 0 ]; then
  echo -e "Error! No arguments supplied\n"
  echo -e $DOCSTRING
  exit 0
fi

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      echo -e $DOCSTRING
      exit 0
      ;;
    -v|--verbose)
      VERBOSE=YES
      shift # past argument
      ;;
    -r|--testrun)
      TESTRUN=YES
      shift # past argument
      ;;
    -c|--context)
      CONTEXT="-c"
      shift # past argument
      ;;
    -s|--show-browser)
      SHOW_BROWSER="-v"
      shift # past argument
      ;;  
    -d|--dry-run)
      DRY="--dry-run"
      shift # past argument
      ;;  
    -t|--threads)
      THREADS="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--outdir)
      OUTDIR="$2"
      shift # past argument
      shift # past value
      ;;
    -p|--tmpdir)
      TMPDIR="$2"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done


if [[ $TESTRUN == YES ]]; then
  echo -e "Start test run\nWait 1-2 min\n"
  mkdir -p $TMPDIR
  $PARALLEL $DRY --jobs $THREADS $PYTHON $SCRIPT -b --input_file {} --out_folder $TMPDIR/output_{/.}-$(date --iso-8601=seconds) ::: $DEFAULT_DATA
  echo -e "\nTest run done. Test output in $TMPDIR/output_ControlFile1/ (directory could be empty)"
  exit 0
fi

DIR=$OUTDIR/run-$(date --iso-8601=seconds)

if [[ $DRY != "--dry-run" ]]; then
  mkdir -p $DIR
  echo -e "Created directory for current run:\n./$DIR\n"
fi

if [[ $VERBOSE == YES ]]; then
  IFS=$'\n'
  echo -e "Running MONSTER pipeline on files:\n${POSITIONAL_ARGS[*]}\n"
fi

echo "Run parallel computing"
$PARALLEL $DRY --jobs $THREADS $PYTHON $SCRIPT $CONTEXT $SHOW_BROWSER -b \
    --hit_size_treshold $hit_size_treshold --hits_nb $hits_nb --input_file {} \
    --out_folder $DIR/output_{/.} ::: ${POSITIONAL_ARGS[*]}

echo "Done"
