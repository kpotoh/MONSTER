#!/bin/bash
# USAGE: sh ./run_pipeline.sh [FILE...]

DOCSTRING="
USAGE: sh ./run_pipeline.sh [OPTION...] [FILE...]\n
\n
Positiosal arguments:\n
FILE - one or many tsv files with a few records (5-10)\n
\n
Options:\n
-p, --python-path [PATH]    path to python executable\n
-t, --threads [INT]         number of threads (default: 0 (all))\n
-h, --help\n
-v, --verbose\n
"
VERBOSE=NO
THREADS=0  # which will run as many jobs in parallel as possible.
PYTHON=/usr/bin/python3
PARALLEL=/usr/bin/parallel
SCRIPT=scripts/pipeline/MONSTER.py
DIR=data/run_$(date --iso-8601=seconds)
POSITIONAL_ARGS=()

if [ $# -eq 0 ]; then
  echo -e "Error! No arguments supplied\n"
  echo -e $DOCSTRING
  exit 0
fi

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      echo $DOCSTRING
      exit 0
      ;;
    -v|--verbose)
      VERBOSE=YES
      shift # past argument
      ;;
    -p|--python-path)
      PYTHON="$2"
      shift # past argument
      shift # past value
      ;;
    -t|--threads)
      THREADS="$2"
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

mkdir -p $DIR
echo -e "Created directory for current run:\n./$DIR\n"
if [[ $VERBOSE -eq YES ]]; then
  IFS=$'\n'
  echo -e "Running MONSTER pipeline on files:\n${POSITIONAL_ARGS[*]}\n"
fi

echo "Start parallel computing"
$PARALLEL --jobs $THREADS $PYTHON $SCRIPT -b --input_file {} --out_folder $DIR/output_{/.} ::: $POSITIONAL_ARGS
echo "Done"
