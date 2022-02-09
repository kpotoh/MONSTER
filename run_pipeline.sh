#!/bin/bash
# USAGE: sh ./run_pipeline.sh [FILE...]

if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "USAGE: sh ./run_pipeline.sh [FILE...]"
    exit 0
fi

SCRIPT=scripts/pipeline/MONSTER.py
DIR=data/run_$(date --iso-8601=seconds)

mkdir $DIR
echo "Created directory for current run:\n$DIR\n"

echo "Start parallel computing"
parallel python3 $SCRIPT -b --input_file {} --out_folder $DIR/output_{/.} ::: $@
echo "Done"
