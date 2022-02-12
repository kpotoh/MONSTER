#!/bin/bash
# USAGE: sh ./run_pipeline.sh [FILE...]

# TODO add here condition on option -h
if [ $# -eq 0 ] 
  then
    echo "Error! No arguments supplied"
    echo "USAGE: sh ./run_pipeline.sh [FILE...]"
    exit 0
fi

PYTHON=/usr/bin/python3
PARALLEL=/usr/bin/parallel
SCRIPT=scripts/pipeline/MONSTER.py
DIR=data/run_$(date --iso-8601=seconds)

mkdir -p $DIR
echo "Created directory for current run:\n$DIR\n"

echo "Start parallel computing"
$PARALLEL $PYTHON $SCRIPT -b --input_file {} --out_folder $DIR/output_{/.} ::: $@
echo "Done"
