#!/usr/bin/env bash

echo Submitting Parallel Script $script
clusterize -l $TIMEUPPERBOUND:00 -m $MEMUPPERBOUND -n $NNODES -p $PPN  -c "time $script & wait"
