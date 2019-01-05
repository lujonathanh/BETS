#!/usr/bin/env bash


module load anaconda

echo "Reading from $scriptlist"

echo "" >> submission_log.txt
echo "******************************" >> submission_log.txt
echo "$(date)" >> submission_log.txt
echo "$(pwd)" >> submission_log.txt
echo "Reading from $scriptlist" >> submission_log.txt

while read script; do
    echo Submitting Parallel Script $script
    clusterize -l $TIMEUPPERBOUND:00 -m $MEMUPPERBOUND -n $NNODES -p $PPN  -c "time $script & wait" >> submission_log.txt
done < $scriptlist