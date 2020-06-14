#!/usr/bin/env bash

echo "Reading from $scriptlist"
while read script; do
    echo Submitting Parallel Script $script
    sh -c "$script"
done < $scriptlist
echo "ALL JOBS FROM $scriptlist COMPLETED. Feel free to continue"
