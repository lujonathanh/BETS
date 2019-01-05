#!/usr/bin/env bash

mkdir output-logs

export PAUSETIME=60

echo Move the output logs every $PAUSETIME seconds
while true; do mv *.out output-logs; sleep $PAUSETIME; done