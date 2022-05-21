#!/bin/sh

echo "Checking #include ordering and spacing"
$1/scripts/code/BclCleaner.py $1 $2
echo "Done checking #include ordering and spacing"

echo "Checking example / class comments and their linkages"
$1/scripts/code/SummarizeExamples.py $1 --screen --errors
echo "Done checking example / class comments and their linkages" 
