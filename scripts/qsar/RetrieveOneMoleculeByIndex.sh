#!/bin/sh

bcl.exe molecule:Properties -input_filenames $1 -input_start $2 -input_max 1 -output $3 > /dev/null
