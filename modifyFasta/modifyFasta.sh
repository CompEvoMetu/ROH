#!/bin/bash

file=$1
filebase=$(basename $file .bed)

python3 get_all.py ${file} > ${filebase}.new.bed
python3 get_seq.py ${file} ${filebase}
