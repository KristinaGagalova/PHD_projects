#!/usr/bin/bash

input_file=$1

while read -r f1 c1; do
  # split the comma delimited field 'c1' into its constituents
  for c in ${c1//|/ }; do
     printf "$f1\t$c\n"
  done
done < $input_file
